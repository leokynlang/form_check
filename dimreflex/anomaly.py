"""
dimreflex.anomaly
=================

后端异常检测 + 异常后的调整程序(可审计的守卫层)。

定位: 套在"log 空间线性拟合"(幂律等)外面的守卫。调用方通常传入
π 空间表(discover._prepare_pi 的产物), 也可以是任意正值列的表——
守卫只关心函数形式, 不关心量纲(量纲合法性已由上游 π 空间保证)。

三层检测:
1. 数据层  非有限值(NaN/inf)、非正值、零方差列 —— 致命, 直接显式失败
2. 残差层  MAD 稳健标准化残差(|z| > threshold 判离群), 对最小二乘
           被离群点带偏这件事本身不敏感
3. 几何层  高杠杆点(hat 矩阵对角 > 2p/n)、设计矩阵条件数(标准化后
           > 30 告警共线性风险)

调整策略(adjustment, 默认 'report_only'):
- 'report_only'  只检测与报告, 不动数据 —— 默认, 符合"不悄悄改数据"哲学
- 'robust'       Huber IRLS 稳健重拟合: 降权而非剔除, 数据一个不丢
- 'ridge'        k-fold 选择正则强度, 稳定共线特征的指数
- 'exclude'      剔除被标记点(离群 ∪ 高杠杆)后重拟合: 每个被剔除点
                 (索引/得分/来源)记录在案; 被标记超过 1/3 样本时拒绝自动
                 剔除并说明原因 —— 那多半是模型形式错了, 不是数据错了
- 'auto'         残差离群选 robust, 仅共线选 ridge, 无需干预则保持原拟合

返回的 before/after 候选、异常清单、审计记录全部显式 —— 调整永远可追溯。
"""


import numpy as np
import pandas as pd

from ._ols import ols as _ols
from ._ols import ols_r2 as _ols_r2
from .families import CandidateModel

_FATAL_TYPES = {
    'missing_column', 'non_numeric', 'non_finite',
    'non_positive', 'zero_variance',
}


# ---------------------------------------------------------------------------
# 第 1 层: 数据健康检查
# ---------------------------------------------------------------------------

def check_data_health(data: pd.DataFrame, target: str,
                      features: list[str]) -> list[dict]:
    """
    数据层体检: 缺列 / 非有限值 / 非正值 / 零方差列。
    返回问题清单, 每项 {'type', ..., 'severity'}; 致命项 severity='fatal'。
    清单为空 = 数据层健康。
    """
    issues: list[dict] = []
    cols = [target, *list(features)]

    missing = [c for c in cols if c not in data.columns]
    if missing:
        issues.append({'type': 'missing_column', 'columns': missing,
                       'severity': 'fatal'})
    present = [c for c in cols if c in data.columns]
    if not present:
        return issues

    values = {}
    non_numeric = []
    for column in present:
        try:
            values[column] = data[column].to_numpy(dtype=float)
        except (TypeError, ValueError):
            non_numeric.append(column)
    if non_numeric:
        issues.append({'type': 'non_numeric', 'columns': non_numeric,
                       'severity': 'fatal'})
    if not values:
        return issues

    non_finite = [c for c, v in values.items() if not np.isfinite(v).all()]
    if non_finite:
        issues.append({'type': 'non_finite', 'columns': non_finite,
                       'severity': 'fatal'})

    non_positive = [c for c, v in values.items() if (v <= 0).any()]
    if non_positive:
        issues.append({'type': 'non_positive', 'columns': non_positive,
                       'severity': 'fatal',
                       'note': 'log 空间拟合要求正值域'})

    zero_var = [c for c, v in values.items()
                if np.nanstd(v) < 1e-15 and np.isfinite(v).all()]
    if zero_var:
        issues.append({'type': 'zero_variance', 'columns': zero_var,
                       'severity': 'fatal',
                       'note': '常数列不携带信息, 且使杠杆值不可辨识'})

    return issues


# ---------------------------------------------------------------------------
# 第 2 层: 残差离群检测 (MAD 稳健 z 分数)
# ---------------------------------------------------------------------------

def detect_residual_outliers(residuals: np.ndarray,
                             threshold: float = 3.5) -> dict:
    """
    MAD 稳健 z 分数: z_i = (r_i − median) / (1.4826·MAD)。
    比 3σ 法稳健——离群点自身几乎不影响中位数与 MAD。
    MAD 退化(残差近乎常数)时显式说明, 不判任何离群。
    """
    r = np.asarray(residuals, dtype=float)
    med = float(np.median(r))
    scale = 1.4826 * float(np.median(np.abs(r - med)))
    if scale < 1e-12:
        return {'indices': [], 'scores': {}, 'threshold': threshold,
                'note': '残差 MAD 退化(近乎常数), 无法判离群'}
    z = (r - med) / scale
    idx = [int(i) for i in np.where(np.abs(z) > threshold)[0]]
    return {'indices': idx,
            'scores': {int(i): float(z[i]) for i in idx},
            'threshold': float(threshold)}


# ---------------------------------------------------------------------------
# 第 3 层: 杠杆点 + 共线性
# ---------------------------------------------------------------------------

def detect_leverage(X: np.ndarray, threshold: float | None = None) -> dict:
    """
    hat 矩阵对角 h_ii = x_i' (X'X)^{-1} x_i, 默认阈值 2p/n(经验规则)。
    高杠杆 = 该点在特征空间中位置极端, 会把拟合面拉向自己——
    与残差离群正交: 杠杆点可以残差很小却支配整条拟合。
    """
    X = np.asarray(X, dtype=float)
    n, p = X.shape
    if n <= p + 1:
        return {'indices': [], 'scores': {}, 'threshold': threshold,
                'note': f'样本 {n} ≤ 参数 {p}+1, 杠杆值不可辨识'}
    XtX_inv = np.linalg.pinv(X.T @ X)
    h = np.einsum('ij,jk,ik->i', X, XtX_inv, X)
    thr = float(threshold) if threshold is not None else 2.0 * p / n
    idx = [int(i) for i in np.where(h > thr)[0]]
    return {'indices': idx,
            'scores': {int(i): float(h[i]) for i in idx},
            'threshold': thr}


def detect_collinearity(X: np.ndarray, cond_threshold: float = 30.0) -> dict | None:
    """
    对特征列(不含常数列)逐列标准化后求条件数——
    标准化消除单位/量级影响, 只留下共线性本身。
    条件数 > 阈值(默认 30)告警: 指数估计不稳, Re/Pr 共享 μ 时常见。
    """
    X = np.asarray(X, dtype=float)
    if X.size == 0:
        return None
    sd = X.std(axis=0)
    sd[sd < 1e-15] = 1.0
    Z = (X - X.mean(axis=0)) / sd
    cond = float(np.linalg.cond(Z))
    if cond > cond_threshold:
        return {'type': 'collinearity', 'severity': 'warning',
                'condition_number': cond, 'threshold': cond_threshold,
                'note': '特征共线性风险: 指数估计可能不稳, 解释时需谨慎'}
    return None


# ---------------------------------------------------------------------------
# 调整程序 1: Huber IRLS 稳健重拟合(降权, 不剔除)
# ---------------------------------------------------------------------------

def robust_irls(X: np.ndarray, y: np.ndarray, k: float = 1.345,
                max_iter: int = 50, tol: float = 1e-10) -> dict:
    """
    Huber IRLS: 标准化残差 |u|≤k 的点权重 1, 其余权重 k/|u|。
    k=1.345 为 95% 效率的标准取值。返回系数、权重、迭代信息。
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    converged = False
    n_iter = 0
    w = np.ones(len(y))
    for _ in range(1, max_iter + 1):
        n_iter += 1
        r = y - X @ coef
        s = 1.4826 * float(np.median(np.abs(r - np.median(r))))
        if s < 1e-12:
            s = 1.0
        u = r / s
        w = np.where(np.abs(u) <= k, 1.0, k / np.maximum(np.abs(u), 1e-300))
        sw = np.sqrt(w)[:, None]
        new_coef, *_ = np.linalg.lstsq(X * sw, y * np.sqrt(w), rcond=None)
        if np.max(np.abs(new_coef - coef)) < tol * (1.0 + np.max(np.abs(coef))):
            coef = new_coef
            converged = True
            break
        coef = new_coef
    return {'coef': coef, 'weights': w, 'n_iter': n_iter,
            'converged': converged}


# ---------------------------------------------------------------------------
# 调整程序 2: 交叉验证岭拟合(共线性稳定化)
# ---------------------------------------------------------------------------

def _ridge_coefficients(X: np.ndarray, y: np.ndarray,
                        alpha: float) -> np.ndarray:
    """Fit ridge on standardized features and return original-scale coefs."""
    features = X[:, 1:]
    mean = features.mean(axis=0)
    scale = features.std(axis=0)
    scale[scale < 1e-15] = 1.0
    Z = (features - mean) / scale
    design = np.column_stack([np.ones(len(y)), Z])
    penalty = np.eye(design.shape[1])
    penalty[0, 0] = 0.0
    beta = np.linalg.pinv(
        design.T @ design + float(alpha) * penalty
    ) @ design.T @ y
    coef = np.empty_like(beta)
    coef[1:] = beta[1:] / scale
    coef[0] = beta[0] - float(np.sum(beta[1:] * mean / scale))
    return coef


def ridge_fit(X: np.ndarray, y: np.ndarray,
              alphas: np.ndarray | None = None,
              n_splits: int = 5,
              random_state: int = 0) -> dict:
    """Choose ridge strength by deterministic k-fold validation."""
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 4:
        raise ValueError("岭调整至少需要 4 个样本")
    if alphas is None:
        alphas = np.logspace(-6, 3, 25)
    alphas = np.asarray(alphas, dtype=float)
    splits = min(max(2, int(n_splits)), n)
    order = np.random.default_rng(random_state).permutation(n)
    folds = np.array_split(order, splits)

    errors = []
    for alpha in alphas:
        squared_error = 0.0
        count = 0
        for valid in folds:
            train = np.setdiff1d(order, valid, assume_unique=True)
            coef = _ridge_coefficients(X[train], y[train], float(alpha))
            residual = y[valid] - X[valid] @ coef
            squared_error += float(np.sum(residual ** 2))
            count += len(valid)
        errors.append(squared_error / max(count, 1))

    best_index = int(np.argmin(errors))
    alpha = float(alphas[best_index])
    coef = _ridge_coefficients(X, y, alpha)
    return {
        'coef': coef,
        'alpha': alpha,
        'cv_rmse_log': float(np.sqrt(errors[best_index])),
        'alpha_grid': [float(value) for value in alphas],
    }


# ---------------------------------------------------------------------------
# 组装: 守卫式拟合入口
# ---------------------------------------------------------------------------

def _logspace_candidate(family: str, X: np.ndarray, log_y: np.ndarray,
                        target: str, features: list[str],
                        coef: np.ndarray, weights: np.ndarray | None = None
                        ) -> CandidateModel:
    r2 = _ols_r2(X, log_y, coef)
    a = float(np.exp(coef[0]))
    exponents = {f: float(b) for f, b in zip(features, coef[1:], strict=True)}
    params = {'prefactor': a, 'exponents': exponents}
    if weights is not None:
        params['n_downweighted'] = int(np.sum(weights < 0.999))
        params['min_weight'] = float(np.min(weights))
    parts = [f"{a:.4g}"] + [f"({f})^{b:.3g}" for f, b in exponents.items()]
    return CandidateModel(
        family=family,
        formula=f"{target} = " + " * ".join(parts),
        params=params,
        success=True,
        r_squared=float(r2),
        n_params=len(features) + 1,
        n_samples=len(log_y),
    )


def fit_with_guard(data: pd.DataFrame, target: str, features: list[str],
                   adjustment: str = 'report_only',
                   mad_threshold: float = 3.5,
                   leverage_factor: float = 2.0,
                   cond_threshold: float = 30.0,
                   huber_k: float = 1.345,
                   exclude_leverage: bool = False) -> dict:
    """
    守卫式拟合: 检测 → (可选)调整 → 重拟合 → 全程留痕。

    参数
    ----
    adjustment: 'report_only'(默认) | 'robust' | 'ridge' | 'exclude' | 'auto'
    mad_threshold:      残差离群的 MAD z 阈值(默认 3.5)
    leverage_factor:    杠杆阈值 = leverage_factor · p/n(默认 2.0)
    cond_threshold:     标准化条件数告警阈值(默认 30)
    huber_k:            Huber 调谐常数(默认 1.345)
    exclude_leverage:   'exclude' 策略是否连同高杠杆点一起剔除。
                        默认 False——残差离群是疑似测量错误(该剔),
                        高杠杆只是条件极端(是设计特征不是错误),
                        默认只报告不剔除, 保住最有信息量的数据点。

    返回 dict:
        success      拟合链路是否完成(≠ 无异常)
        health       数据层问题清单(有 fatal 项则 success=False)
        anomalies    {residual_outliers, leverage, collinearity}
        before       OLS 基准候选(family='power_law')
        after        调整后候选(report_only 时为 None)
        comparison   before/after 的 R² 与系数位移(有 after 时)
        audit        调整审计(剔除清单/降权统计/拒绝原因)
    """
    if adjustment not in ('report_only', 'robust', 'ridge', 'exclude', 'auto'):
        return {'success': False,
                'reason': f"未知调整策略 '{adjustment}',"
                          "可选: report_only / robust / ridge / exclude / auto"}

    health = check_data_health(data, target, features)
    if any(i['severity'] == 'fatal' for i in health):
        return {'success': False, 'health': health,
                'reason': '数据层存在致命问题, 守卫拒绝拟合(先修数据再看异常)'}

    y = data[target].to_numpy(dtype=float)
    X_feat = np.column_stack([np.log(data[f].to_numpy(dtype=float))
                              for f in features])
    log_y = np.log(y)
    X = np.column_stack([np.ones(len(y)), X_feat])
    n, p = X.shape

    # ---- 基准: 普通 OLS(统计口径统一走共享引擎) ----
    coef_ols = _ols(X, log_y)["coef"]
    before = _logspace_candidate('power_law', X, log_y, target, features, coef_ols)

    # ---- 三层检测 ----
    resid = log_y - X @ coef_ols
    outliers = detect_residual_outliers(resid, threshold=mad_threshold)
    lev = detect_leverage(X, threshold=leverage_factor * p / n)
    collin = detect_collinearity(X_feat, cond_threshold=cond_threshold)

    detected = sorted(set(outliers['indices']) | set(lev['indices']))
    flagged = sorted(set(outliers['indices']) |
                     (set(lev['indices']) if exclude_leverage else set()))
    applied_adjustment = adjustment
    if adjustment == 'auto':
        if outliers['indices']:
            applied_adjustment = 'robust'
            decision = '检测到残差离群点:自动采用 Huber 稳健调整'
        elif collin is not None:
            applied_adjustment = 'ridge'
            decision = '检测到特征共线性:自动采用交叉验证岭调整'
        else:
            applied_adjustment = 'report_only'
            decision = '未发现需要后端干预的问题:保留原始拟合'
    else:
        decision = f"调用方指定调整策略:{adjustment}"

    result = {
        'success': True,
        'adjustment': adjustment,
        'applied_adjustment': applied_adjustment,
        'health': health,
        'anomalies': {'residual_outliers': outliers,
                      'leverage': lev,
                      'collinearity': collin},
        'before': before.to_dict(),
        'after': None,
        'comparison': None,
        'audit': {'decision': decision,
                  'flagged_indices': flagged,
                  'flagged_sources': {str(i): (
                      'residual' if i in outliers['indices'] else '') +
                      ('+leverage' if i in lev['indices'] else '')
                      for i in flagged}},
    }

    if not detected and collin is None and applied_adjustment == 'report_only':
        result['audit']['note'] = '未检测到离群/高杠杆点, 无需调整'
        return result

    if applied_adjustment == 'report_only':
        result['audit']['note'] = "report_only 策略: 只检测与报告, 不调整"
        return result

    if applied_adjustment == 'exclude' and not flagged:
        result['audit']['note'] = (
            f"检测到 {len(detected)} 个高杠杆点(条件极端, 非数据错误), "
            "默认只报告不剔除; 如需连同剔除请传 exclude_leverage=True。")
        return result
    # adjustment == 'robust' 永不早退: 全局重拟合策略, 杠杆点若残差小
    # 会被自然保留权重 ≈ 1, 无需点级剔除语义参与

    # ---- 调整策略 ----
    if applied_adjustment == 'robust':
        ir = robust_irls(X, log_y, k=huber_k)
        after = _logspace_candidate('power_law_robust', X, log_y,
                                    target, features, ir['coef'],
                                    weights=ir['weights'])
        result['after'] = after.to_dict()
        result['audit']['irls'] = {'n_iter': ir['n_iter'],
                                   'converged': ir['converged'],
                                   'n_downweighted': int(np.sum(ir['weights'] < 0.999))}

    elif applied_adjustment == 'ridge':
        ridge = ridge_fit(X, log_y)
        after = _logspace_candidate(
            'power_law_ridge', X, log_y, target, features, ridge['coef'])
        after.params['ridge_alpha'] = ridge['alpha']
        after.params['cv_rmse_log'] = ridge['cv_rmse_log']
        result['after'] = after.to_dict()
        result['audit']['ridge'] = {
            'alpha': ridge['alpha'],
            'cv_rmse_log': ridge['cv_rmse_log'],
            'selection': 'deterministic_kfold',
        }

    elif applied_adjustment == 'exclude':
        if len(flagged) > max(1, n // 3):
            result['audit']['refused'] = (
                f"被标记点 {len(flagged)} 个 > 样本 {n} 的 1/3, "
                "拒绝自动剔除——这多半是模型形式错了, 不是数据错了。"
                "请改用 adjustment='robust', 或回到函数族选择层面。")
            return result
        keep = np.array([i not in set(flagged) for i in range(n)])
        Xk, yk = X[keep], log_y[keep]
        coef_k = _ols(Xk, yk)["coef"]
        after = _logspace_candidate('power_law_excluded', Xk, yk,
                                    target, features, coef_k)
        after.n_samples = int(keep.sum())
        result['after'] = after.to_dict()
        result['audit']['excluded'] = [
            {'index': int(i),
             'target_value': float(data[target].to_numpy(dtype=float)[i]),
             'source': result['audit']['flagged_sources'][str(i)],
             'residual_z': outliers['scores'].get(i),
             'leverage': lev['scores'].get(i)}
            for i in flagged]

    # ---- before/after 对比 ----
    if result['after'] is not None:
        aft = result['after']
        shift = {f: abs(aft['params']['exponents'][f]
                        - before.params['exponents'][f])
                 for f in features}
        result['comparison'] = {
            'r2_before': before.r_squared,
            'r2_after': aft['r_squared'],
            'max_exponent_shift': float(max(shift.values())) if shift else 0.0,
            'exponent_shift': shift,
        }

    return result
