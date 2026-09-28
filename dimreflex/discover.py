"""
dimreflex.discover
==================

发现层:在无量纲化的数据上寻找候选函数关系。

当前实现:
- log-log 线性回归 (适合幂律形式,如 Nu = a * Re^b * Pr^c)
- 输出候选函数 + R^2 + 系数置信区间
- method='multi' / discover_multi_pi: 在 π 空间运行全部注册函数族
  (函数族定义见 families.py)
"""


import numpy as np
import pandas as pd

from ._ols import ols as _ols
from .autodiscovery import assess_family as _assess_family
from .autodiscovery import discover_auto as _discover_auto
from .families import discover_multi


def compute_groups(data: pd.DataFrame, groups: list[dict[str, float]]) -> pd.DataFrame:
    """
    根据 dimensionless_groups() 返回的群定义,从原始数据计算无量纲数值。

    每个 group 是 {var_name: exponent},无量纲值 = ∏ data[var]^exponent
    """
    out = pd.DataFrame(index=data.index)
    # 零点取负指数会产生 inf(奇异群)——这是合法情形,由 _prepare_pi 的
    # 逐特征降级分类处理,此处不告警
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        for i, group in enumerate(groups):
            col = np.ones(len(data))
            for var, exp in group.items():
                col = col * (data[var].values ** exp)
            # 给列起一个可读的名字
            name = '_'.join(f"{var}^{exp:g}" for var, exp in group.items())
            out[f"PI_{i}_{name}"] = col
    return out


def fit_power_law(data: pd.DataFrame,
                  target: str,
                  features: list[str]) -> dict:
    """
    在 log-log 空间拟合: log(target) = log(a) + b1*log(f1) + b2*log(f2) + ...

    这相当于发现幂律关系 target = a * f1^b1 * f2^b2 * ...

    要求 target 和所有 features 都为正(否则 log 没意义)。
    """
    y = data[target].values
    X = data[features].values

    if (not np.isfinite(y).all()) or (not np.isfinite(X).all()):
        return {
            'success': False,
            'reason': '数据中存在非有限值(NaN/inf),无法在 log 空间拟合幂律',
        }

    if (y <= 0).any() or (X <= 0).any():
        return {
            'success': False,
            'reason': '存在非正值,无法在 log 空间拟合幂律',
        }

    log_y = np.log(y)
    log_X = np.log(X)

    # 加一列常数项(对应 log(a));OLS 引擎负责 lstsq/R²/标准误等全部统计
    A = np.column_stack([np.ones(len(log_y)), log_X])
    fit = _ols(A, log_y)

    coef = fit["coef"]
    r2 = fit["r_squared"]
    a = np.exp(coef[0])
    exponents = coef[1:]
    n, p = A.shape
    se = fit["std_errors"]

    # 构建可读的公式字符串
    formula_parts = [f"{a:.4g}"]
    for f, b in zip(features, exponents, strict=True):
        formula_parts.append(f"({f})^{b:.3g}")
    formula = f"{target} = " + " * ".join(formula_parts)

    return {
        'success': True,
        'formula': formula,
        'prefactor': a,
        'exponents': dict(zip(features, exponents, strict=True)),
        'r_squared': r2,
        'std_errors': dict(zip(['log_a', *features], se, strict=True)),
        'design_rank': int(fit["rank"]),
        'design_columns': int(p),
        'identifiable': bool(fit["rank"] == p),
        'condition_number': float(fit["condition_number"]),
        'n_samples': n,
    }


def _prepare_pi(data: pd.DataFrame, system, target: str) -> dict:
    """
    端到端流水线的公共前半段(discover 与 discover_multi_pi 共用):

    1. 数据就绪检查(check_data)
    2. 生成无量纲群(自动对 target 做目标隔离:唯一因变量群含 target 且指数为
       +1,其余群不含 target——保证关系式物理可解释,不会"用 target 预测 target")
    3. 计算 π 数值表,定位唯一含 target 的群作为 Y、其余群作为 X

    返回 {'report', 'groups', 'pi_data', 'y_col', 'x_cols', 'note',
          'dropped_features', 'shared_vars'}。
    无法继续时 note 非空,其余字段为默认值——调用方据此显式降级。
    """
    report = system.check_data(data)
    out = {'report': report, 'groups': [], 'pi_data': None,
           'y_col': None, 'x_cols': [], 'note': None,
           'dropped_features': [], 'shared_vars': {}}

    if report['level'] < 2:
        out['note'] = (
            f"数据等级 Level {report['level']},不足以进行函数发现。"
            f"框架不会输出未经验证的拟合结果。")
        return out

    groups = system.dimensionless_groups(target=target)
    out['groups'] = groups

    if not groups:
        out['note'] = "未找到无量纲群——可能变量已经全部独立,或量纲矩阵满秩。"
        return out

    pi_data = compute_groups(data, groups)
    out['pi_data'] = pi_data

    # 目标隔离后应恰好一个群含 target,此处保留防御
    target_groups = [i for i, g in enumerate(groups) if target in g]
    if not target_groups:
        out['note'] = f"目标变量 '{target}' 未出现在任何无量纲群中。"
        return out

    y_idx = target_groups[0]
    out['y_col'] = pi_data.columns[y_idx]
    out['x_cols'] = [c for i, c in enumerate(pi_data.columns) if i != y_idx]

    if not out['x_cols']:
        out['note'] = "只有一个无量纲群,无法进行回归。"
        return out

    # -- 逐特征降级:含非正值/零点奇异的特征群显式剔除,而非让整个流程失败 --
    # 零点往往是合法物理量(如攻角 alpha=0 deg),让一个群毁掉全部函数族
    # 粒度过粗;改为剔除该群、留痕、用剩余群继续拟合。
    # 刻意的不对称:源数据含 NaN/inf 时不剔除,留给拟合层整体失败——
    # 零是物理现实(含 0^负指数→inf 的奇异群),NaN 是数据质量事故,
    # 后者框架保持最大保守。
    admissible, dropped = [], []
    for col in out['x_cols']:
        col_idx = int(pi_data.columns.get_loc(col))
        src = data[list(groups[col_idx])].values
        if not np.isfinite(src).all():
            admissible.append(col)  # 源数据事故:不剔除,走整体失败
            continue
        values = pi_data[col].values
        if (not np.isfinite(values).all()) or (values <= 0).any():
            dropped.append(col)
        else:
            admissible.append(col)
    out['dropped_features'] = dropped
    out['x_cols'] = admissible
    if dropped:
        report['warnings'].append(
            f"特征群 {dropped} 含非正值或在零点处无定义,"
            f"无法在 log 空间使用,已显式剔除;剩余 {len(admissible)} 个"
            f"特征群继续拟合(被剔除群的影响会被剩余群部分吸收,"
            f"解读指数时请注意)。")
    if not admissible:
        out['note'] = ("全部特征群都含非正值,无法在 log 空间回归。"
                       "框架不会输出未经验证的拟合结果。")
        return out

    # -- 共享变量检测:目标群与特征群共享变量会在数学上抬高 R² ------------
    # 这是量纲分析的经典伪相关陷阱(如目标群 p/(f·μ) 与含 f 的特征群,
    # 或教科书情形 Nu 与 Re 共享 D)。框架不禁止这类拟合,但必须显式
    # 标注,不把可能虚高的 R² 静默输出。
    y_vars = set(groups[y_idx]) - {target}
    shared = {}
    for col in out['x_cols']:
        col_idx = int(pi_data.columns.get_loc(col))
        common = sorted(y_vars & set(groups[col_idx]))
        if common:
            shared[col] = common
    out['shared_vars'] = shared
    if shared:
        report['warnings'].append(
            f"目标群与特征群共享变量 {shared}:共享变量会在数学上抬高 R²"
            f"(π 群伪相关),该 R² 不宜与无共享变量的拟合直接比较。")

    # -- 形式检查(v0.5.0):幂律形式是否适合当前数据 ------------------------
    # 钩子放在公共前半段,discover/discover_multi_pi/discover_auto_pi
    # 三个 π 入口共享,防止绕行(首版只接 discover() 就漏了另外两个)。
    # 负对照(validate_friction.py)的预期修复:幂律路径不得对形式失配静默。
    out['form_check'] = form_check(pi_data, out['y_col'], out['x_cols'])
    _apply_form_check(report, out['form_check'])

    return out


_LEVEL3_CV_R2 = 0.9


def _cv_r_squared(data: pd.DataFrame, target: str,
                  features: list[str], family: str):
    """对指定函数族做确定性 k-fold CV,返回 out-of-fold log 空间 R²。

    Level 3 的统一判据由此而来:不再看样本内 R²,而是复用 autodiscovery
    的同一套确定性折叠与预测契约,保证 power_law/multi/auto 三入口
    口径一致。验证未完成/失败时返回 None(不升级、不硬判)。
    """
    try:
        assessment = _assess_family(data, target, features, family)
    except Exception:  # noqa: BLE001 — CV 失败按"验证不可用"降级,不升级
        return None
    if not assessment.get("success"):
        return None
    return assessment.get("validation", {}).get("cv_r_squared")


def _promote_level_3(report: dict, cv_r_squared,
                     in_sample_r_squared=None) -> None:
    """用 out-of-fold CV R² 统一提升数据等级到 Level 3。

    CV 不可用或未达阈值时绝不升级;若样本内 R² 高而 CV 低,把矛盾写进
    warnings(这正是负对照暴露的"统计精确 ≠ 形式正确"情形)。
    """
    if cv_r_squared is None or not np.isfinite(cv_r_squared):
        return
    if cv_r_squared > _LEVEL3_CV_R2:
        report['level'] = 3
        report['warnings'].append(
            f"已生成经确定性 k-fold 验证的候选经验关联式 (Level 3, "
            f"CV R²={float(cv_r_squared):.4f})。"
            "进一步提升至 Level 4 仍需跨工况外部验证数据。")
        return
    if in_sample_r_squared is not None and in_sample_r_squared > _LEVEL3_CV_R2:
        report['warnings'].append(
            f"样本内 R²={float(in_sample_r_squared):.4f} 看似达标,但 "
            f"k-fold CV R²={float(cv_r_squared):.4f} ≤ {_LEVEL3_CV_R2:g}, "
            "候选未升级 Level 3(可能过拟合或形式失配)。")


def form_check(data: pd.DataFrame, target: str, features: list[str],
               n_boot: int = 200, random_state: int = 42) -> dict:
    """幂律形式合适性检查:线性假设的偏残差曲率检验(v0.5.0)。

    [动机] 框架原有检查全部关于"数据够不够"(L0-L3),没有关于
    "幂律形式合不合适"的。负对照(validate_friction.py, Prandtl-Karman
    对数律)实锤:幂律路径会对非幂律真值静默输出 R²=0.99/Level 3,
    而该模型外推误差可达 -20%。本检查补上"形式合不合适"这一维。

    [方法] log-log 空间里幂律=直线。先做常规线性拟合,再给设计矩阵
    加一列二次项(以 x 均值为中心,避免与线性项共线),其系数 c2 就是
    "线性拟合之后残余的弯曲量"(偏残差曲率)——对多特征同样成立。
    H0(真值是幂律)下 c2=0;关键是判"离 0 多远才算远":

    [阈值为什么用残差 bootstrap] 教训(负对照 C5):曲率统计量的噪声
    尺度是数据设计的函数(σ ∝ x 量程跨度的平方,实测 5.9 倍),
    写死常数阈值必然在某些设计上误报、另一些上漏报。因此对当前数据
    做残差重采样:在线性拟合的预测值上加"打乱的残差"重建 n_boot 个
    伪数据集——它们继承了本数据的 x 网格、样本量、噪声水平,唯一
    丢掉的是真实曲率。伪数据集的 c2 分布=本设计的零分布,
    t = (c2_obs - mean_null) / std_null。|t|>3 报警。

    残差独立重采样(非整行重抽):残差在小样本下近似同方差(log 空间
    乘性噪声=加性),打乱残差保住零分布又不引入重复行。

    返回(同构候选 schema 风格,失败也是结构化结果,绝不抛异常):
        {'checked': bool, 'alarm': bool, 'curvature': float,
         't_stat': float, 'n_samples': int, 'n_boot': int,
         'reason': str}   # checked=False 时携带原因
    """
    out = {'checked': False, 'alarm': False, 'curvature': 0.0,
           't_stat': 0.0, 'n_samples': 0, 'n_boot': int(n_boot),
           'reason': ''}
    y = data[target].values
    X = data[features].values
    if (not np.isfinite(y).all()) or (not np.isfinite(X).all()):
        out['reason'] = '存在非有限值,形式检查不适用(数据层失败已另行报告)'
        return out
    if (y <= 0).any() or (X <= 0).any():
        out['reason'] = '存在非正值,形式检查不适用(log 空间不可达)'
        return out
    n = len(y)
    if n < 12:
        out['reason'] = f'样本量 n={n} 过少,曲率检验无功效'
        return out

    lx = np.log(X)
    ly = np.log(y)
    # 二次列必须**平方**:曾误写成 (lx-x_mean)(未平方),那只是截距与
    # x 列的线性组合,设计矩阵秩亏,lstsq 最小范数解给出任意系数,
    # 零分布精确骑在观测值上,t 恒为 0——检查永远沉默。调试指纹:
    # null.mean ≈ c2_obs 即"统计量列不携带独立信息"。
    if features and len(features) == 1:
        quad = (lx[:, 0] - lx[:, 0].mean()) ** 2
    else:
        quad = ((lx - lx.mean(axis=0)) ** 2).mean(axis=1)

    def _curv(response, design):
        coef, _, _, _ = np.linalg.lstsq(design, response, rcond=None)
        return float(coef[-1])          # 最后一列 = 二次项

    base = np.column_stack([np.ones(n), lx])
    design = np.column_stack([base, quad])
    c2_obs = _curv(ly, design)

    pred = base @ np.linalg.lstsq(base, ly, rcond=None)[0]
    resid = ly - pred
    rng = np.random.default_rng(random_state)
    null_c2 = np.empty(n_boot)
    for i in range(n_boot):
        boot_y = pred + rng.permutation(resid)
        null_c2[i] = _curv(boot_y, design)

    sd = float(null_c2.std())
    t_val = float((c2_obs - null_c2.mean()) / sd) if sd > 0 else 0.0
    out.update({'checked': True, 'curvature': c2_obs,
                't_stat': t_val, 'n_samples': int(n)})
    if abs(t_val) > 3.0:
        out['alarm'] = True
        out['reason'] = (
            f"偏残差曲率 t={t_val:+.1f} 超出本数据 bootstrap 零分布 3σ"
            f"(n_boot={n_boot}):幂律形式对当前数据可能不适,建议运行"
            f" method='auto' 比较其他函数族。")
    else:
        out['reason'] = (
            f"偏残差曲率 t={t_val:+.1f} 在本数据 bootstrap 零分布 3σ 内"
            f"(n_boot={n_boot}):未发现反对幂律形式的证据。")
    return out


def _apply_form_check(report: dict, payload: dict) -> None:
    """把形式检查结果挂到 report.warnings(报警才写,留痕不噪)。"""
    if payload.get('checked') and payload.get('alarm'):
        report['warnings'].append(
            "形式失配警告(log-log 线性检查): " + payload['reason'])



def discover(data: pd.DataFrame,
             system,
             target: str,
             method: str = 'power_law') -> dict:
    """
    端到端发现流程:
    1. 检查数据可用性
    2. 生成无量纲群(自动目标隔离)
    3. 在无量纲空间拟合候选函数

    method='power_law' (默认): 仅幂律拟合,返回 'fit' 字典(与历史行为一致)。
    method='multi': 在 π 空间运行全部注册函数族,返回 'candidates' 候选列表
                   (CandidateModel),并把 R² 最高的成功候选放入 'fit'。
    method='auto': k-fold 验证全部函数族,用 AICc/BIC 约束复杂度,自主推荐
                   模型并按需执行可审计的稳健/岭后端调整。

    返回:
        {
            'data_report': System.check_data 的结果
            'groups': 找到的无量纲群定义
            'pi_data': 计算出的无量纲数值表
            'dropped_features': 因含非正值被显式剔除的特征群(留痕)
            'shared_vars': 目标群与特征群共享的变量(R² 可能虚高的提示)
            'form_check': 幂律形式合适性检查结果(偏残差曲率 t 检验,
                         报警时 data_report.warnings 同步留痕)
            'fit': 拟合结果(如果数据足够;method='multi' 时为最优候选字典)
            'candidates': 仅 method='multi' —— 全部函数族候选
        }
    """
    prep = _prepare_pi(data, system, target)
    report = prep['report']
    result = {'data_report': report, 'groups': prep['groups'],
              'pi_data': prep['pi_data'],
              'dropped_features': prep['dropped_features'],
              'shared_vars': prep['shared_vars']}

    if prep['note']:
        result['fit'] = None
        result['note'] = prep['note']
        return result

    # 形式检查已在 _prepare_pi 完成,这里透传给调用方
    fc = result['form_check'] = prep.get('form_check')

    if method == 'power_law':
        fit_result = fit_power_law(prep['pi_data'], prep['y_col'], prep['x_cols'])
        result['fit'] = fit_result
        if fit_result.get('success'):
            cv_r2 = _cv_r_squared(prep['pi_data'], prep['y_col'],
                                 prep['x_cols'], 'power_law')
            _promote_level_3(report, cv_r2,
                             in_sample_r_squared=fit_result.get('r_squared'))
    elif method == 'multi':
        candidates = discover_multi(prep['pi_data'], target=prep['y_col'],
                                    features=prep['x_cols'])
        result['candidates'] = candidates
        best = next((c for c in candidates if c.success), None)
        result['fit'] = best.to_dict() if best is not None else None
        if best is not None:
            cv_r2 = _cv_r_squared(prep['pi_data'], prep['y_col'],
                                  prep['x_cols'], best.family)
            _promote_level_3(report, cv_r2,
                             in_sample_r_squared=best.r_squared)
    elif method == 'auto':
        autonomous = _discover_auto(
            prep['pi_data'], target=prep['y_col'],
            features=prep['x_cols'])
        # auto 路径:推荐了幂律族但形式检查报警时,把矛盾显式写进
        # 推荐理由——用户不该自己去挖 form_check 才知道推荐有争议。
        rec = autonomous.get('recommended') or {}
        if (rec.get('family') == 'power_law' and fc.get('alarm')
                and isinstance(rec.get('selection_reason'), str)):
            rec['selection_reason'] += (
                " [形式检查报警: " + fc['reason'] + "]")
        result['autonomous'] = autonomous
        result['assessments'] = autonomous.get('assessments', [])
        result['fit'] = autonomous.get('final_model')
        if autonomous.get('success'):
            cv_r2 = autonomous['recommended']['validation']['cv_r_squared']
            _promote_level_3(report, cv_r2)
    else:
        result['fit'] = {'success': False, 'reason': f"未知方法: {method}"}

    return result


def discover_auto_pi(data: pd.DataFrame,
                     system,
                     target: str,
                     families: list[str] | None = None,
                     n_splits: int = 5,
                     random_state: int = 42,
                     adjust_backend: bool = True) -> dict:
    """Run validated autonomous discovery in dimensionless pi space."""
    prep = _prepare_pi(data, system, target)
    result = {
        'data_report': prep['report'],
        'groups': prep['groups'],
        'pi_data': prep['pi_data'],
        'dropped_features': prep['dropped_features'],
        'shared_vars': prep['shared_vars'],
        'form_check': prep.get('form_check'),
        'autonomous': None,
        'fit': None,
    }
    if prep['note']:
        result['note'] = prep['note']
        return result

    autonomous = _discover_auto(
        prep['pi_data'], target=prep['y_col'], features=prep['x_cols'],
        families=families, n_splits=n_splits,
        random_state=random_state, adjust_backend=adjust_backend)
    result['autonomous'] = autonomous
    result['fit'] = autonomous.get('final_model')
    if autonomous.get('success'):
        cv_r2 = autonomous['recommended']['validation']['cv_r_squared']
        _promote_level_3(prep['report'], cv_r2)
    return result


def discover_multi_pi(data: pd.DataFrame,
                      system,
                      target: str,
                      families: list[str] | None = None) -> dict:
    """
    端到端发现流程(多函数族版):在无量纲空间运行注册函数族。

    这是函数族漏斗与量纲感知引擎的合流点:所有候选都在 π 空间生成,
    π 空间中的运算符(+、×、^ …)天然量纲合法,量纲一致性免费获得。

    返回:
        {
            'data_report': check_data 结果,
            'groups': 无量纲群定义(已做目标隔离),
            'pi_data': π 数值表,
            'dropped_features': 因含非正值被显式剔除的特征群(留痕),
            'shared_vars': 目标群与特征群共享的变量(R² 可能虚高的提示),
            'candidates': 同构 CandidateModel 列表(成功者按 R² 降序在前),
            'best': R² 最高的成功候选字典(无成功候选时为 None),
        }
    """
    prep = _prepare_pi(data, system, target)
    report = prep['report']
    result = {'data_report': report, 'groups': prep['groups'],
              'pi_data': prep['pi_data'],
              'dropped_features': prep['dropped_features'],
              'shared_vars': prep['shared_vars'],
              'form_check': prep.get('form_check'),
              'candidates': [], 'best': None}

    if prep['note']:
        result['note'] = prep['note']
        return result

    candidates = discover_multi(prep['pi_data'], target=prep['y_col'],
                                features=prep['x_cols'], families=families)
    result['candidates'] = candidates

    best = next((c for c in candidates if c.success), None)
    if best is not None:
        result['best'] = best.to_dict()
        cv_r2 = _cv_r_squared(prep['pi_data'], prep['y_col'],
                              prep['x_cols'], best.family)
        _promote_level_3(report, cv_r2,
                         in_sample_r_squared=best.r_squared)

    return result
