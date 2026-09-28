"""
dimreflex.families
==================

函数族基底类 + 统一候选模型 schema + 函数族注册表。

对应《项目评价与改进建议》2.1(函数族三步走)与第四节(落地建议):

- CandidateModel: 所有函数族的同构输出 —— {family, formula, params,
  success, r_squared, n_params, n_samples, reason},各函数族/发现层
  产出同构候选,上层(排序、验证、报告)无需关心具体家族。
- FunctionFamily: 基底类(模板方法模式)。子类只实现 _fit();
  fit() 统一兜底——捕获一切数学异常转为显式失败、为无原因的失败补默认
  原因。这把"数据不足时显式降级,绝不硬输出"的框架哲学从调用方责任
  变成了基类保证:任何函数族实现都不可能向调用方抛数学异常。
- PowerLawFamily: 现有 fit_power_law 的插件化迁移,第一个函数族。

未来扩展点(2.1 三步走,均继承 FunctionFamily):
- PiecewisePowerLawFamily  —— 分段幂律(变点检测)
- RationalLogFamily        —— log 空间有理函数(Gnielinski 类非纯幂律)
- SymbolicRegressionFamily —— PySR 符号回归(算子限 {+, ×, ^})
"""

from dataclasses import asdict, dataclass, fields

import numpy as np
import pandas as pd

from ._ols import ols as _ols

# ---------------------------------------------------------------------------
# 统一候选模型 schema
# ---------------------------------------------------------------------------

@dataclass
class CandidateModel:
    """
    统一候选模型 —— 所有函数族的同构输出。

    成功时:formula 是可读公式,params 携带系数(结构由家族自定义),
    r_squared / n_params / n_samples 填实,supports 为空。
    失败时:success=False,reason 说明原因,绝不抛异常、绝不硬输出。
    """
    family: str
    formula: str
    params: dict
    success: bool
    r_squared: float = float('nan')
    n_params: int = 0
    n_samples: int = 0
    reason: str = ""

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> 'CandidateModel':
        # 宽松反序列化:忽略未知字段(向前兼容),缺省字段用默认值
        known = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in d.items() if k in known})


# ---------------------------------------------------------------------------
# 基底类
# ---------------------------------------------------------------------------

class FunctionFamily:
    """
    所有函数族的基底类(模板方法模式)。

    契约:
        - 子类实现 _fit(),返回 CandidateModel;允许抛任何异常。
        - fit() 是对外唯一入口:统一兜底,把子类的崩溃、NaN、非正值等
          一切失败都折叠为 success=False 的候选模型。

    这样上层(discover_multi / 注册表循环)可以无差别地遍历所有函数族:
    一个坏家族只会贡献一个带 reason 的失败候选,不会中断整个发现流程。
    """
    name: str = "base"

    def fit(self, data: pd.DataFrame, target: str,
            features: list[str]) -> CandidateModel:
        """模板方法:调用 _fit 并兜底。绝不向调用方抛异常。"""
        try:
            cm = self._fit(data, target=target, features=features)
        except Exception as exc:  # noqa: BLE001 —— 兜底是本方法的职责
            return CandidateModel(
                family=self.name, formula="", params={},
                success=False,
                reason=f"{type(exc).__name__}: {exc}",
                n_samples=len(data),
            )
        if not cm.success and not cm.reason:
            cm.reason = "拟合失败(家族未提供原因)"
        return cm

    def _fit(self, data: pd.DataFrame, target: str,
             features: list[str]) -> CandidateModel:
        raise NotImplementedError


def _failure(family: str, reason: str, n_samples: int) -> CandidateModel:
    """构造显式失败候选的小工厂,避免各家族重复样板。"""
    return CandidateModel(family=family, formula="", params={},
                          success=False, reason=reason, n_samples=n_samples)


# ---------------------------------------------------------------------------
# 函数族 1:幂律(现有逻辑的插件化迁移)
# ---------------------------------------------------------------------------

class PowerLawFamily(FunctionFamily):
    """
    幂律族:y = a · ∏ x_i^{b_i}

    log-log 线性最小二乘,适合 Dittus-Boelter 类关联式。
    是原 discover.py fit_power_law 的插件化迁移,行为保持一致。
    """
    name = "power_law"

    def _fit(self, data, target, features) -> CandidateModel:
        y = data[target].values
        X = data[features].values

        if (not np.isfinite(y).all()) or (not np.isfinite(X).all()):
            return _failure(self.name,
                            '数据中存在非有限值(NaN/inf),无法在 log 空间拟合幂律',
                            len(data))

        if (y <= 0).any() or (X <= 0).any():
            return _failure(self.name,
                            '存在非正值,无法在 log 空间拟合幂律',
                            len(data))

        log_y = np.log(y)
        log_X = np.log(X)

        # 加一列常数项(对应 log(a)),最小二乘;统计计算统一走共享引擎
        A = np.column_stack([np.ones(len(log_y)), log_X])
        fit = _ols(A, log_y)
        coef = fit["coef"]
        r2 = fit["r_squared"]
        a = float(np.exp(coef[0]))
        exponents = {f: float(b) for f, b in zip(features, coef[1:], strict=True)}
        n, p = A.shape
        se = fit["std_errors"]
        rank = fit["rank"]

        parts = [f"{a:.4g}"] + [f"({f})^{b:.3g}" for f, b in exponents.items()]
        formula = f"{target} = " + " * ".join(parts)

        return CandidateModel(
            family=self.name,
            formula=formula,
            params={
                'prefactor': a,
                'exponents': exponents,
                'std_errors': {k: float(v) for k, v
                               in zip(['log_a', *list(features)], se,
                                      strict=True)},
                'design_rank': int(rank),
                'design_columns': int(p),
                'identifiable': bool(rank == p),
                'condition_number': float(np.linalg.cond(A)),
            },
            success=True,
            r_squared=float(r2),
            n_params=len(features) + 1,
            n_samples=n,
        )


# ---------------------------------------------------------------------------
# 函数族 2:分段幂律(2.1 第一步:变点检测)
# ---------------------------------------------------------------------------

def _logspace_ols(log_y: np.ndarray, A: np.ndarray) -> dict:
    """一次 log 空间最小二乘,返回系数与 R²。供分段搜索内层循环复用。"""
    fit = _ols(A, log_y)
    return {'coef': fit["coef"], 'ss_res': fit["ss_res"], 'r2': fit["r_squared"]}


class PiecewisePowerLawFamily(FunctionFamily):
    """
    分段幂律族:沿第一个特征(工况变量)检测变点,每段是完整的多特征幂律。

        y = a_k · ∏_i x_i^{b_{k,i}} ,  分段边界由 features[0] 的变点决定

    对应 2.1 第一步:识别工况转变(如层流/湍流 Re≈2300)导致的
    分段关联式。实现为穷举单变点搜索(不引入 pwlf 依赖):

        1. 对每个内部候选切点(排序索引网格,每段至少 m 个样本,
           m = max(min_segment_frac·n, p+2) 保证段内设计矩阵可辨识):
           两段各自做 log-log 线性拟合(设计矩阵 [1, log x_1, …, log x_p]),
           累计残差平方和;
        2. 取 SSE 最小的切点,与"不分段"的单幂律做 AIC 比较
           (AIC = n·ln(SSE/n) + 2k):ΔAIC > 10 才分段——这是主门槛,
           AIC 惩罚项本身抑制过拟合;ΔR²(≥ min_r2_gain,默认 0.01)仅作
           防退化下限。宽动态范围数据上 log 空间 R² 天然偏高,
           绝对 ΔR² 不能作为主门槛;
        3. 门槛不达标 → 如实报告单段——绝不虚报变点。

    变点的可识别性:数据只能把过渡约束在"左段最后一点"与"右段第一点"
    之间的间隙里。报告的 breakpoints 取该间隙的几何中点,并把间隙本身
    存入 params['breakpoint_interval'],不假装拥有区间外的分辨率。

    诚实声明(体现在 params['note']):
        - 变点仅沿第一特征搜索(工况变量);其余特征不参与分段边界判定,
          以完整幂律指数进入每段(段间指数允许不同)。
    """
    name = "piecewise_power_law"

    MIN_SAMPLES = 10  # 低于此样本数连单段拟合都不可信,显式失败

    def __init__(self,
                 max_breaks: int = 1,
                 min_r2_gain: float = 0.01,
                 min_segment_frac: float = 0.2):
        self.max_breaks = max_breaks
        self.min_r2_gain = min_r2_gain
        self.min_segment_frac = min_segment_frac

    # -- 内部:在索引区间 [lo, hi) 上拟合 log(y) = c0 + Σ ci·log(x_i) ------

    def _fit_segment(self, log_X, log_y, lo, hi) -> dict:
        A = np.column_stack([np.ones(hi - lo), log_X[lo:hi]])
        out = _logspace_ols(log_y[lo:hi], A)
        out['a'] = float(np.exp(out['coef'][0]))
        out['lo'], out['hi'] = lo, hi
        return out

    def _fit(self, data, target, features) -> CandidateModel:
        n = len(data)
        if not features:
            return _failure(self.name, '至少需要一个特征列', n)
        p = len(features)
        y = data[target].values.astype(float)
        X = data[features].values.astype(float)

        # 统一防线:非有限值 → 显式失败(与全框架一致)
        if (not np.isfinite(y).all()) or (not np.isfinite(X).all()):
            return _failure(self.name,
                            '数据中存在非有限值(NaN/inf),无法拟合分段幂律', n)
        if (y <= 0).any() or (X <= 0).any():
            return _failure(self.name,
                            '存在非正值,无法在 log 空间拟合分段幂律', n)

        order = np.argsort(X[:, 0])
        X, y = X[order], y[order]
        log_X, log_y = np.log(X), np.log(y)

        if n < max(self.MIN_SAMPLES, p + 3):
            return _failure(self.name,
                            f'样本数 {n} < {max(self.MIN_SAMPLES, p + 3)},'
                            f'不足以支持 {p} 特征的分段拟合(每段至少需 '
                            f'{p + 2} 个样本)', n)

        note = None

        # ---- 基准:单段完整幂律(k = p+1) ----
        base = self._fit_segment(log_X, log_y, 0, n)
        base_aic = n * np.log(base['ss_res'] / n) + 2 * (p + 1)

        if self.max_breaks < 1:
            return self._single_segment_result(features, target, base, n, note)

        # ---- 变点网格:每段至少 m 个样本;小样本逐索引穷举,大样本封顶 200 ----
        m = max(int(self.min_segment_frac * n), p + 2)
        lo_i, hi_i = m, n - m          # 切点索引允许的范围
        if hi_i - lo_i < 3:
            note = (f"样本量不足以在每段保留 {m} 个样本的同时搜索变点,"
                    "如实报告为单段幂律。")
            return self._single_segment_result(features, target, base, n, note)
        best = None
        if hi_i - lo_i <= 200:
            cand = np.arange(lo_i, hi_i)
        else:
            cand = np.unique(np.linspace(lo_i, hi_i, 200).astype(int))
        two_seg_k = 2 * (p + 1) + 1    # 两段各 p+1 参数 + 1 切点
        for c in cand:
            left = self._fit_segment(log_X, log_y, 0, c)
            right = self._fit_segment(log_X, log_y, c, n)
            sse = left['ss_res'] + right['ss_res']
            if not np.isfinite(sse) or sse <= 0:
                continue
            aic = n * np.log(sse / n) + 2 * two_seg_k
            if best is None or aic < best['aic']:
                best = {'aic': aic, 'cut': int(c),
                        'left': left, 'right': right, 'sse': sse}

        if best is None:
            return self._single_segment_result(features, target, base, n, note)

        # ---- 双准则:ΔAIC > 10 且 ΔR² 达标,否则如实报告单段 ----
        ss_tot = float(np.sum((log_y - log_y.mean()) ** 2))
        d_aic = base_aic - best['aic']
        d_r2 = (1 - best['sse'] / ss_tot) - base['r2']
        if d_aic <= 10 or d_r2 < self.min_r2_gain:
            if note is None:
                note = (f"变点增益未达阈值(ΔAIC={d_aic:.1f}, "
                        f"ΔR²={d_r2:.3f}),如实报告为单段幂律。")
            return self._single_segment_result(features, target, base, n, note)

        bp_cut = int(best['cut'])
        x1 = X[:, 0]
        bp = float(np.sqrt(x1[bp_cut - 1] * x1[bp_cut]))  # 间隙的几何中点
        bp_lo, bp_hi = float(x1[bp_cut - 1]), float(x1[bp_cut])
        L, R = best['left'], best['right']
        feat0 = features[0]

        def seg_formula(seg) -> str:
            parts = [f"{seg['a']:.4g}"] + [
                f"({f})^{b:.3g}"
                for f, b in zip(features, seg['coef'][1:], strict=True)]
            return " * ".join(parts)

        def seg_params(seg, range_str) -> dict:
            entry = {
                'range': range_str, 'prefactor': seg['a'],
                'exponents': {f: float(b) for f, b
                              in zip(features, seg['coef'][1:], strict=True)},
            }
            if p == 1:
                # 向后兼容:单特征保留旧 schema 的标量 exponent 键
                entry['exponent'] = float(seg['coef'][1])
            return entry

        formula = (f"{target} = {seg_formula(L)},  {feat0} < {bp:.4g};\n"
                   f"{target} = {seg_formula(R)},  {feat0} ≥ {bp:.4g}")
        r2 = 1 - best['sse'] / ss_tot

        params = {
            'n_segments': 2,
            'breakpoints': [bp],
            'breakpoint_interval': [(bp_lo, bp_hi)],
            'segments': [
                seg_params(L, f"{feat0} < {bp:.4g}"),
                seg_params(R, f"{feat0} ≥ {bp:.4g}"),
            ],
            'delta_aic': float(d_aic),
        }
        if note:
            params['note'] = note
        return CandidateModel(
            family=self.name, formula=formula, params=params,
            success=True, r_squared=float(r2), n_params=two_seg_k, n_samples=n,
        )

    def _single_segment_result(self, features, target, base, n, note
                               ) -> CandidateModel:
        parts = [f"{base['a']:.4g}"] + [
            f"({f})^{b:.3g}"
            for f, b in zip(features, base['coef'][1:], strict=True)]
        formula = f"{target} = " + " * ".join(parts)
        params = {
            'n_segments': 1,
            'breakpoints': [],
            'segments': [{
                'range': '全部', 'prefactor': base['a'],
                'exponents': {f: float(b) for f, b
                              in zip(features, base['coef'][1:], strict=True)},
            }],
        }
        if len(features) == 1:
            params['segments'][0]['exponent'] = float(base['coef'][1])
        if note:
            params['note'] = note
        return CandidateModel(
            family=self.name, formula=formula, params=params,
            success=True, r_squared=float(base['r2']),
            n_params=len(features) + 1, n_samples=n,
        )


# ---------------------------------------------------------------------------
# 函数族 3:log 空间有理函数(2.1 第二步,Gnielinski 类非纯幂律)
# ---------------------------------------------------------------------------

class RationalLogFamily(FunctionFamily):
    """
    log 空间有理函数族(多特征版):

        log10(y) = num(L) / den(L) + Σ_j b_j·Z_j

        num(L)  = c0 + c1·L + … + c_{p-1}·L^{p-1}
        den(L)  = 1 + q1·L + … + q_d·L^d
        L       = log10(x_1) − L0    (主特征,中心化让条件数友好)
        Z_j     = log10(x_j) − c_j   (其余特征的中心化 log,幂律协变量项)

    分母常数项固定为 1 保证可辨识(有理函数分子分母同乘常数不改变函数)。
    单特征(p=1)退化为原 MVP——num/den 一次性线性化最小二乘,数值路径
    与历史版本逐位一致。

    多特征拟合:分块交替最小二乘——
      (a) 固定协变量系数 b,对残差 log10y − Z·b 拟合有理部分(沿用
          一次性线性化:[V, −r·W]·[c; q] = r);
      (b) 固定有理部分,对中心化残差做 Z 上的 OLS 更新 b(残差均值交还
          给分子常数项 c0 吸收)。
    每块都是条件精确最小二乘,SSE 单调不增;迭代至收敛或上限(默认 50 轮),
    全程确定性、无随机初始化。

    已知限制(诚实声明):交替最小二乘收敛到局部平稳点而非全局最优;
    线性化在噪声较大时有偏(errors-in-variables)。den_degree=0 退化为
    log 多项式 + 线性协变量(幂律是其特例)。分母在数据上过小/发散时
    显式失败,不硬输出。
    """
    name = "rational_log"

    def __init__(self, num_degree: int = 2, den_degree: int = 1,
                 max_iter: int = 50, tol: float = 1e-10):
        self.num_degree = num_degree      # 分子系数个数(阶数 = num_degree-1)
        self.den_degree = den_degree      # 分母系数个数(分母常数 1 之外)
        self.max_iter = max_iter          # 交替最小二乘的迭代上限
        self.tol = tol                    # SSE 相对收敛容差

    def _fit(self, data, target, features) -> CandidateModel:
        n = len(data)
        if not features:
            return _failure(self.name, '至少需要一个特征列', n)
        p = len(features)
        y = data[target].values.astype(float)
        X = data[features].values.astype(float)

        if (not np.isfinite(y).all()) or (not np.isfinite(X).all()):
            return _failure(self.name,
                            '数据中存在非有限值(NaN/inf),无法拟合有理函数', n)
        if (y <= 0).any() or (X <= 0).any():
            return _failure(self.name,
                            '存在非正值,无法在 log 空间拟合有理函数', n)

        n_params = self.num_degree + self.den_degree + (p - 1)
        if n < n_params + 2:
            return _failure(self.name,
                            f'样本数 {n} 不足以支撑 {self.num_degree}+'
                            f'{self.den_degree}+{p - 1} 参数的有理拟合'
                            f'(至少需 {n_params + 2} 个样本)', n)

        log10y = np.log10(y)
        log10X = np.log10(X)
        L0 = float(log10X[:, 0].mean())
        L = log10X[:, 0] - L0

        if p > 1:
            Z_centers = log10X[:, 1:].mean(axis=0)
            Z = log10X[:, 1:] - Z_centers
        else:
            Z_centers = np.zeros(0)
            Z = np.zeros((n, 0))

        V = np.vander(L, N=self.num_degree, increasing=True)      # n × p
        # W[:, k-1] = L^k,供分母 1 + Σ q_k·L^k 使用
        W = (L[:, None] ** np.arange(1, self.den_degree + 1)
             if self.den_degree else np.zeros((n, 0)))

        def fit_rational(response):
            """一次性线性化:(1 + W·q)·r̂ = V·c → [V, −r·W]·[c; q] = r。"""
            A = np.hstack([V, -W * response[:, None]])
            sol, _, _, _ = np.linalg.lstsq(A, response, rcond=None)
            if not np.isfinite(sol).all():
                return None
            c, q = sol[:self.num_degree], sol[self.num_degree:]
            den = 1.0 + (W @ q if self.den_degree else 0.0)
            if ((not np.isfinite(den).all())
                    or float(np.min(np.abs(den))) <= 0.05):
                return None
            return c, q, den

        c = q = den = None
        b = np.zeros(p - 1)
        converged = True
        if p == 1:
            # 单特征:原 MVP 数值路径,保持逐位一致
            rat = fit_rational(log10y)
            if rat is None:
                return _failure(self.name, '线性化求解发散,拟合失败', n)
            c, q, den = rat
            pred = (V @ c) / den
            sse = float(np.sum((log10y - pred) ** 2))
        else:
            # 多特征:分块交替最小二乘
            sse = np.inf
            converged = False
            for _ in range(self.max_iter):
                rat = fit_rational(log10y - Z @ b)
                if rat is None:
                    if not np.isfinite(sse):
                        return _failure(self.name,
                                        '线性化求解发散,拟合失败', n)
                    break  # 迭代中发散:保留上一轮有效解
                c, q, den = rat
                resid = log10y - (V @ c) / den
                # 残差中心化:均值由分子常数项 c0 在下一轮吸收
                b_new, _, _, _ = np.linalg.lstsq(
                    Z, resid - resid.mean(), rcond=None)
                if not np.isfinite(b_new).all():
                    break
                b = b_new
                new_sse = float(np.sum(
                    (log10y - (V @ c) / den - Z @ b) ** 2))
                if abs(sse - new_sse) <= self.tol * (1.0 + new_sse):
                    sse = new_sse
                    converged = True
                    break
                sse = new_sse
            if c is None:
                return _failure(self.name, '线性化求解发散,拟合失败', n)
            pred = (V @ c) / den + Z @ b

        ss_res = sse
        ss_tot = float(np.sum((log10y - log10y.mean()) ** 2))
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
        aic = float(n * np.log(ss_res / n) + 2 * n_params)

        num_terms = " + ".join(f"{ci:+.4g}·L^{k}" if k else f"{ci:+.4g}"
                               for k, ci in enumerate(c))
        if self.den_degree:
            den_terms = " + ".join(f"{qi:+.4g}·L^{k}"
                                   for k, qi in enumerate(q, start=1))
            expr = f"({num_terms}) / (1 {den_terms})"
        else:
            expr = f"({num_terms})"
        if p > 1:
            extra_terms = " + ".join(
                f"{bi:+.4g}·(log10({f}) − {cj:.3g})"
                for f, bi, cj in zip(features[1:], b, Z_centers, strict=True))
            expr = f"{expr} + [{extra_terms}]"

        def predict_log10(Lv, Zv=None, _c=c, _q=q, _p=self.num_degree,
                          _d=self.den_degree, _b=b):
            # 标量入 → 标量出;数组入 → 数组出。Zv 为其余特征的中心化
            # log10 值(一维行向量;单特征模型调用时省略)。
            scalar_input = np.ndim(Lv) == 0
            Lv_arr = np.atleast_1d(np.asarray(Lv, dtype=float))
            out = np.vander(Lv_arr, N=_p, increasing=True) @ _c
            if _d:
                out = out / (1.0 + sum(qi * Lv_arr ** k
                                       for k, qi in enumerate(_q, start=1)))
            if Zv is not None and len(_b):
                zv = np.atleast_1d(np.asarray(Zv, dtype=float))
                out = out + zv @ _b
            return float(out[0]) if scalar_input else out

        params = {
            'num_degree': self.num_degree,
            'den_degree': self.den_degree,
            'num_coefs': [float(v) for v in c],
            'den_coefs': [float(v) for v in np.atleast_1d(q)],
            'log_center': L0,
            'aic': aic,
            'predict_log10': predict_log10,
        }
        if p > 1:
            params['extra_features'] = list(features[1:])
            params['extra_coefs'] = [float(v) for v in b]
            params['extra_centers'] = [float(v) for v in Z_centers]
            params['alternating_ls_converged'] = converged

        return CandidateModel(
            family=self.name,
            formula=(f"log10({target}) = {expr},  "
                     f"L = log10({features[0]}) − {L0:.3g}"),
            params=params,
            success=True,
            r_squared=float(r2),
            n_params=n_params,
            n_samples=n,
        )


# ---------------------------------------------------------------------------
# 函数族 4:符号回归占位(2.1 第三步;显式失败优于假装实现)
# ---------------------------------------------------------------------------

class SymbolicRegressionFamily(FunctionFamily):
    """
    符号回归族占位(2.1 第三步:接入 PySR,算子限 {+, ×, ^})。

    PySR 未安装时显式失败并说明安装方式;已安装时也如实报告
    本占位尚未实现搜索配置——绝不假装成功。这符合框架
    "数据不足/能力不足时显式降级,绝不硬输出"的哲学。
    """
    name = "symbolic"

    def _fit(self, data, target, features) -> CandidateModel:
        try:
            import pysr  # noqa: F401
            available = True
        except ImportError:
            available = False
        if not available:
            return _failure(
                self.name,
                'PySR 未安装:符号回归(改进建议 2.1 第三步)依赖 pysr 包,'
                '请 `pip install pysr` 后重试;未安装时本家族显式失败',
                len(data))
        return _failure(
            self.name,
            'PySR 已安装,但 SymbolicRegressionFamily 尚为实现完整搜索'
            '配置的占位(算子集 {+, ×, ^} 与 parsimony 约束待定)',
            len(data))


# ---------------------------------------------------------------------------
# 注册表 + 多族发现(文档第四节:注册表模式,各层产出同构候选)
# ---------------------------------------------------------------------------

FAMILIES: dict[str, type] = {}


def register_family(family_cls: type) -> type:
    """把函数族类注册进 FAMILIES(以类属性 name 为键)。支持装饰器用法。"""
    key = getattr(family_cls, 'name', None)
    if not key:
        raise ValueError(f"{family_cls.__name__} 缺少 name 类属性,无法注册")
    FAMILIES[key] = family_cls
    return family_cls


for _builtin in (PowerLawFamily, PiecewisePowerLawFamily,
                 RationalLogFamily, SymbolicRegressionFamily):
    register_family(_builtin)


def discover_multi(data: pd.DataFrame,
                   target: str,
                   features: list[str],
                   families: list[str] | None = None) -> list[CandidateModel]:
    """
    对所有(或指定)注册函数族运行发现,返回同构候选列表。

    任何家族的崩溃都已被 FunctionFamily.fit() 兜底为失败候选,
    因此本函数绝不因单个家族的问题中断。排序:成功的候选按
    R² 降序在前,失败候选(带 reason)排在后面供诊断。
    """
    names = families if families is not None else list(FAMILIES)
    candidates: list[CandidateModel] = []
    for name in names:
        if name not in FAMILIES:
            candidates.append(CandidateModel(
                family=name, formula="", params={}, success=False,
                reason=f"未注册的函数族 '{name}'", n_samples=len(data)))
            continue
        candidates.append(FAMILIES[name]().fit(data, target, features))

    successes = [c for c in candidates
                 if c.success and np.isfinite(c.r_squared)]
    failures = [c for c in candidates
                if not (c.success and np.isfinite(c.r_squared))]
    successes.sort(key=lambda c: c.r_squared, reverse=True)
    return successes + failures
