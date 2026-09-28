"""scenarios.py — 随机/系统误差场景库(全部确定性,seed 全程传递)。

设计原则:
- 每个场景 = 一个 python 函数 seed->DataFrame(target, x1, [x2]) + 元数据。
- H0 场景(target 是特征的严格单项幂律,噪声层任意):H0_BASELINE,
  H0_MULTIVAR, H0_CORR, H0_MULTIFIT(三 H0 均为严格单项式真值)。
- H1 场景(系统误差):QUAD_LIN, CUBIC_LIN, INTERACT, THRESHOLD,
  EXP_MODEL, LOG_MODEL, RATIONAL, OMITTED_M、H1_HEAVYTAIL,
  H1_OUTLIER, H1_OUTLIER_GRP, H1_VAR_GROWING(方差层 H1),
  DILUTION(稀释), OPPOSING(反向曲率)。
- 真函数经 y = exp(f_true(lx)) 进入;噪声以相对尺度加在 ly 上,
  乘性进入 y(相对噪声语义与 validate_friction 一致)。
"""
import numpy as np
import pandas as pd


def _finalize(ly, rng, rel):
    """噪声层:ly += rel * 噪声,再转回原空间。返回 (y, y_raw_clean)。"""
    ly_c = ly + rel * rng.standard_normal(len(ly))
    y = np.exp(ly_c)
    y_clean = np.exp(ly)
    return y, y_clean

def _domain(n, seed, x1_lo=0.1, x1_hi=10.0, x2_lo=0.05, x2_hi=2.0):
    """固定 log 均匀域(量程跨度 ~4.6,与 pilot/validate_friction 同量级)。
    同一 seed 的 x1/x2 网格在同 n 下全场景一致(除共线场景),噪声流独立。"""
    rng = np.random.default_rng(seed * 7919 + 104729)
    lx1 = rng.uniform(np.log(x1_lo), np.log(x1_hi), n)
    lx2 = rng.uniform(np.log(x2_lo), np.log(x2_hi), n)
    return lx1, lx2, rng

def _mkdf(y, lx1, lx2=None):
    if lx2 is None:
        return pd.DataFrame({"target": y, "x1": np.exp(lx1)})
    return pd.DataFrame({"target": y, "x1": np.exp(lx1), "x2": np.exp(lx2)})

# --------------------------------------------------------------- 随机误差(H0)
def h0_baseline(seed, n, rel):
    """严格幂律 y = 1.5 * x1^0.7,高斯相对噪声。"""
    lx1, _, rng = _domain(n, seed)
    ly = np.log(1.5) + 0.7 * lx1
    y, _ = _finalize(ly, rng, rel)
    return _mkdf(y, lx1), "x1"

def h0_multivar(seed, n, rel):
    """y = 1.2 * x1^0.5 * x2^0.9,高斯相对噪声。"""
    lx1, lx2, rng = _domain(n, seed)
    ly = np.log(1.2) + 0.5 * lx1 + 0.9 * lx2
    y, _ = _finalize(ly, rng, rel)
    return _mkdf(y, lx1, lx2), "x1,x2"

def h0_corr(seed, n, rel):
    """y = 1.2 * x1^0.5 * x2^0.9,特征相关 rho=0.9(设计层 H0,系统误差为 0)。"""
    rng0 = np.random.default_rng(seed * 7919 + 104729)
    z1 = rng0.uniform(np.log(0.1), np.log(10.0), n)
    z2 = rng0.uniform(np.log(0.05), np.log(2.0), n)
    lz = np.log(0.5) + 0.5 * (z1 + z2)             # 相关化:对数空间线性混合
    lx1, lx2 = 2 * z1 - lz, 2 * z2 - lz
    rng = np.random.default_rng((seed + 31) * 15485863 + 7)
    ly = np.log(1.2) + 0.5 * lx1 + 0.9 * lx2
    y, _ = _finalize(ly, rng, rel)
    return _mkdf(y, lx1, lx2), "x1,x2"

SCENARIOS = {}
def _reg(name):
    def deco(fn):
        SCENARIOS[name] = fn
        return fn
    return deco

# ----------------------------------------------------------- 系统误差(H1)
@_reg("QUAD_LIN")
def quad_lin(seed, n, rel):
    """y = 1.5 * x1^0.7 * (1 + 0.10 * z1^2), z = x/x_ref 归一。
    单变量二次曲率:log 空间 ly = a + b*lx1 + c*(lx1-mean)^2, c = log(1+0.10)-类。"""
    lx1, _, rng = _domain(n, seed)
    z2 = (lx1 - lx1.mean()) ** 2
    ly = np.log(1.5) + 0.7 * lx1 + 0.10 * z2
    y, _ = _finalize(ly, rng, rel)
    return _mkdf(y, lx1), "x1"

@_reg("CUBIC_LIN")
def cubic_lin(seed, n, rel):
    """log 空间三次曲率:ly = a + b*lx1 + c3*(lx1-mean)^3。"""
    lx1, _, rng = _domain(n, seed)
    s = lx1.std()
    ly = np.log(1.5) + 0.7 * lx1 + 0.05 * ((lx1 - lx1.mean()) / s) ** 3
    y, _ = _finalize(ly, rng, rel)
    return _mkdf(y, lx1), "x1"

@_reg("INTERACT")
def interact(seed, n, rel):
    """y = 1.2 * x1^0.5 * x2^0.9 * (1 + 0.30 * z1 * z2), z=x/x_med。
    乘法交互项:log 空间产生 lx1*lx2 交叉项(两种 RESET 不建模,power 缺口)。"""
    lx1, lx2, rng = _domain(n, seed)
    cross = (lx1 - np.median(lx1)) * (lx2 - np.median(lx2))
    ly = np.log(1.2) + 0.5 * lx1 + 0.9 * lx2 + 0.30 * cross / cross.std()
    y, _ = _finalize(ly, rng, rel)
    return _mkdf(y, lx1, lx2), "x1,x2"

@_reg("THRESHOLD")
def threshold(seed, n, rel):
    """阈值结构:lx1 < m 时指数 0.4,>= m 时指数 1.0(分段幂律)。"""
    lx1, _, rng = _domain(n, seed)
    m = 0.5 * (np.log(0.1) + np.log(10.0))
    expo = np.where(lx1 < m, 0.4, 1.0)
    ly = np.log(1.5) + expo * lx1
    y, _ = _finalize(ly, rng, rel)
    return _mkdf(y, lx1), "x1"

@_reg("EXP_MODEL")
def exp_model(seed, n, rel):
    """指数:y = 1.0 * exp(0.5 * x1)(log 空间 = exp(x1),强凸)。"""
    lx1, _, rng = _domain(n, seed)
    ly = 0.5 * np.exp(lx1)
    y, _ = _finalize(ly, rng, rel)
    return _mkdf(y, lx1), "x1"

@_reg("LOG_MODEL")
def log_model(seed, n, rel):
    """对数:y = 0.8 + 1.1 * ln(x1),域平移 +2.4 使 x1 ∈ [1.1, 110](y 严格正)。"""
    lx1, _, rng = _domain(n, seed)
    lx1 = lx1 + 2.4
    ly = np.log(0.8 + 1.1 * lx1)
    y, _ = _finalize(ly, rng, rel)
    return _mkdf(y, lx1), "x1"

@_reg("RATIONAL")
def rational(seed, n, rel):
    """有理:y = x1 / (0.2 + 0.1 * x1)(饱和曲线,log 空间凹)。"""
    lx1, _, rng = _domain(n, seed)
    x1 = np.exp(lx1)
    ly = np.log(x1 / (0.2 + 0.1 * x1))
    y, _ = _finalize(ly, rng, rel)
    return _mkdf(y, lx1), "x1"

@_reg("OMITTED_M")
def omitted_m(seed, n, rel):
    """遗漏变量(H1):x2 与 x1 相关且二次进入真值,拟合面只见 x1。
    注:x2 若独立线性进入,E[ly|x1] 仍线性(截距吸收)——不构成形式失配,
    这是方法学要点,报告里明说。此处曲率强度与 QUAD_LIN 同量级(0.075)。"""
    lx1, _, rng = _domain(n, seed)
    lx2 = 0.5 + 0.4 * lx1 + 0.3 * rng.standard_normal(n)
    ly = np.log(1.5) + 0.7 * lx1 + 0.4 * lx2 ** 2
    y, _ = _finalize(ly, rng, rel)
    return pd.DataFrame({"target": y, "x1": np.exp(lx1)}), "x1"

@_reg("DILUTION")
def dilution(seed, n, rel):
    """一个变量曲率被稀释:x1 有曲率 c=0.10,x2 无 → form_check 的平均平方列
    携带的信号减半;检验它是否仍保有功效。"""
    lx1, lx2, rng = _domain(n, seed)
    z2 = (lx1 - lx1.mean()) ** 2
    ly = np.log(1.2) + 0.5 * lx1 + 0.9 * lx2 + 0.10 * z2
    y, _ = _finalize(ly, rng, rel)
    return _mkdf(y, lx1, lx2), "x1,x2"

@_reg("OPPOSING")
def opposing(seed, n, rel):
    """两个变量曲率方向相反:+0.12 与 -0.12。两特征同域独立(log[0.1,10])
    → (lx1,lx2) 可交换,平均曲率列 (z1+z2)/2 的总体系数精确为 0:
    form_check 的检验统计量期望携带零信号,而真实失配存在——对抗性场景。"""
    rng = np.random.default_rng(seed * 7919 + 104729)
    lx1 = rng.uniform(np.log(0.1), np.log(10.0), n)
    lx2 = rng.uniform(np.log(0.1), np.log(10.0), n)
    z1 = (lx1 - lx1.mean()) ** 2
    z2 = (lx2 - lx2.mean()) ** 2
    ly = np.log(1.2) + 0.5 * lx1 + 0.9 * lx2 + 0.12 * z1 - 0.12 * z2
    y, _ = _finalize(ly, rng, rel)
    return _mkdf(y, lx1, lx2), "x1,x2"

# --------------------------------------------------- 噪声层 H1(结构幂律 + 非高斯)
from scipy import stats as _stats


def _noise_layer(seed, n, kind):
    """返回 ly 尺度的噪声向量(均值 0),不改变均值结构。"""
    rng = np.random.default_rng(seed * 104729 + 13007)
    if kind == "t3":
        return _stats.t(3).rvs(n, random_state=rng) / np.sqrt(3.0)
    if kind == "var_growing":
        return None                                  # 特殊处理,见 H1_VAR_GROWING
    return rng.standard_normal(n)


@_reg("H1_HEAVYTAIL")
def h1_heavytail(seed, n, rel):
    """H1(噪声层):结构=严格幂律,噪声=Student-t(3) 重尾(标准化方差 1)。"""
    lx1, _, _rng = _domain(n, seed)
    eps = _noise_layer(seed, n, "t3")
    ly = np.log(1.5) + 0.7 * lx1 + rel * eps
    y = np.exp(ly)
    return _mkdf(y, lx1), "x1"

@_reg("H1_VAR_GROWING")
def h1_var_growing(seed, n, rel):
    """H1(噪声层):结构=严格幂律,SD 随 log 预测值线性增大(0.3x→1.7x rel)。"""
    lx1, _, rng = _domain(n, seed)
    ly_c = np.log(1.5) + 0.7 * lx1
    span = lx1.max() - lx1.min()
    sd_vec = rel * (0.3 + 1.4 * (lx1 - lx1.min()) / span)   # 单调于 lx1 → 随预测增大
    ly = ly_c + sd_vec * rng.standard_normal(n)
    y = np.exp(ly)
    return _mkdf(y, lx1), "x1"

@_reg("H1_AUTOCORR")
def h1_autocorr(seed, n, rel):
    """H1(噪声层):结构=严格幂律,AR(1) 自相关噪声 rho=0.6(平稳,标准化)。
    打破排列零分布的可交换性假设 → 测各法 FPR 失真。"""
    lx1, _, rng = _domain(n, seed)
    rho = 0.6
    eps_white = rng.standard_normal(n + 200)
    eps_ar = np.empty(n + 200)
    eps_ar[0] = eps_white[0] / np.sqrt(1 - rho ** 2)
    for i in range(1, n + 200):
        eps_ar[i] = rho * eps_ar[i - 1] + eps_white[i]
    eps_ar = eps_ar[200:] / np.sqrt(1.0 / (1 - rho ** 2))   # 弃 burn-in,单位方差
    ly = np.log(1.5) + 0.7 * lx1 + rel * eps_ar
    y = np.exp(ly)
    return _mkdf(y, lx1), "x1"

@_reg("H1_OUTLIER")
def h1_outlier(seed, n, rel):
    """H1(噪声层):严格幂律 + 2% 单点异常(ly 加 8σ 单侧冲击)。"""
    lx1, _, rng = _domain(n, seed)
    ly_c = np.log(1.5) + 0.7 * lx1
    eps = rng.standard_normal(n)
    k = max(1, int(0.02 * n))
    idx = rng.choice(n, size=k, replace=False)
    eps[idx] += 8.0
    ly = ly_c + rel * eps
    y = np.exp(ly)
    return _mkdf(y, lx1), "x1"

@_reg("H1_OUTLIER_GRP")
def h1_outlier_grp(seed, n, rel):
    """H1(噪声层):严格幂律 + 一组成组异常(连续 6% 区段,ly 平移 +5σ)。"""
    lx1, _, rng = _domain(n, seed)
    ly_c = np.log(1.5) + 0.7 * lx1
    eps = rng.standard_normal(n)
    k = max(2, int(0.06 * n))
    start = int(rng.integers(0, n - k))
    eps[start:start + k] += 5.0
    ly = ly_c + rel * eps
    y = np.exp(ly)
    return _mkdf(y, lx1), "x1"

@_reg("COLLINEAR_CURV")
def collinear_curv(seed, n, rel):
    """特征高度共线(rho≈0.95)+ x1 曲率:c=0.10。
    H1 双重压力:曲率信号 + 设计矩阵病态(测 rank/cond 失败登记与 p 值稳定性)。"""
    rng0 = np.random.default_rng(seed * 7919 + 104729)
    z1 = rng0.uniform(np.log(0.1), np.log(10.0), n)
    z2 = rng0.uniform(np.log(0.05), np.log(2.0), n)
    lz = np.log(0.5) + 0.925 * (z1 + z2)
    lx1, lx2 = 2 * z1 - lz, 2 * z2 - lz
    rng = np.random.default_rng((seed + 77) * 15485863 + 7)
    zq = (lx1 - lx1.mean()) ** 2
    ly = np.log(1.2) + 0.5 * lx1 + 0.9 * lx2 + 0.10 * zq
    y, _ = _finalize(ly, rng, rel)
    return _mkdf(y, lx1, lx2), "x1,x2"

# ---------------------------------------------------------------- 注册与元数据
SCENARIOS["H0_BASELINE"] = h0_baseline
SCENARIOS["H0_MULTIVAR"] = h0_multivar
SCENARIOS["H0_CORR"] = h0_corr

_M = (np.log(0.1) + np.log(10.0)) / 2.0      # x1 域中心(确定性)
_S = (np.log(10.0) - np.log(0.1)) / np.sqrt(12.0)

def _t_pl(X, ctx):        return np.log(1.5) + 0.7 * X[:, 0]
def _t_mv(X, ctx):        return np.log(1.2) + 0.5 * X[:, 0] + 0.9 * X[:, 1]
def _t_quad(X, ctx):      return _t_pl(X, ctx) + 0.10 * (X[:, 0] - ctx["m1"]) ** 2
def _t_cubic(X, ctx):     return _t_pl(X, ctx) + 0.05 * ((X[:, 0] - ctx["m1"]) / ctx["s1"]) ** 3
def _t_inter(X, ctx):
    cross = (X[:, 0] - ctx["med"][0]) * (X[:, 1] - ctx["med"][1])
    return _t_mv(X, ctx) + 0.30 * cross / ctx["cross_sd"]
def _t_thresh(X, ctx):
    expo = np.where(X[:, 0] < _M, 0.4, 1.0)   # 生成器用固定域中点,非样本均值
    return np.log(1.5) + expo * X[:, 0]
def _t_exp(X, ctx):       return 0.5 * np.exp(X[:, 0])
def _t_logm(X, ctx):      return np.log(0.8 + 1.1 * (X[:, 0] + 2.4))
def _t_rat(X, ctx):       return -np.log(0.2 + 0.1 * np.exp(X[:, 0])) + X[:, 0]
def _t_omit(X, ctx):
    lx2 = 0.5 + 0.4 * X[:, 0]
    return np.log(1.5) + 0.7 * X[:, 0] + 0.4 * lx2 ** 2
def _t_dil(X, ctx):       return _t_mv(X, ctx) + 0.10 * (X[:, 0] - ctx["m1"]) ** 2
def _t_opp(X, ctx):
    return (_t_mv(X, ctx) + 0.12 * (X[:, 0] - ctx["m1"]) ** 2
            - 0.12 * (X[:, 1] - ctx["m2"]) ** 2)

# name -> (类别 H0/H1, 族, 真值函数或 None, 说明)
SCENARIO_META = {
    "H0_BASELINE":    ("H0", "monomial",      _t_pl,   "y=1.5 x1^0.7 高斯"),
    "H0_MULTIVAR":    ("H0", "monomial",      _t_mv,   "y=1.2 x1^.5 x2^.9 高斯"),
    "H0_CORR":        ("H0", "monomial",      _t_mv,   "同上,rho=0.9 相关设计"),
    "H1_HEAVYTAIL":   ("H1", "noise-heavy",   _t_pl,   "结构幂律,t(3) 噪声"),
    "H1_VAR_GROWING": ("H1", "noise-hetero",  _t_pl,   "结构幂律,SD 随预测增长"),
    "H1_AUTOCORR":    ("H1", "noise-autocorr", _t_pl,  "结构幂律,AR(1) ρ=0.6 噪声"),
    "H1_OUTLIER":     ("H1", "noise-outlier", _t_pl,   "结构幂律,2% 单点 8σ"),
    "H1_OUTLIER_GRP": ("H1", "noise-outlier", _t_pl,   "结构幂律,6% 成组 +5σ"),
    "QUAD_LIN":       ("H1", "quadratic",     _t_quad,  "log 空间二次 c=+0.10"),
    "CUBIC_LIN":      ("H1", "cubic",         _t_cubic, "log 空间三次(标准化)0.05"),
    "INTERACT":       ("H1", "interaction",   _t_inter, "lx1*lx2 交叉项 0.30sd"),
    "THRESHOLD":      ("H1", "piecewise",     _t_thresh, "指数 0.4→1.0 阈值"),
    "EXP_MODEL":      ("H1", "exponential",   _t_exp,   "y=exp(0.5 x1)"),
    "LOG_MODEL":      ("H1", "logarithmic",   _t_logm,  "y=0.8+1.1 ln x1"),
    "RATIONAL":       ("H1", "rational",      _t_rat,   "y=x1/(0.2+0.1x1) 饱和"),
    "OMITTED_M":      ("H1", "omitted-var",   _t_omit,  "x2 相关且二次进入,被遗漏"),
    "DILUTION":       ("H1", "dilution",      _t_dil,   "x1 曲率 0.10,x2 无(平均稀释)"),
    "OPPOSING":       ("H1", "opposing",      _t_opp,   "+0.12 与 -0.12 反向曲率"),
    "COLLINEAR_CURV": ("H1", "collinear-quad", _t_dil,  "rho≈0.95 + x1 曲率 0.10"),
}

# LOG_MODEL 域平移:x1 = exp(lx+0.8),lx∈[ln0.1,ln10] → 真值函数按平移后 lx 表达
def _t_logm_fixed(X, ctx):
    return np.log(0.8 + 1.1 * (X[:, 0] + 0.8))
SCENARIO_META["LOG_MODEL"] = ("H1", "logarithmic", _t_logm_fixed, "y=0.8+1.1 ln x1")
