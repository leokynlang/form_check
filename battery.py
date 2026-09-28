# -*- coding: utf-8 -*-
"""battery.py — 多函数族"形式合格性电池"正式版(实验层,核心代码零改动)。

[机制] 每个可线性化函数族 F 一组坐标(u(x), v(y)),F 的检验 H0 =
"条件均值在 F 坐标里线性"。统计量 = [1,u,u²] 设计下二次列系数;
零分布 = wild bootstrap(Rademacher 翻转残差)。

[为什么 wild 不是排列] 探针 1(runs/battery_probe,排列)实锤:
v(y)=y 与 v(y)=1/y 坐标下乘性噪声残差异方差,排列可交换性破裂 →
rational 假升级 30-34%。wild 修复到 4.5-9%(runs/battery_probe_wild)。
power/exponential 族(log 变换同方差化)两法等价——坐标决定零分布合法性,
不是实现细节。

[语义] 输出候选列表 viable={F: p_F>=α} 而非单一赢家(项目宪法)。
三种错误分开报告,不合并:
- 假排除(真族∉viable / 空集假升级):错过可用形式,代价=升级到贵路径;
- 假包含(错族∈viable):给出不一致候选,代价=下游选错;
- 电池无该族(truth_family=none):升级是正确结局,escalation=正确负。
合成 OPPOSING 已知 power 族会假包含(fc 同源盲区),正式跑量化它。

[多重性] 每族各自 α=0.05,不做族间校正:viable 是"各自检验不拒绝"的
可行集,非拒绝语义下多重性不膨胀第一类错误。max-T 只属于"选最小 p
当路由"的用法,那不在本电池语义内。
"""
import numpy as np

B_WILD = 999
ALPHA = 0.05

# 族注册:u 作用在每个特征列上(逐列),v 作用在 y 上。
FAMILIES = {
    "power":       (np.log,   np.log),
    "exponential": (lambda x: x, np.log),
    "logarithmic": (np.log,   lambda y: y),
    "rational":    (lambda x: 1.0 / x, lambda y: 1.0 / y),
}
FAM_NAMES = list(FAMILIES)

# 合成场景真值族注册(先于正式跑写死)。
# power*   = 严格单项式(结构),含噪声层场景(检验的是噪声鲁棒性);
# quad/cubic/dilution/collinear/omitted/threshold/interact/opposing
#          = 无电池族(线性化加性形式装不下)→ 升级是正确结局,
#            同时量化"假包含"(错族留在 viable)的族级分布。
SYNTHETIC_TRUTH = {
    "H0_BASELINE": "power", "H0_MULTIVAR": "power", "H0_CORR": "power",
    "H1_HEAVYTAIL": "power", "H1_VAR_GROWING": "power",
    "H1_AUTOCORR": "power", "H1_OUTLIER": "power", "H1_OUTLIER_GRP": "power",
    "EXP_MODEL": "exponential", "LOG_MODEL": "logarithmic",
    "RATIONAL": "rational",
    "QUAD_LIN": None, "CUBIC_LIN": None, "DILUTION": None,
    "COLLINEAR_CURV": None, "OMITTED_M": None, "THRESHOLD": None,
    "INTERACT": None, "OPPOSING": None,
}


def _standardize_cols(U):
    mu = U.mean(axis=0)
    sd = U.std(axis=0)
    sd = np.where(sd > 0, sd, 1.0)
    return (U - mu) / sd


def spec_test_wild(X, y, uf, vf, seed, b_wild=B_WILD):
    """单族形式合格性检验(多特征):u=uf(X) 逐列,quad=标准化 u 列平方的行均值。

    设计 [1, U_std, quad];统计量=二次列系数;wild 零分布 =
    pred + resid*Rademacher,经 pinv 末行向量化(c2 是 v 的线性泛函)。
    失败显式登记,不抛异常,不出垃圾 p 值。
    """
    out = {"p": np.nan, "t": np.nan, "fail": None}
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    if (not np.isfinite(X).all()) or (not np.isfinite(y).all()):
        out["fail"] = "nonfinite"
        return out
    if (X <= 0).any() or (y <= 0).any():
        out["fail"] = "nonpositive"
        return out
    try:
        U_raw = uf(X)                          # 逐列变换(n×p)
        v = vf(y)
    except (ValueError, FloatingPointError, OverflowError):
        out["fail"] = "transform_failed"
        return out
    if (not np.isfinite(U_raw).all()) or (not np.isfinite(v).all()):
        out["fail"] = "transform_nonfinite"
        return out
    U = _standardize_cols(U_raw)
    n = len(v)
    base = np.column_stack([np.ones(n), U])
    quad = (U ** 2).mean(axis=1)
    design = np.column_stack([base, quad])
    coef, _, _, _ = np.linalg.lstsq(design, v, rcond=None)
    t_obs = float(coef[-1])
    b, _, _, _ = np.linalg.lstsq(base, v, rcond=None)
    pred = base @ b
    resid = v - pred
    if float(resid @ resid) <= 1e-12 * float(v @ v):
        out["fail"] = "degenerate_zero_residual"
        return out
    rng = np.random.default_rng(seed)
    pinv_row = np.linalg.pinv(design)[-1, :]
    c2_base = float(pinv_row @ pred)
    signs = rng.choice((-1.0, 1.0), size=(b_wild, n))
    null = c2_base + (resid[None, :] * signs) @ pinv_row
    mu, sd = float(null.mean()), float(null.std())
    if not np.isfinite(sd) or sd <= 0:
        out["fail"] = "degenerate_null"
        return out
    p = (1.0 + float(np.sum(np.abs(null - mu) >= abs(t_obs - mu)))) / (b_wild + 1.0)
    out.update({"p": float(p), "t": float((t_obs - mu) / sd)})
    return out


def battery_run(df, target, features, seed, b_wild=B_WILD):
    """对一个数据集跑全部四族,返回宽行(p/rej/fail_族 + viable 字符串)。"""
    X = df[list(features)].values.astype(float)
    y = df[target].values.astype(float)
    row = {}
    viable = []
    for i, fam in enumerate(FAM_NAMES):
        uf, vf = FAMILIES[fam]
        bs = int(seed) * 1_000_003 + i * 7919 + 17
        r = spec_test_wild(X, y, uf, vf, bs, b_wild=b_wild)
        row[f"p_{fam}"] = r["p"]
        row[f"rej_{fam}"] = bool(np.isfinite(r["p"]) and r["p"] < ALPHA)
        row[f"fail_{fam}"] = r["fail"]
        if (not r["fail"]) and r["p"] >= ALPHA:
            viable.append(fam)
    row["viable"] = ",".join(viable)
    row["n_viable"] = len(viable)
    return row
