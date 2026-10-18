"""common.py — form_check vs Ramsey RESET 公共实验台(只读 dimreflex)。

纪律(不可违背):
- 方法 A 调用 discover.py::form_check 本体,不复制不暗改;其 |t|>3 原生阈值
  仅作"原生报警率"留痕,一切决策用包装层算的经验 p 值(同一排列零分布),
  四种方法 thus 在同一 alpha 尺度上比较(修掉旧 pilot 直接比 |t|>3 vs p<0.05 的不公平)。
- 任何失败(秩亏/非正值/样本量/数值/统计量失配/完美拟合退化)都写进行内
  failure 字段,绝不静默丢弃。
- 本目录外零写入。
"""
import hashlib
import os
import platform
import subprocess
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
# 2026-09-05 拆分:vendored 包在本目录 ./dimreflex/;_bootstrap 强制其优先
REPO = HERE
PKG_DIR = os.path.join(HERE, "dimreflex")
import _bootstrap  # noqa: E402

from dimreflex import form_check as dimreflex_form_check
from dimreflex._ols import ols as _dimreflex_ols

ALPHA = 0.05          # 所有方法统一的显著性水平
B_PERM = 999          # 包装层排列零分布次数(p 值分辨率 1/1000)
_WMARK = object()     # 占位,防误用

# validate_friction.py 已移至 validation/(2026-09-04 目录整理),路径同步更新;
# 历史快照中的哈希语义不变(同一文件内容)。
_CORE_FILES = ["dimreflex/discover.py", "dimreflex/_ols.py",
               "validation/validate_friction.py", "dimreflex/__init__.py"]


# ---------------------------------------------------------------- 环境快照
def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def capture_env():
    import dimreflex
    import scipy
    meta = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "pandas": pd.__version__,
        "dimreflex_version": getattr(dimreflex, "__version__", "unknown"),
        "platform": f"{platform.system()} {platform.release()} {platform.machine()}",
        "core_file_sha256": {},
        "alpha": ALPHA,
        "b_perm": B_PERM,
        "form_check_signature": (
            "discover.py::form_check(data, target, features, n_boot, random_state)"
            " — 调用本体,未修改"),
    }
    for name in _CORE_FILES:
        p = os.path.join(REPO, name)
        meta["core_file_sha256"][name] = sha256_file(p) if os.path.exists(p) else "MISSING"
    try:
        head = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=REPO,
            capture_output=True, text=True, timeout=15, check=False)
        meta["git_head"] = head.stdout.strip()
        dirty = subprocess.run(
            ["git", "status", "--porcelain"], cwd=REPO,
            capture_output=True, text=True, timeout=15, check=False)
        meta["git_dirty_entries"] = len(dirty.stdout.strip().splitlines())
    except Exception as exc:  # noqa: BLE001 — git 缺失/超时统一记入审计元数据,不中断实验
        meta["git_head"] = f"unavailable: {exc}"
    meta["scenarios_sha256"] = sha256_file(os.path.join(HERE, "scenarios.py"))
    meta["common_sha256"] = sha256_file(os.path.abspath(__file__))
    return meta


# ---------------------------------------------------------------- OLS 小件
def _ols(X, y):
    """最小二乘:返回 coef, resid, rss, hat 对角线。X 需已含截距列。"""
    coef, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ coef
    rss = float(resid @ resid)
    XtX_inv = np.linalg.pinv(X.T @ X)
    hat = np.einsum("ij,jk,ik->i", X, XtX_inv, X)
    return coef, resid, rss, hat, XtX_inv


def _classify_fc_reason(reason):
    if "非有限" in reason:
        return "checked_false_nonfinite"
    if "非正值" in reason:
        return "checked_false_nonpositive"
    if "过少" in reason:
        return "checked_false_small_n"
    return "checked_false_other:" + reason[:60]


# ---------------------------------------------------------------- 方法 A
def form_check_with_p(data, target, features, seed, n_boot=200, b_perm=B_PERM):
    """调用 dimreflex.form_check 本体,再按其排列统计量算经验 p 值。

    统计量 = 偏残差曲率 c2(form_check 返回的 'curvature')。
    零分布 = 线性拟合残差排列重建伪数据集(form_check 文档记载的方案),
    由包装层以独立流实现,次数 B_PERM(分辨率 1/(B_PERM+1))。
    p_emp = (1 + #{|c2_perm - mean| >= |c2_obs - mean|}) / (B_PERM + 1)。

    一致性自检:包装层用同一设计重算 c2_obs,与 form_check 返回的 curvature
    比对,不符 → failure='statistic_mismatch'(意味着 discover.py 实现漂移,
    排列零分布失效),该数据集记失败,不 silently 出 p 值。
    """
    t0 = time.perf_counter()
    fc = dimreflex_form_check(data, target=target, features=list(features),
                              n_boot=n_boot, random_state=seed)
    out = {"p_fc": np.nan, "rej_fc": False, "fc_t_native": fc.get("t_stat", np.nan),
           "fc_alarm_native": bool(fc.get("alarm", False)),
           "curvature": fc.get("curvature", np.nan),
           "fail_fc": None, "time_fc": time.perf_counter() - t0}
    if not fc.get("checked", False):
        out["fail_fc"] = _classify_fc_reason(str(fc.get("reason", "")))
        return out

    y = data[target].values
    Xraw = data[list(features)].values
    lx = np.log(Xraw)
    ly = np.log(y)
    n = len(ly)
    if len(features) == 1:
        quad = (lx[:, 0] - lx[:, 0].mean()) ** 2
    else:
        quad = ((lx - lx.mean(axis=0)) ** 2).mean(axis=1)
    base = np.column_stack([np.ones(n), lx])
    design = np.column_stack([base, quad])

    coef, _, _, _, _ = _ols(design, ly)
    c2_wrap = float(coef[-1])
    if not np.isclose(c2_wrap, fc["curvature"], rtol=1e-6, atol=1e-10):
        out["fail_fc"] = "statistic_mismatch"
        return out

    # 排列零分布:线性拟合的预测 + 打乱的残差(与 form_check 文档一致)。
    # c2 是 y 的线性泛函:c2(pred+perm) = c2(pred) + pinv行 @ perm,
    # 故用 pinv 最后一行把 B_PERM 次排列向量化(与逐次 lstsq 差 ~1e-13 相对)。
    bcoef, _, _, _, _ = _ols(base, ly)
    pred = base @ bcoef
    resid = ly - pred
    ly_ss = float(ly @ ly)
    if float(resid @ resid) <= 1e-12 * ly_ss:
        out["fail_fc"] = "degenerate_zero_residual"   # 零噪声+精确幂律:检验不可运行
        return out
    rng = np.random.default_rng(int(seed) * 1_000_003 + 17)
    pinv_row = np.linalg.pinv(design)[-1, :]          # (n,): c2 = pinv_row @ y
    c2_base = float(pinv_row @ pred)
    idx = np.argsort(rng.random((b_perm, n)), axis=1)
    null_c2 = c2_base + resid[idx] @ pinv_row
    mu, sd = float(null_c2.mean()), float(null_c2.std())
    if not np.isfinite(sd) or sd <= 0:
        out["fail_fc"] = "degenerate_zero_residual"
        return out
    t_obs = (c2_wrap - mu) / sd
    p_emp = (1.0 + float(np.sum(np.abs(null_c2 - mu) >= abs(c2_wrap - mu)))) / (b_perm + 1.0)
    out.update({"p_fc": float(p_emp), "rej_fc": bool(p_emp < ALPHA),
                "fc_t_native": float(fc["t_stat"]),
                "fc_t_emp": float(t_obs), "curvature": c2_wrap})
    return out


# ---------------------------------------------------------------- 方法 B/C/D
def reset_test(Xlog, ly, powers=(2, 3), variant="fitted", vcov="classical"):
    """Ramsey RESET(log 空间)。Xlog: n×p 的 log 特征矩阵(不含截距)。

    variant='fitted': 加入 yhat_r^k(k ∈ powers) — Ramsey(1969) 教科书形式。
    variant='exog':   加入每个特征列的 k 次幂 — exog 变体。
    vcov='classical': 经典 F = [(RSS_r-RSS_u)/q] / [RSS_u/(n-k)]。
    vcov='hc3':       HC3 夹心协方差的 Wald 统计量 W,参考分布 F(q, n-k)
                      (W/q, 与 car::linearHypothesis 的稳健检验同口径)。
    返回 dict(p, rej, stat, q, k, fail, cond)。任何失败写 fail,不抛异常。
    """
    n = ly.shape[0]
    X_r = np.column_stack([np.ones(n), Xlog])
    p_feat = X_r.shape[1]
    ly_ss = float(ly @ ly)
    try:
        if np.linalg.matrix_rank(X_r) < p_feat:
            return {"p": np.nan, "rej": False, "stat": np.nan, "q": 0,
                    "k": 0, "fail": "rank_deficient", "cond": np.inf}
        b_r, _, rss_r, _, _ = _ols(X_r, ly)
        yhat = X_r @ b_r
        if variant == "fitted":
            added = np.column_stack([yhat ** k for k in powers])
        else:
            added = np.column_stack([Xlog[:, j] ** k
                                     for j in range(Xlog.shape[1]) for k in powers])
        X_u = np.column_stack([X_r, added])
        q = added.shape[1]
        k = X_u.shape[1]
        if n - k < 1:
            return {"p": np.nan, "rej": False, "stat": np.nan, "q": q, "k": k,
                    "fail": "insufficient_df", "cond": np.inf}
        if np.linalg.matrix_rank(X_u) < k:
            return {"p": np.nan, "rej": False, "stat": np.nan, "q": q, "k": k,
                    "fail": "rank_deficient", "cond": np.inf}
        b_u, resid_u, rss_u, hat_u, XtX_inv = _ols(X_u, ly)
        cond = float(np.linalg.cond(X_u))
        if float(resid_u @ resid_u) <= 1e-12 * ly_ss:
            return {"p": np.nan, "rej": False, "stat": np.nan, "q": q, "k": k,
                    "fail": "degenerate_zero_residual", "cond": cond}
        if vcov == "classical":
            if not (rss_u > 0):
                return {"p": np.nan, "rej": False, "stat": np.nan, "q": q, "k": k,
                        "fail": "degenerate_perfect_fit", "cond": cond}
            F = ((rss_r - rss_u) / q) / (rss_u / (n - k))
            F = max(float(F), 0.0)
            if not np.isfinite(F):
                return {"p": np.nan, "rej": False, "stat": np.nan, "q": q, "k": k,
                        "fail": "numerical_nonfinite", "cond": cond}
            pval = float(stats.f.sf(F, q, n - k))
            return {"p": pval, "rej": bool(pval < ALPHA), "stat": F, "q": q,
                    "k": k, "fail": None, "cond": cond}
        # HC3: V = (X'X)^-1 [Σ (e_i/(1-h_i))^2 x_i x_i'] (X'X)^-1
        one_m_h = 1.0 - hat_u
        if np.any(one_m_h <= 1e-12):
            return {"p": np.nan, "rej": False, "stat": np.nan, "q": q, "k": k,
                    "fail": "hc3_leverage_degenerate", "cond": cond}
        eu = resid_u / one_m_h
        meat = X_u.T @ (X_u * (eu ** 2)[:, None])
        V = XtX_inv @ meat @ XtX_inv
        R = np.zeros((q, k))
        R[:, -q:] = np.eye(q)
        Vq = V[np.ix_(range(k - q, k), range(k - q, k))]
        try:
            W = float(b_u[-q:] @ np.linalg.solve(Vq, b_u[-q:]))
        except np.linalg.LinAlgError:
            return {"p": np.nan, "rej": False, "stat": np.nan, "q": q, "k": k,
                    "fail": "hc3_singular_covariance", "cond": cond}
        W = max(W, 0.0)
        if not np.isfinite(W):
            return {"p": np.nan, "rej": False, "stat": np.nan, "q": q, "k": k,
                    "fail": "numerical_nonfinite", "cond": cond}
        pval = float(stats.f.sf(W / q, q, n - k))
        return {"p": pval, "rej": bool(pval < ALPHA), "stat": W / q, "q": q,
                "k": k, "fail": None, "cond": cond}
    except np.linalg.LinAlgError:
        return {"p": np.nan, "rej": False, "stat": np.nan, "q": 0, "k": 0,
                "fail": "numerical_linalg", "cond": np.inf}


_RESET_TAGS = [("reset_fitted_F", {"variant": "fitted", "vcov": "classical"}),
               ("reset_fitted_HC3", {"variant": "fitted", "vcov": "hc3"}),
               ("reset_exog_HC3", {"variant": "exog", "vcov": "hc3"})]


def run_all_methods(df, target, features, seed, n_boot=200):
    """对一个数据集跑 A+B+C+D,返回宽行 dict(失败字段逐方法留痕)。"""
    Xlog = np.log(df[list(features)].values)
    ly = np.log(df[target].values)
    row = {}
    row.update(form_check_with_p(df, target, features, seed, n_boot=n_boot))
    for name, kw in _RESET_TAGS:
        t0 = time.perf_counter()
        r = reset_test(Xlog, ly, powers=(2, 3), **kw)
        dt = time.perf_counter() - t0
        key = name.split("_", 1)[1].lower()          # fitted_f / fitted_hc3 / exog_hc3
        row[f"p_{key}"] = r["p"]
        row[f"rej_{key}"] = r["rej"]
        row[f"stat_{key}"] = r["stat"]
        row[f"fail_{key}"] = r["fail"]
        row[f"cond_{key}"] = r["cond"]
        row[f"time_{key}"] = dt
    return row


# ---------------------------------------------------------------- 拟合优度
def fitting_quality(Xlog, ly):
    """受限(基)幂律 OLS 的拟合优度 + RESET 统计量中间量。

    与 A/B/C/D 四门无关,只描述"受限幂律拟合本身多好"——即
    fitted-but-wrong 现象学的自变量侧:拟合优度高不代表形式正确。

    指标(log 空间):
    - r2 / r2_adj:决定系数与自由度修正;
    - aic / bic:以 n·ln(RSS/n)+k·penalty 的高斯似然口径(常数项略去,
      仅用于跨模型/跨场景相对比较);
    - rmse_log:残差均方根;sigma2:无偏残差方差(引擎给出);
    - quality_cond:受限设计条件数;n_par:参数数(含截距)。

    RESET 中间量(fitted-F 的机制侧):
    - yhat_sd / resid_sd:受限预测与残差的尺度;
    - aug_r2_yhat2 / aug_r2_yhat3:把 (yhat-mean)^k 投影到基设计上的 R²(k=2,3)
      —— 幂次增广列已被基张成的份额。越接近 1,说明该增广列相对基近乎冗余,
      fitted-RESET 的该方向携带的独立信号越少(解释近共线盲区的直接量);
    - corr_yhat2_resid / corr_yhat3_resid:|corr(yhat^k, 残差)|,残差中对应
      幂次弯折与拟合值的线性关联强度 —— fitted-RESET 能否借该列检出。

    失败(秩亏/退化)写 fail_quality,不抛异常;与 reset_test 口径一致。
    调用 dimreflex._ols.ols 单一引擎,不复制统计数学。
    """
    Xlog = np.asarray(Xlog, dtype=float)
    ly = np.asarray(ly, dtype=float)
    n = ly.shape[0]
    design = np.column_stack([np.ones(n), Xlog])
    p = int(design.shape[1])
    out = {"r2": np.nan, "r2_adj": np.nan, "aic": np.nan, "bic": np.nan,
           "rmse_log": np.nan, "sigma2": np.nan, "quality_cond": np.nan,
           "n_par": p, "yhat_sd": np.nan, "resid_sd": np.nan,
           "aug_r2_yhat2": np.nan, "aug_r2_yhat3": np.nan,
           "corr_yhat2_resid": np.nan, "corr_yhat3_resid": np.nan,
           "fail_quality": None}
    try:
        if np.linalg.matrix_rank(design) < p:
            out["fail_quality"] = "rank_deficient"
            return out
        d = _dimreflex_ols(design, ly)
        rss, tss = float(d["ss_res"]), float(d["ss_tot"])
        out["r2"] = float(d["r_squared"])
        out["sigma2"] = float(d["sigma2"])
        out["quality_cond"] = float(d["condition_number"])
        if not np.isfinite(rss) or not np.isfinite(tss) or tss <= 0:
            out["fail_quality"] = "degenerate_ss"
            return out
        if rss <= 1e-12 * tss:
            out["fail_quality"] = "degenerate_zero_residual"
            return out
        out["r2_adj"] = (1.0 - (1.0 - out["r2"]) * (n - 1) / (n - p)
                         if n > p else np.nan)
        out["rmse_log"] = float(np.sqrt(rss / n))
        out["aic"] = float(n * np.log(rss / n) + 2 * p)
        out["bic"] = float(n * np.log(rss / n) + p * np.log(n))
        yhat = np.asarray(d["prediction_log"], dtype=float)
        resid = ly - yhat
        yhat_sd = float(yhat.std())
        out["yhat_sd"] = yhat_sd
        out["resid_sd"] = float(resid.std())
        yc = yhat - yhat.mean()
        for k, aug_key, corr_key in ((2, "aug_r2_yhat2", "corr_yhat2_resid"),
                                     (3, "aug_r2_yhat3", "corr_yhat3_resid")):
            col = yc ** k
            ck, _, _, _ = np.linalg.lstsq(design, col, rcond=None)
            ss_res_k = float(((col - design @ ck) ** 2).sum())
            ss_tot_k = float(((col - col.mean()) ** 2).sum())
            out[aug_key] = (1.0 - ss_res_k / ss_tot_k) if ss_tot_k > 0 else np.nan
            if resid.std() > 0 and col.std() > 0:
                out[corr_key] = float(abs(np.corrcoef(col, resid)[0, 1]))
    except np.linalg.LinAlgError:
        out["fail_quality"] = "numerical_linalg"
    return out


# ---------------------------------------------------------------- 外推误差
def extrapolation_error(df, target, features, truth_logmean, far=3.0, m_pts=12):
    """训练域外的真值 vs 幂律预测(log 域外推 ×frac 倍量程)。

    truth_logmean: callable(lx_matrix)->log 空间真均值(各场景注册的真函数)。
    返回 (中位相对误差, 最远点相对误差),均为相对真值的 |pred-true|/true。
    """
    Xlog_tr = np.log(df[list(features)].values)
    lo, hi = Xlog_tr.min(axis=0), Xlog_tr.max(axis=0)      # 已是 log 空间
    m_pts_lo = m_pts // 2
    dlog = float(np.log(far))
    grid_lo = np.linspace(lo[0] - dlog, lo[0], m_pts_lo)   # 低端外推(log 平移)
    grid_hi = np.linspace(hi[0], hi[0] + dlog,
                          m_pts - m_pts_lo)                # 高端外推
    lx1 = np.concatenate([grid_lo, grid_hi])
    Xe = np.tile(np.median(Xlog_tr, axis=0), (len(lx1), 1))
    Xe[:, 0] = lx1
    ye_true = truth_logmean(Xe)
    Xtr = np.column_stack([np.ones(len(Xlog_tr)), Xlog_tr])
    b, _, _, _, _ = _ols(Xtr, np.log(df[target].values))       # 训练域内拟合
    ye_pred = np.column_stack([np.ones(len(Xe)), Xe]) @ b      # 外推点上预测
    rel = np.abs(np.exp(ye_pred - ye_true) - 1.0)   # log 差稳定形式(防溢出)
    return float(np.median(rel)), float(rel[-1])
