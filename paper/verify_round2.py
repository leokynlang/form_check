# -*- coding: utf-8 -*-
"""verify_round2.py — 第二轮审稿意见的数据核验。

P1. D(exog-HC3) 独立实现验证: 覆盖全部 5 个高条件数任务 + 其余任务分层抽样;
    三条路径对比(我们 X'X pinv / 我们 SVD / statsmodels),含统计量与 p 差异。
P2. .067 复核的精确计数与 Wilson 区间。
P3. S2 H1 侧配对任务级 bootstrap: Δ(D−A), Δ(D−B), Δ(D−C) 区间与 D 获胜比例;
    以及 gate_cmp 上 Δ(F−A)(F 无法与 D 配对,见正文说明)。
P4. 失效/良态任务的 cond(X_aug) 与 cond(X'X) 实测(供数值机制段引用)。

写入 paper/round2_results.json。
"""
import json
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats as st

warnings.filterwarnings("ignore")
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from srsd_data import registry_with_meta, make_dataset_fn  # noqa: E402
import common  # noqa: E402

OUT = {}

# ---------------------------------------------------------------- 公共小件
def svd_hc3_wald(y, Xlog, variant):
    """我们自己的 SVD 路径: hat 与三明治都经 np.linalg.pinv(X_u)(SVD)。"""
    n = y.shape[0]
    Xr = np.column_stack([np.ones(n), Xlog])
    yhat = Xr @ np.linalg.lstsq(Xr, y, rcond=None)[0]
    if variant == "fitted":
        added = np.column_stack([yhat ** 2, yhat ** 3])
    else:
        added = np.column_stack([Xlog ** 2, Xlog ** 3])
    Xu = np.column_stack([Xr, added])
    q = added.shape[1]
    b_u, _, _, _, _ = common._ols_svd(Xu, y) if hasattr(common, "_ols_svd") \
        else (np.linalg.lstsq(Xu, y, rcond=None)[0],) * 5
    P = np.linalg.pinv(Xu)                      # SVD 路径
    G = P @ P.T                                 # = (X'X)^-1, 经 X 的 SVD
    resid = y - Xu @ b_u
    hat = np.einsum("ij,jk,ik->i", Xu, G, Xu)
    eu = resid / (1.0 - hat)
    meat = Xu.T @ (Xu * (eu ** 2)[:, None])
    V = G @ meat @ G
    k = Xu.shape[1]
    Vq = V[np.ix_(range(k - q, k), range(k - q, k))]
    W = float(b_u[-q:] @ np.linalg.solve(Vq, b_u[-q:]))
    p = float(st.f.sf(max(W, 0.0) / q, q, n - k))
    return {"stat": W / q, "p": p, "Vq": Vq}


def sm_hc3_wald(y, Xlog, variant):
    """statsmodels 参考路径。"""
    n = y.shape[0]
    Xr = sm.add_constant(Xlog)
    yhat = Xr @ np.linalg.lstsq(Xr, y, rcond=None)[0]
    if variant == "fitted":
        added = np.column_stack([yhat ** 2, yhat ** 3])
    else:
        added = np.column_stack([Xlog ** 2, Xlog ** 3])
    Xu = sm.add_constant(added) if False else np.column_stack([Xr, added])
    m = sm.OLS(y, Xu).fit(cov_type="HC3")
    q = added.shape[1]
    R = np.zeros((q, Xu.shape[1]))
    R[:, -q:] = np.eye(q)
    w = m.wald_test(R, scalar=True, use_f=True)
    Vq = np.asarray(m.cov_params())[-q:, -q:]
    return {"stat": float(w.statistic), "p": float(w.pvalue), "Vq": Vq}


def cond_pair(Xlog, variant):
    n = Xlog.shape[0]
    Xr = np.column_stack([np.ones(n), Xlog])
    if variant == "fitted":
        yh = Xr @ np.linalg.lstsq(Xr, np.log1p(np.abs(Xlog).mean(
            axis=1)), rcond=None)[0]  # 占位: cond 与 y 无关,用设计本身
    aug_exog = np.column_stack([Xr, Xlog ** 2, Xlog ** 3])
    c_aug = float(np.linalg.cond(aug_exog))
    c_xtx = float(np.linalg.cond(aug_exog.T @ aug_exog))
    return c_aug, c_xtx


# ---------------------------------------------------------------- P1: D 验证
print("== P1: D(exog-HC3) 独立验证 ==")
HIGH = ["feynman-i.34.27", "feynman-i.38.12", "feynman-ii.27.18",
        "feynman-iii.12.43", "feynman-iii.15.14"]
plan = {}
for t in registry_with_meta():
    name = t.get("eq_name")
    if not name:
        continue
    if name in HIGH:
        plan[name] = [(200, 0.01, s) for s in range(20)]
    elif name.startswith("feynman-") and ".34." not in name:
        pass
# 其余 8 个 H0 任务分层: n×rel×seed
h0_all = [t.get("eq_name") for t in registry_with_meta()
          if t.get("klass", t.get("eq_class", "")) == "H0"
          or (t.get("eq_name") or "").startswith("feynman-")]
srsd = pd.read_parquet(os.path.join(ROOT, "runs", "srsd_full",
                                    "raw_srsd.parquet"))
h0_tasks = sorted(srsd[srsd["klass"] == "H0"]["task"].unique())
cls_by_task = {}
for t in registry_with_meta():
    if t.get("eq_name") in h0_tasks:
        cls_by_task[t["eq_name"]] = t["cls"]
strat = {}
for task in h0_tasks:
    if task in HIGH:
        strat[task] = [(200, 0.01, s) for s in range(20)]
    else:
        strat[task] = [(n, rel, s) for n in (50, 200)
                       for rel in (0.005, 0.03) for s in range(3)]

P1 = {}
cls_map = {}
for t in registry_with_meta():
    if t.get("eq_name") in strat:
        cls_map[t["eq_name"]] = t["cls"]
for task, grid in strat.items():
    fn = make_dataset_fn(cls_map[task])
    rej = {"xx": 0, "svd": 0, "sm": 0}
    dstat, dp, dcov = [], [], []
    conds = []
    n_tot = 0
    n_pair = 0
    for (n, rel, seed) in grid:
        df, feats = fn(seed=seed, n=n, rel=rel)
        fl = feats.split(",") if isinstance(feats, str) else list(feats)
        y = np.log(df["target"].values)
        Xlog = np.log(df[fl].values)
        row = common.run_all_methods(df, "target", fl, seed=seed)
        rej["xx"] += bool(row["rej_exog_hc3"])
        try:
            sv = svd_hc3_wald(y, Xlog, "exog")
            smr = sm_hc3_wald(y, Xlog, "exog")
            rej["svd"] += (sv["p"] < 0.05)
            rej["sm"] += (smr["p"] < 0.05)
            n_pair = n_pair + 1 if (np.isfinite(sv["stat"])
                                    and np.isfinite(smr["stat"])) else n_pair
            if (np.isfinite(sv["stat"]) and np.isfinite(smr["stat"])
                    and smr["stat"] > 0):
                dstat.append(abs(sv["stat"] - smr["stat"])
                             / max(smr["stat"], 1e-300))
                dp.append(abs(sv["p"] - smr["p"]))
                denom = np.linalg.norm(smr["Vq"])
                if denom > 0:
                    dcov.append(float(np.linalg.norm(sv["Vq"] - smr["Vq"])
                                      / denom))
        except Exception as exc:  # noqa: BLE001
            print('EXC', task, (n, rel, seed), repr(exc)[:120])
            dp.append(np.nan)
        n_tot += 1
        c_aug, c_xtx = cond_pair(Xlog, "exog")
        conds.append((c_aug, c_xtx))
    P1[task] = {
        "n_fits": n_tot,
        "FPR_our_xx_path": round(rej["xx"] / n_tot, 4),
        "FPR_our_svd_path": round(rej["svd"] / n_tot, 4),
        "FPR_statsmodels": round(rej["sm"] / n_tot, 4),
        "max_rel_stat_diff_svd_vs_sm": round(float(np.nanmax(dstat)), 4)
        if dstat else None,
        "max_abs_p_diff_svd_vs_sm": round(float(np.nanmax(dp)), 5)
        if dp else None,
        "max_cov_frobenius_rel_diff": round(float(np.nanmax(dcov)), 4)
        if dcov else None,
        "n_paired_fits": n_pair,
        "cond_X_median": round(float(np.median([c[0] for c in conds])), 2),
        "cond_XtX_median": round(float(np.median([c[1] for c in conds])), 2),
    }
    print(task, P1[task])
OUT["P1_D_verify"] = P1

# ---------------------------------------------------------------- P2: Wilson
print("== P2: .067 precision ==")
k, nn = 8, 120
ph = k / nn
z = 1.959963985
den = 1 + z ** 2 / nn
cen = (ph + z ** 2 / (2 * nn)) / den
half = z * np.sqrt(ph * (1 - ph) / nn + z ** 2 / (4 * nn ** 2)) / den
OUT["P2_wilson"] = {"k": k, "n": nn, "phat": round(ph, 4),
                    "wilson95": [round(cen - half, 4), round(cen + half, 4)]}
print(OUT["P2_wilson"])

# ---------------------------------------------------------------- P3: H1 配对 bootstrap
print("== P3: paired task bootstrap ==")
s = srsd
h0 = s[(s["klass"] == "H0") & (s["rel"] >= 0.005)]
h1 = s[(s["klass"] == "H1") & (s["rel"] >= 0.005)]
METHODS = {"A": "p_fc", "B": "p_fitted_f", "C": "p_fitted_hc3",
           "D": "p_exog_hc3"}
thr = {}
for m, col in METHODS.items():
    for n in (50, 200, 1000):
        pool = h0[h0["n"] == n][col].dropna().values
        thr[(m, n)] = float(np.quantile(pool, 0.05))
tasks = sorted(h1["task"].unique())
P3 = {"n_tasks": len(tasks), "deltas": {}, "win_fraction": {}}
for a, b in [("D", "A"), ("D", "B"), ("D", "C")]:
    per = {}
    for n in (50, 200, 1000):
        sel = h1[h1["n"] == n]
        pa, pb = [], []
        for t in tasks:
            g = sel[sel["task"] == t]
            pa.append(float((g[METHODS[a]] <= thr[(a, n)]).mean()))
            pb.append(float((g[METHODS[b]] <= thr[(b, n)]).mean()))
        pa, pb = np.array(pa), np.array(pb)
        rng = np.random.default_rng(20260908)
        boots = [float((pa[pick] - pb[pick]).mean())
                 for pick in (rng.integers(0, len(tasks), len(tasks))
                              for _ in range(10000))]
        per[f"n{n}"] = {
            "mean_D": round(float(pa.mean()), 3),
            "mean_other": round(float(pb.mean()), 3),
            "delta": round(float((pa - pb).mean()), 3),
            "ci95": [round(float(np.quantile(boots, 0.025)), 3),
                     round(float(np.quantile(boots, 0.975)), 3)],
            "D_wins_tasks": round(float((pa > pb).mean()), 3)}
    P3["deltas"][f"{a}-{b}"] = per
# F vs A (gate_cmp, n=50/200)
g = pd.read_csv(os.path.join(ROOT, "runs", "gate_cmp", "results.csv"))
gs = g[(g["bench"] == "srsd") & (g["klass"] == "H1")
       & (~g["group"].str.startswith("noise"))]
def calib(arr):
    arr = np.asarray([v for v in arr if np.isfinite(v)])
    cand = np.unique(np.concatenate([arr, [arr.max() + 1e-12]]))
    for t in cand:
        if float((arr > t).mean()) <= 0.05:
            return float(t)
    return float(arr.max())
gh0 = g[(g["bench"] == "srsd") & (g["klass"] == "H0")]
per = {}
for n in (50, 200):
    m_star = calib(gh0[gh0["n"] == n]["delta_aicc"].tolist())
    sel = gs[gs["n"] == n]
    fpa, apa = [], []
    for t in sorted(sel["task"].unique()):
        gt = sel[sel["task"] == t]
        fpa.append(float((gt["delta_aicc"] > m_star).mean()))
        apa.append(float(gt["rej_fc"].mean()))
    fpa, apa = np.array(fpa), np.array(apa)
    rng = np.random.default_rng(20260908)
    boots = [float((fpa[pick] - apa[pick]).mean())
             for pick in (rng.integers(0, len(fpa), len(fpa))
                          for _ in range(10000))]
    per[f"n{n}"] = {"mean_F": round(float(fpa.mean()), 3),
                    "mean_A": round(float(apa.mean()), 3),
                    "delta": round(float((fpa - apa).mean()), 3),
                    "ci95": [round(float(np.quantile(boots, 0.025)), 3),
                             round(float(np.quantile(boots, 0.975)), 3)],
                    "F_wins_tasks": round(float((fpa > apa).mean()), 3)}
P3["deltas"]["F-A"] = per
OUT["P3_paired"] = P3
print(json.dumps(P3["deltas"]["D-A"], indent=1))
print(json.dumps(P3["deltas"]["F-A"], indent=1))

with open(os.path.join(HERE, "round2_results.json"), "w",
          encoding="utf-8") as fh:
    json.dump(OUT, fh, indent=1, ensure_ascii=False)
print("saved -> paper/round2_results.json")
