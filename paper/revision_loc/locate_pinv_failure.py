"""locate_pinv_failure.py — HC3 失效归因定位(审稿意见第4条的实验核实)。

三个实验:
  E1 复现:失败任务(i.34.27, iii.12.43)+健康对照(i.12.1),n=200, rel=1%,
     60 seeds,同一数据喂四条 HC3 路径:
       harness(pinv(X'X) 夹心,Gram hat)
       statsmodels(SVD 路径)
       svd-native(自写:SVD hat + pinv(X)pinv(X)^T 协方差)
       classical-F(参照,应校准)
  E2 载体分解:i.34.27 上 20 seeds,hat 来源 × 协方差来源 2×2 交叉,
     定位失效在杠杆 hat、夹心协方差、还是两者。
  E3 条件数扫描:全部 13 个 H0 任务各 10 seeds,记录 cond(X_aug)、
     cond(X'X)、pinv 截断偏差、log10 采样域宽,harness-C 拒绝率对条件数。

纪律:只读 run 数据生成器;一切结果写本目录 results_localization.json。
"""
import json
import os
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, r"D:\form_check")

from scipy import stats as sstats  # noqa: E402
import statsmodels.api as smapi    # noqa: E402

from common import _ols, reset_test            # noqa: E402  harness 本体路径
from srsd_data import make_dataset_fn, registry_with_meta  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ALPHA = 0.05
N, REL = 200, 0.01


# ---------------------------------------------------------------- SVD-native 路径
def svd_ols(X, y):
    """coef/resid/hat/XtX_inv 全部从 X 本身的 SVD 导出,不组 Gram 矩阵。"""
    U, s, Vt = np.linalg.svd(X, full_matrices=False)
    s_safe = np.where(s > 0, s, 1.0)
    coef = Vt.T @ ((U.T @ y) / s_safe)
    resid = y - X @ coef
    hat = np.einsum("ij,ij->i", U, U)              # H = U U^T 的对角
    XtX_inv = (Vt.T / (s_safe ** 2)) @ Vt          # = pinv(X) @ pinv(X)^T
    return coef, resid, hat, XtX_inv


def hc3_wald_from(X, coef, resid, hat, XtX_inv, q):
    """给定 hat/XtX_inv 来源的 HC3-Wald(q 个新增列的联合检验)。"""
    n, k = X.shape
    one_m_h = 1.0 - hat
    if np.any(one_m_h <= 1e-12):
        return np.nan, "lev_degen"
    eu = resid / one_m_h
    meat = X.T @ (X * (eu ** 2)[:, None])
    V = XtX_inv @ meat @ XtX_inv
    Vq = V[np.ix_(range(k - q, k), range(k - q, k))]
    try:
        W = float(coef[-q:] @ np.linalg.solve(Vq, coef[-q:]))
    except np.linalg.LinAlgError:
        return np.nan, "singular"
    W = max(W, 0.0)
    if not np.isfinite(W):
        return np.nan, "nonfinite"
    return float(sstats.f.sf(W / q, q, n - k)), None


def augmented(Xlog, ly, powers=(2, 3), variant="fitted"):
    """构造与 harness reset_test 完全一致的无约束设计,返回 X_u, q。"""
    n = ly.shape[0]
    X_r = np.column_stack([np.ones(n), Xlog])
    b_r, _, _, _, _ = _ols(X_r, ly)
    yhat = X_r @ b_r
    if variant == "fitted":
        added = np.column_stack([yhat ** k for k in powers])
    else:
        added = np.column_stack([Xlog[:, j] ** k
                                 for j in range(Xlog.shape[1]) for k in powers])
    return np.column_stack([X_r, added]), added.shape[1]


def statsmodels_hc3_p(X_u, ly, q):
    res = smapi.OLS(ly, X_u).fit(cov_type="HC3")
    R = np.hstack([np.zeros((q, X_u.shape[1] - q)), np.eye(q)])
    return float(res.wald_test(R, use_f=True).pvalue)


def base_seed(s, n_idx=1, rel_idx=2):              # run_srsd 的规则,n=200,rel=1%
    return 200_000 + s * 97 + n_idx * 13 + rel_idx * 5_000_000


def gen(feats_fn, s):
    df, feats = feats_fn(seed=base_seed(s), n=N, rel=REL)
    cols = feats.split(",") if isinstance(feats, str) else list(feats)
    return np.log(df[cols].values), np.log(df["target"].values)


# ================================================================ E1 复现
def task_cls(reg, tid):
    hits = [r["cls"] for r in reg if r["eq_name"].endswith(tid)]
    if not hits:
        raise KeyError(f"task {tid} not in registry")
    return hits[0]


def run_e1(reg, tasks=("i.34.27", "iii.12.43", "i.12.1"), seeds=range(60)):
    out = {}
    for tid in tasks:
        cls = task_cls(reg, tid)
        fn = make_dataset_fn(cls)
        rej = {"harness_hc3": 0, "classical_F": 0, "statsmodels": 0,
               "svd_native": 0, "n": 0, "fails": []}
        p_store = []
        for s in seeds:
            try:
                Xlog, ly = gen(fn, s)
            except ValueError as exc:
                rej["fails"].append(f"seed{s}:{exc}")
                continue
            X_u, q = augmented(Xlog, ly)
            r_h = reset_test(Xlog, ly, variant="fitted", vcov="hc3")
            r_f = reset_test(Xlog, ly, variant="fitted", vcov="classical")
            p_sm = statsmodels_hc3_p(X_u, ly, q)
            coef, resid, hat_s, XtX_s = svd_ols(X_u, ly)
            p_sv, _ = hc3_wald_from(X_u, coef, resid, hat_s, XtX_s, q)
            rej["harness_hc3"] += bool(r_h["rej"])
            rej["classical_F"] += bool(r_f["rej"])
            rej["statsmodels"] += p_sm < ALPHA
            rej["svd_native"] += p_sv < ALPHA
            rej["n"] += 1
            p_store.append({"seed": s, "p_harness": r_h["p"], "p_F": r_f["p"],
                            "p_sm": p_sm, "p_svd": p_sv,
                            "cond_X": r_h["cond"]})
        out[tid] = {**{k: v for k, v in rej.items() if k != "fails"},
                    "rate": {m: (rej[m] / rej["n"]) for m in
                             ("harness_hc3", "classical_F", "statsmodels",
                              "svd_native")},
                    "p_detail": p_store, "gen_fails": rej["fails"]}
    return out


# ================================================================ E2 载体分解
def run_e2(reg, tid="i.34.27", seeds=range(20)):
    cls = task_cls(reg, tid)
    fn = make_dataset_fn(cls)
    paths = {"gram_hat+gram_cov": [], "svd_hat+gram_cov": [],
             "gram_hat+svd_cov": [], "svd_hat+svd_cov": [],
             "statsmodels": []}
    diag = []
    for s in seeds:
        Xlog, ly = gen(fn, s)
        X_u, q = augmented(Xlog, ly)
        coef, resid, _, _ = svd_ols(X_u, ly)          # 拟合本身 SVD(与 lstsq 一致)
        _, _, _, _, XtX_G = _ols(X_u, ly)
        hat_G = np.einsum("ij,jk,ik->i", X_u, XtX_G, X_u)
        _, _, hat_S, XtX_S = svd_ols(X_u, ly)
        for name, hat, cov in (("gram_hat+gram_cov", hat_G, XtX_G),
                               ("svd_hat+gram_cov", hat_S, XtX_G),
                               ("gram_hat+svd_cov", hat_G, XtX_S),
                               ("svd_hat+svd_cov", hat_S, XtX_S)):
            p, _ = hc3_wald_from(X_u, coef, resid, hat, cov, q)
            paths[name].append(p)
        paths["statsmodels"].append(statsmodels_hc3_p(X_u, ly, q))
        G = X_u.T @ X_u
        P_G = np.linalg.pinv(G)
        diag.append({
            "seed": s,
            "log10_cond_X": float(np.log10(np.linalg.cond(X_u))),
            "log10_cond_Gram": float(np.log10(np.linalg.cond(G))),
            "gram_truncation_maxdev": float(np.max(np.abs(P_G @ G - np.eye(G.shape[0])))),
            "svd_truncation_maxdev": float(np.max(np.abs(
                np.linalg.pinv(X_u) @ X_u - np.eye(X_u.shape[1])))),
            "cov_relerr_gram_vs_svd": float(np.linalg.norm(P_G - XtX_S) /
                                            np.linalg.norm(XtX_S)),
            "hat_maxabsdiff": float(np.max(np.abs(hat_G - hat_S))),
        })
    summ = {}
    for name, ps in paths.items():
        ps = np.asarray(ps, float)
        summ[name] = {"median_p": float(np.nanmedian(ps)),
                      "rej_rate": float(np.nanmean(ps < ALPHA))}
    return {"paths": summ, "diagnostics": diag}


# ================================================================ E3 条件数扫描
H0_TASKS = ["i.12.1", "i.29.4", "i.34.27", "i.38.12", "i.39.10", "i.39.22",
            "i.43.31", "i.47.23", "ii.27.16", "ii.27.18", "ii.8.31",
            "iii.12.43", "iii.15.14"]


def run_e3(reg, seeds=range(10)):
    rows = []
    for tid in H0_TASKS:
        cls = task_cls(reg, tid)
        fn = make_dataset_fn(cls)
        inst = cls()
        widths = [float(np.log10(max(float(so.max_value), 1e-300)) -
                        np.log10(max(abs(float(so.min_value)), 1e-300)))
                  for so in inst.sampling_objs]
        n_feat = len(widths)
        rej = 0
        conds = []
        for s in seeds:
            try:
                Xlog, ly = gen(fn, s)
            except ValueError:
                continue
            X_u, q = augmented(Xlog, ly)
            r_h = reset_test(Xlog, ly, variant="fitted", vcov="hc3")
            rej += bool(r_h["rej"])
            conds.append(r_h["cond"])
        rows.append({"task": tid, "n_feat": n_feat,
                     "min_log10_domain_width": float(min(widths)),
                     "median_log10_cond_X": float(np.median(np.log10(conds))) if conds else None,
                     "rej_rate_harness_hc3": rej / max(len(conds), 1)})
    # 点二列相关:拒绝率(0/1 聚合到任务) vs 条件数
    xs = np.array([r["median_log10_cond_X"] for r in rows if r["median_log10_cond_X"]])
    ys = np.array([r["rej_rate_harness_hc3"] for r in rows if r["median_log10_cond_X"]])
    r_pb = float(np.corrcoef(xs, ys)[0, 1])
    return {"rows": rows, "point_biserial_r_cond_vs_rej": r_pb}


def main():
    reg = registry_with_meta()
    results = {"env": {"numpy": np.__version__},
               "E1_reproduction": run_e1(reg),
               "E2_carrier_decomposition": run_e2(reg),
               "E3_condition_scan": run_e3(reg)}
    path = os.path.join(HERE, "results_localization.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=1, ensure_ascii=False)
    # 控制台摘要
    print("=== E1 复现(n=200, rel=1%, 60 seeds) ===")
    for tid, blk in results["E1_reproduction"].items():
        rate = blk["rate"]
        print(f"{tid:>10}: harness_HC3 {rate['harness_hc3']:.3f} | "
              f"classical_F {rate['classical_F']:.3f} | "
              f"statsmodels {rate['statsmodels']:.3f} | "
              f"svd_native {rate['svd_native']:.3f}  (n={blk['n']})")
    print("=== E2 载体分解(i.34.27, 20 seeds) ===")
    for name, st in results["E2_carrier_decomposition"]["paths"].items():
        print(f"{name:>20}: median p {st['median_p']:.4g}  rej {st['rej_rate']:.2f}")
    d = results["E2_carrier_decomposition"]["diagnostics"]
    print(f"  median log10 cond(X) = {np.median([x['log10_cond_X'] for x in d]):.2f}, "
          f"log10 cond(X'X) = {np.median([x['log10_cond_Gram'] for x in d]):.2f}")
    print(f"  Gram 截断最大偏差 = {np.median([x['gram_truncation_maxdev'] for x in d]):.3g}, "
          f"SVD 截断最大偏差 = {np.median([x['svd_truncation_maxdev'] for x in d]):.3g}")
    print(f"  协方差 Gram-vs-SVD 相对误差 = "
          f"{np.median([x['cov_relerr_gram_vs_svd'] for x in d]):.3g}")
    print("=== E3 条件数扫描(13 任务 × 10 seeds) ===")
    e3 = results["E3_condition_scan"]
    for r in e3["rows"]:
        print(f"{r['task']:>10}: cond 10^{r['median_log10_cond_X']:.1f}  "
              f"min_log10_width {r['min_log10_domain_width']:.2f}  "
              f"feat {r['n_feat']}  rej {r['rej_rate_harness_hc3']:.2f}")
    print(f"点二列 r(条件数, 拒绝率) = {e3['point_biserial_r_cond_vs_rej']:.3f}")
    print("JSON ->", path)


if __name__ == "__main__":
    main()
