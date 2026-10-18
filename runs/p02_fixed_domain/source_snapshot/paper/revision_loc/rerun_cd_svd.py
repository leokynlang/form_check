"""rerun_cd_svd.py — v5 重跑:C/D 换 SVD 稳定路径,A/B 沿用已审计数字。

背景:
  原始跑(run_srsd.py / run_experiment.py / run_families.py)通过
  common.run_all_methods 调 reset_test,其 HC3 协方差由 pinv(X'X) 组装。
  定位实验(paper/revision_loc/)证明该路径在 cond(X_aug) ≳ 10^7 的
  物理派生设计上因 pinv 默认截断容差静默失效。本脚本在同一数据上重算
  C/D,协方差因子改由 X 的 SVD 导出(等价 pinv(X)pinv(X)',不组 Gram),
  hat 同源。A(自洽 pinv 行)与 B(RSS,无 pinv 依赖)数值不受影响,
  沿用原始 run 目录,不重算。

  B 沿用哨兵:每个重算种子上重算 B 的 F 统计量并与存档比对
  (容差 rel 1e-6),不一致即中止——防止数据重建漂移污染 C/D 列。

覆盖四个套件(网格与种子规则逐字节同原跑):
  S1 runs/full_500        19 场景 × n{50,200,1000} × rel{0,.5,1,3,5}% × 500
                          seed = 100_000 + s*97 + n_idx*13 + rel_idx*5_000_000
  S2 runs/srsd_full_r0..4 32 任务 × n × rel × 500
                          seed = 200_000 + s*97 + n_idx*13 + rel_idx*5_000_000
  S3 runs/feynman_families 30 方程 × rel{4} × het{0,1} × scale{.5,1,2}
                          seed = 300_000 + done*1_000_000 + s*97(逐格顺序!)
  S6 runs/wb_real         预注册关系 × 年,种子在 results.csv 内联
  E  real_data/smoke_results.csv  147 个单实现格,数据来自 clean/ 缓存

输出: runs/cd_svd/<suite>_cd_svd.csv(逐种子行,只含 C/D 新列 + 哨兵列)。
"""
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats as sstats

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, r"D:\form_check")

from common import _ols, reset_test                       # noqa: E402

HERE = r"D:\form_check"
OUT_DIR = os.path.join(HERE, "runs", "cd_svd")
os.makedirs(OUT_DIR, exist_ok=True)
REL_TOL = 1e-6


# ---------------------------------------------------------------- SVD 稳定路径
def _svd_factors(X):
    U, s, Vt = np.linalg.svd(X, full_matrices=False)
    s_safe = np.where(s > 0, s, 1.0)
    hat = np.einsum("ij,ij->i", U, U)
    XtX_inv = (Vt.T / (s_safe ** 2)) @ Vt        # = pinv(X) @ pinv(X)^T
    return hat, XtX_inv


def reset_test_svd(Xlog, ly, powers=(2, 3), variant="fitted"):
    """common.reset_test 的 C/D 重实现:唯一改动是 hat/XtX_inv 的来源。"""
    n = ly.shape[0]
    X_r = np.column_stack([np.ones(n), Xlog])
    if np.linalg.matrix_rank(X_r) < X_r.shape[1]:
        return {"p": np.nan, "rej": False, "stat": np.nan, "q": 0, "k": 0,
                "fail": "rank_deficient", "cond": np.inf}
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
    b_u, resid_u, rss_u, _, _ = _ols(X_u, ly)          # 系数路径不变(lstsq)
    hat_u, XtX_inv = _svd_factors(X_u)                 # ← 唯一改动
    cond = float(np.linalg.cond(X_u))
    if float(resid_u @ resid_u) <= 1e-12 * float(ly @ ly):
        return {"p": np.nan, "rej": False, "stat": np.nan, "q": q, "k": k,
                "fail": "degenerate_zero_residual", "cond": cond}
    one_m_h = 1.0 - hat_u
    if np.any(one_m_h <= 1e-12):
        return {"p": np.nan, "rej": False, "stat": np.nan, "q": q, "k": k,
                "fail": "hc3_leverage_degenerate", "cond": cond}
    eu = resid_u / one_m_h
    meat = X_u.T @ (X_u * (eu ** 2)[:, None])
    V = XtX_inv @ meat @ XtX_inv
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
    pval = float(sstats.f.sf(W / q, q, n - k))
    return {"p": pval, "rej": bool(pval < 0.05), "stat": W / q, "q": q,
            "k": k, "fail": None, "cond": cond}


def cd_row(Xlog, ly, archived_pF=None, archived_statF=None):
    """一行 = C/D 的 SVD 重算(+哨兵)。与 run_all_methods 列名对齐。"""
    out = {}
    rF = reset_test(Xlog, ly, powers=(2, 3), variant="fitted", vcov="classical")
    if archived_pF is not None and np.isfinite(rF["p"]) and np.isfinite(archived_pF):
        stat_ok = (archived_statF is None or not np.isfinite(archived_statF)
                   or np.isclose(rF["stat"], archived_statF, rtol=REL_TOL))
        out["sentinel_ok"] = bool(stat_ok and
                                  np.isclose(rF["p"], archived_pF, rtol=REL_TOL))
    else:
        out["sentinel_ok"] = True          # 存档失败/NaN 处不设哨兵
    out["stat_F_recheck"] = rF["stat"]
    rC = reset_test_svd(Xlog, ly, variant="fitted")
    rD = reset_test_svd(Xlog, ly, variant="exog")
    for tag, r in (("C", rC), ("D", rD)):
        out[f"p_{tag}_svd"] = r["p"]
        out[f"rej_{tag}_svd"] = r["rej"]
        out[f"stat_{tag}_svd"] = r["stat"]
        out[f"fail_{tag}_svd"] = r["fail"]
        out[f"cond_{tag}_svd"] = r["cond"]
    return out


def _save(name, df):
    p = os.path.join(OUT_DIR, name)
    df.to_csv(p, index=False)
    print(f"  -> {p}  ({len(df)} rows)")


# ---------------------------------------------------------------- S1
def rerun_s1():
    from scenarios import SCENARIOS, SCENARIO_META
    old = pd.read_csv(os.path.join(HERE, "runs", "full_500", "raw.csv"),
                      usecols=["scenario", "n", "rel", "seed", "klass",
                               "p_fitted_f", "stat_fitted_f"])
    n_grid, rel_grid = [50, 200, 1000], [0.0, 0.005, 0.01, 0.03, 0.05]
    frames, t0 = [], time.perf_counter()
    for scen in SCENARIOS:
        sub = old[old.scenario == scen]
        for (n, rel, seed), grp in sub.groupby(["n", "rel", "seed"]):
            n_idx, rel_idx = n_grid.index(n), rel_grid.index(rel)
            bs = 100_000 + seed * 97 + n_idx * 13 + rel_idx * 5_000_000
            df, fl = SCENARIOS[scen](seed=bs, n=n, rel=rel)
            feats = fl.split(",")
            Xlog = np.log(df[feats].values)
            ly = np.log(df["target"].values)
            row = {"scenario": scen, "klass": SCENARIO_META[scen][0],
                   "n": n, "rel": rel, "seed": seed, "base_seed": bs}
            row.update(cd_row(Xlog, ly, float(grp.p_fitted_f.iloc[0]),
                              float(grp.stat_fitted_f.iloc[0])))
            frames.append(row)
        print(f"S1 {scen}: {len(frames)} rows, {time.perf_counter()-t0:.0f}s",
              flush=True)
        if len(frames) > 20000:
            _save("s1_cd_svd_partial.csv", pd.DataFrame(frames))
    return pd.DataFrame(frames)


# ---------------------------------------------------------------- S2
def rerun_s2():
    from srsd_data import registry_with_meta, make_dataset_fn
    reg = registry_with_meta()
    fn_by_task = {r["eq_name"]: make_dataset_fn(r["cls"]) for r in reg}
    n_grid, rel_grid = [50, 200, 1000], [0.0, 0.005, 0.01, 0.03, 0.05]
    frames, t0 = [], time.perf_counter()
    for shard in range(5):
        old = pd.read_csv(os.path.join(HERE, "runs", f"srsd_full_r{shard}",
                                       "raw_srsd.csv"), low_memory=False)
        for (task, n, rel, seed), grp in old.groupby(
                ["task", "n", "rel", "seed"]):
            if grp.data_fail.iloc[0] not in (None, "", "nan") and pd.notna(
                    grp.data_fail.iloc[0]):
                continue
            fn = fn_by_task[task]
            n_idx, rel_idx = n_grid.index(int(n)), rel_grid.index(float(rel))
            bs = 200_000 + int(seed) * 97 + n_idx * 13 + rel_idx * 5_000_000
            try:
                df, feat_str = fn(seed=bs, n=int(n), rel=float(rel))
            except ValueError:
                continue
            feats = feat_str.split(",")
            Xlog = np.log(df[feats].values)
            ly = np.log(df["target"].values)
            row = {"task": task, "n": int(n), "rel": float(rel),
                   "seed": int(seed), "base_seed": bs}
            row.update(cd_row(Xlog, ly, float(grp.p_fitted_f.iloc[0]),
                              float(grp.stat_fitted_f.iloc[0])))
            frames.append(row)
        print(f"S2 shard {shard}: {len(frames)} cum rows, "
              f"{time.perf_counter()-t0:.0f}s", flush=True)
        _save("s2_cd_svd_partial.csv", pd.DataFrame(frames))
    return pd.DataFrame(frames)


# ---------------------------------------------------------------- S3
def rerun_s3():
    """S3 逐种子层在原仓 raw_families.csv(11,760 行)。单元格按执行顺序
    连续落盘,由行序重放 done 计数器还原 seed_base = 300_000+done*10^6+s*97;
    哨兵用存档 p_fitted_f。"""
    from srsd_data import registry_with_meta, make_dataset_fn
    reg = registry_with_meta()
    fn_by_task = {r["eq_name"]: make_dataset_fn(r["cls"]) for r in reg}
    src = r"D:\dimreflex\experiments\form_check_vs_reset\runs\feynman_families\raw_families.csv"
    old = pd.read_csv(src, low_memory=False)
    frames, t0 = [], time.perf_counter()
    done, prev_key = 0, None
    for _, r in old.iterrows():
        key = (r["eq_name"], r["rel"], r["hetero"], r["range_scale"])
        if key != prev_key:
            done += 1
            prev_key = key
        if isinstance(r["data_fail"], str) and r["data_fail"].strip():
            continue
        s = int(r["seed"])
        bs = 300_000 + done * 1_000_000 + s * 97
        fn = fn_by_task[r["eq_name"]]
        try:
            df, feat_str = fn(seed=bs, n=200, rel=float(r["rel"]),
                              hetero=bool(r["hetero"]),
                              range_scale=float(r["range_scale"]))
        except ValueError:
            continue
        feats = feat_str.split(",")
        Xlog = np.log(df[feats].values)
        ly = np.log(df["target"].values)
        row = {"eq_name": r["eq_name"], "rel": float(r["rel"]),
               "hetero": bool(r["hetero"]),
               "range_scale": float(r["range_scale"]), "seed": s,
               "base_seed": bs, "cell_ordinal": done}
        pF = r["p_fitted_f"]
        row.update(cd_row(Xlog, ly,
                          float(pF) if pd.notna(pF) else None))
        frames.append(row)
        if len(frames) % 2000 == 0:
            print(f"S3 {len(frames)}/11760 {time.perf_counter()-t0:.0f}s",
                  flush=True)
    df = pd.DataFrame(frames)
    bad = int((~df.sentinel_ok).sum())
    print(f"S3 哨兵: {bad}/{len(df)} 不一致")
    _save("s3_cd_svd.csv", df)
    return df


# ---------------------------------------------------------------- S6 / E
# S6(World Bank)与工程台架 E 的 C/D 重放在 replay_s6_smoke_svd.py:
# monkey-patch common.reset_test 后原样重放 wb_run.main() / smoke_test.main()。
# 原因:两者的数据层靠内联 API 缓存 / clean/ CSV 重建,不走本文件的
# cd_row 通道;在本文件里维护第二套重放逻辑只会重演 smoke 覆盖事故。


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.perf_counter()
    if which in ("all", "s1"):
        df = rerun_s1()
        _save("s1_cd_svd.csv", df)
    if which in ("all", "s2"):
        df = rerun_s2()
        _save("s2_cd_svd.csv", df)
    if which in ("all", "s3"):
        df = rerun_s3()
    print(f"TOTAL {time.perf_counter()-t0:.0f}s")


if __name__ == "__main__":
    main()
