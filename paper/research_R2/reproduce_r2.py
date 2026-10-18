"""reproduce_r2.py — 用 research_R2/audit_results.py 的确切口径复现 Table 15。

口径（来自 audit_results.py）：
  valid = sentinel_ok == True 且 task != feynman-i.30.5
  外推池 = valid 中 rel>0 且 r2 有限 且 extrap_far_relerr 有限
  screen = r2 >= .99
  danger = extrap_far_relerr > .10
  certificate = screen 且四门 p 有限 且无 typed decline 且无 native-alpha 拒绝

对照：我 2026-10-18 的第一版复现（reproduce_c1.py）未排除 i.30.5，
      也未用 sentinel_ok 过滤。
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = Path(__file__).resolve().parent
assert OUT_DIR == ROOT / "paper" / "research_R2", "Unexpected script directory"
ERR_T = 0.10


def build(use_sentinel=True, drop_i305=True):
    raw = pd.read_parquet(ROOT / "runs/srsd_full/raw_srsd.parquet")
    q = pd.read_parquet(ROOT / "runs/srsd_full/quality_recompute.parquet")
    cd = pd.read_csv(ROOT / "runs/cd_svd/s2_cd_svd.csv", low_memory=False)
    df = raw.merge(q[["task", "n", "rel", "seed", "r2"]],
                   on=["task", "n", "rel", "seed"], how="left")
    df = df.merge(cd[["task", "n", "rel", "seed", "sentinel_ok",
                      "p_C_svd", "p_D_svd", "rej_C_svd", "rej_D_svd",
                      "fail_C_svd", "fail_D_svd"]],
                  on=["task", "n", "rel", "seed"], how="left")

    if use_sentinel:
        df = df[df.sentinel_ok == True]  # noqa: E712
    if drop_i305:
        df = df[df.task != "feynman-i.30.5"]

    rej = df[["rej_fc", "rej_fitted_f", "rej_C_svd", "rej_D_svd"]]
    fail = df[["fail_fc", "fail_fitted_f", "fail_C_svd", "fail_D_svd"]].notna()
    p_ok = np.isfinite(df[["p_fc", "p_fitted_f", "p_C_svd", "p_D_svd"]]).all(axis=1)
    # certificate：四门 p 有限、无 typed decline、无 native-alpha 拒绝
    df["_cert"] = p_ok & (~fail.any(axis=1)) & (~rej.fillna(True).any(axis=1))

    pool = df[(df.klass == "H1") & (df.rel > 0)
              & np.isfinite(df.r2) & np.isfinite(df.extrap_far_relerr)]
    return df, pool


def table15(pool):
    rows = []
    for n in (50, 200, 1000):
        g = pool[pool.n == n]
        S = g[g.r2 >= 0.99]
        Sd = S[S.extrap_far_relerr > ERR_T]
        Ld = Sd[Sd._cert]
        Da = g[g.extrap_far_relerr > ERR_T]
        rows.append({"n": n, "records": len(g), "S": len(S), "Sd": len(Sd),
                     "Ld": len(Ld), "Da": len(Da),
                     "Ld/S": round(len(Ld) / len(S), 5) if len(S) else np.nan,
                     "Ld/Sd": round(len(Ld) / len(Sd), 5) if len(Sd) else np.nan,
                     "Ld/Da": round(len(Ld) / len(Da), 5) if len(Da) else np.nan})
    t = pd.DataFrame(rows)
    pooled = {"S": int(t.S.sum()), "Sd": int(t.Sd.sum()),
              "Ld": int(t.Ld.sum()), "Da": int(t.Da.sum())}
    pooled["Ld/S"] = round(pooled["Ld"] / pooled["S"], 5)
    pooled["Ld/Sd"] = round(pooled["Ld"] / pooled["Sd"], 5)
    pooled["Ld/Da"] = round(pooled["Ld"] / pooled["Da"], 5)
    return t, pooled


def gate_only_n50(pool):
    """R2 §4.10 文字：n=50 有 966 个危险记录四门全绿（不限 R²）。"""
    g = pool[(pool.n == 50) & (pool.extrap_far_relerr > ERR_T)]
    return int(g._cert.sum())


if __name__ == "__main__":
    print("=" * 72)
    print("A. R2 口径（sentinel_ok + 排除 i.30.5）")
    _, pool = build(True, True)
    t, pooled = table15(pool)
    print(t.to_string(index=False))
    print("pooled:", pooled)
    print(f"n=50 门卫全绿的危险记录（不限 R²）: {gate_only_n50(pool)}")
    print(f"任务数: {pool.task.nunique()}  每 n 记录: "
          f"{pool[pool.n == 50].shape[0]}")

    print()
    print("=" * 72)
    print("B. 我的第一版口径（无 sentinel 过滤、不排 i.30.5）")
    _, pool2 = build(False, False)
    t2, pooled2 = table15(pool2)
    print(t2.to_string(index=False))
    print("pooled:", pooled2)
    print(f"n=50 门卫全绿的危险记录: {gate_only_n50(pool2)}")
    print(f"任务数: {pool2.task.nunique()}  每 n 记录: "
          f"{pool2[pool2.n == 50].shape[0]}")

    print()
    print("=" * 72)
    print("C. 差异来源：i.30.5 与 i.26.2 的记录数")
    _, p_all = build(False, False)
    for task in ("feynman-i.30.5", "feynman-i.26.2"):
        s = p_all[p_all.task == task]
        print(f"  {task}: 记录 {len(s)}  n=50 {len(s[s.n == 50])}  "
              f"有限far {int(np.isfinite(s.extrap_far_relerr).sum())}")