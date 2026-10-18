"""threshold_decomposition.py — R2 口径下的阈值扫描 + 任务级分解。

R2 §4.10 已报总体放行率随阈值的变化（41.30%/15.85%/6.40% at 2/10/20%），
但未做任务级分解。本脚本补那一步：481 个放行记录在不同阈值下由哪些
任务贡献，检验 bonus.15 的贡献是否为阈值切分效应。

口径与 research_R2/audit_results.py 完全一致。
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = Path(__file__).resolve().parent
assert OUT_DIR == ROOT / "paper" / "research_R2", "Unexpected script directory"
OUT = OUT_DIR / "evidence" / "threshold_decomposition.json"


def build():
    raw = pd.read_parquet(ROOT / "runs/srsd_full/raw_srsd.parquet")
    q = pd.read_parquet(ROOT / "runs/srsd_full/quality_recompute.parquet")
    cd = pd.read_csv(ROOT / "runs/cd_svd/s2_cd_svd.csv", low_memory=False)
    df = raw.merge(q[["task", "n", "rel", "seed", "r2"]],
                   on=["task", "n", "rel", "seed"], how="left")
    df = df.merge(cd[["task", "n", "rel", "seed", "sentinel_ok",
                      "rej_C_svd", "rej_D_svd", "fail_C_svd", "fail_D_svd",
                      "p_C_svd", "p_D_svd"]],
                  on=["task", "n", "rel", "seed"], how="left")
    df = df[(df.sentinel_ok == True) & (df.task != "feynman-i.30.5")]  # noqa: E712
    rej = df[["rej_fc", "rej_fitted_f", "rej_C_svd", "rej_D_svd"]]
    fail = df[["fail_fc", "fail_fitted_f", "fail_C_svd", "fail_D_svd"]].notna()
    p_ok = np.isfinite(df[["p_fc", "p_fitted_f", "p_C_svd", "p_D_svd"]]).all(axis=1)
    df["_cert"] = p_ok & (~fail.any(axis=1)) & (~rej.fillna(True).any(axis=1))
    pool = df[(df.klass == "H1") & (df.rel > 0)
              & np.isfinite(df.r2) & np.isfinite(df.extrap_far_relerr)]
    return df, pool


def sweep(pool, n_fixed=None):
    out = []
    for th in (0.02, 0.05, 0.10, 0.15, 0.20, 0.50):
        p = pool if n_fixed is None else pool[pool.n == n_fixed]
        S = p[p.r2 >= 0.99]
        Sd = S[S.extrap_far_relerr > th]
        Ld = Sd[Sd._cert]
        by_task = Ld.groupby("task").size().sort_values(ascending=False)
        b15 = int(by_task.get("feynman-bonus.15", 0))
        b5 = int(by_task.get("feynman-bonus.5", 0))
        out.append({
            "threshold": th,
            "S": int(len(S)), "Sd": int(len(Sd)),
            "Sd_frac": round(len(Sd) / len(S), 5),
            "Ld": int(len(Ld)),
            "Ld_over_S": round(len(Ld) / len(S), 5) if len(S) else None,
            "Ld_over_Sd": round(len(Ld) / len(Sd), 5) if len(Sd) else None,
            "bonus15": b15, "bonus5": b5,
            "bonus15_share": round(b15 / len(Ld), 4) if len(Ld) else None,
            "bonus5_share": round(b5 / len(Ld), 4) if len(Ld) else None,
            "other": int(len(Ld) - b15 - b5),
            "by_task": {k: int(v) for k, v in by_task.items()},
        })
    return out


def far_quantiles(pool, task):
    s = pool[(pool.task == task) & (pool.r2 >= 0.99)]
    v = s.extrap_far_relerr.dropna()
    if len(v) == 0:
        return None
    return {"n": int(len(v)),
            "q01": round(float(v.quantile(.01)), 4),
            "q25": round(float(v.quantile(.25)), 4),
            "median": round(float(v.median()), 4),
            "q75": round(float(v.quantile(.75)), 4),
            "q90": round(float(v.quantile(.90)), 4),
            "q99": round(float(v.quantile(.99)), 4),
            "frac_gt_10pct": round(float((v > .10).mean()), 4),
            "frac_gt_15pct": round(float((v > .15).mean()), 4)}


if __name__ == "__main__":
    df, pool = build()
    res = {
        "pool": {"tasks": int(pool.task.nunique()),
                 "records_per_n": int((pool.n == 50).sum()),
                 "total": int(len(pool))},
        "pooled_sweep": sweep(pool),
        "n50_sweep": sweep(pool, 50),
        "n200_sweep": sweep(pool, 200),
        "n1000_sweep": sweep(pool, 1000),
        "far_quantiles": {t: far_quantiles(pool, t)
                          for t in ("feynman-bonus.15", "feynman-bonus.5",
                                    "feynman-ii.38.14", "feynman-i.27.6")},
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2, default=str)

    print(f'池: {res["pool"]["tasks"]} 任务 x {res["pool"]["records_per_n"]}/n')
    print()
    print("=== pooled 阈值扫描（R2 口径）===")
    print(f'{"th":>5} {"S":>7} {"Sd":>6} {"Sd_frac":>8} {"Ld":>5} '
          f'{"Ld/Sd":>8} {"b15":>5} {"b15sh":>7} {"b5":>5} {"other":>5}')
    for e in res["pooled_sweep"]:
        print(f'{e["threshold"]:>5.2f} {e["S"]:>7} {e["Sd"]:>6} {e["Sd_frac"]:>8.4f} '
              f'{e["Ld"]:>5} {str(e["Ld_over_Sd"]):>8} {e["bonus15"]:>5} '
              f'{str(e["bonus15_share"]):>7} {e["bonus5"]:>5} {e["other"]:>5}')
    print()
    print("=== n=50 阈值扫描（与 R2 §4.10 同档）===")
    for e in res["n50_sweep"]:
        print(f'{e["threshold"]:>5.2f} Sd_frac={e["Sd_frac"]:.4f} '
              f'Ld/Sd={e["Ld_over_Sd"]} Ld={e["Ld"]} '
              f'bonus15={e["bonus15"]} share={e["bonus15_share"]}')
    print()
    print("=== far 分布分位 ===")
    for t, q in res["far_quantiles"].items():
        if q:
            print(f'  {t:20s} n={q["n"]:>5} med={q["median"]:.4f} '
                  f'q25={q["q25"]:.4f} q75={q["q75"]:.4f} q90={q["q90"]:.4f} '
                  f'>10%={q["frac_gt_10pct"]:.4f} >15%={q["frac_gt_15pct"]:.4f}')
    print(f"\n已写入 {OUT}")