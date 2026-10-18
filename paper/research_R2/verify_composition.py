"""verify_composition.py — 为 manuscript.md 的修订逐条验证要写进论文的数字。

口径与 reproduce_r2.py 完全一致（sentinel_ok==True 且 task != feynman-i.30.5）。
所有输出写入 verify_composition.json，供逐字核对。
"""
import json
import sys
from pathlib import Path

import numpy as np

OUT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT_DIR))
from reproduce_r2 import ERR_T, build  # noqa: E402

OUT = OUT_DIR / "evidence" / "verify_composition.json"


def main():
    _, pool = build(True, True)
    res = {}
    D = {}

    # 1. 按任务的 Da / Ld（pooled，R2 口径）
    per_task = []
    for task, g in pool.groupby("task"):
        Da = g[g.extrap_far_relerr > ERR_T]
        Ld = Da[Da._cert & (Da.r2 >= 0.99)]
        per_task.append({
            "task": task.replace("feynman-", ""),
            "records": int(len(g)),
            "S": int((g.r2 >= 0.99).sum()),
            "Sd": int(((g.r2 >= 0.99) & (g.extrap_far_relerr > ERR_T)).sum()),
            "Ld": int(len(Ld)),
            "Da": int(len(Da)),
        })
    per_task.sort(key=lambda r: -r["Ld"])
    res["per_task_pooled"] = per_task
    D["A. 按任务 pooled（R2 口径）"] = per_task

    # 2. bonus 两任务的分 n 贡献
    bonus = {}
    for task in ("feynman-bonus.15", "feynman-bonus.5"):
        g = pool[pool.task == task]
        row = {}
        for n in (50, 200, 1000):
            gn = g[g.n == n]
            row[n] = {
                "records": int(len(gn)),
                "Sd": int(((gn.r2 >= 0.99) & (gn.extrap_far_relerr > ERR_T)).sum()),
                "Ld": int((gn._cert & (gn.r2 >= 0.99)
                           & (gn.extrap_far_relerr > ERR_T)).sum()),
            }
        bonus[task.replace("feynman-", "")] = row
    res["bonus_by_n"] = bonus
    D["B. bonus 两任务分 n"] = bonus

    # 3. 排除 bonus.15 后的 Ld/Sd
    excl = {}
    for n in (50, 200, 1000):
        g = pool[(pool.n == n) & (pool.task != "feynman-bonus.15")]
        S = g[g.r2 >= 0.99]
        Sd = S[S.extrap_far_relerr > ERR_T]
        Ld = Sd[Sd._cert]
        excl[n] = {"S": int(len(S)), "Sd": int(len(Sd)), "Ld": int(len(Ld)),
                   "Ld/Sd": round(len(Ld) / len(Sd), 5) if len(Sd) else None}
    res["excluding_bonus15"] = excl
    D["C. 排除 bonus.15 后"] = excl

    # 4. bonus.15 far 分布（R²≥0.99 层）
    dist = {}
    for task in ("feynman-bonus.15", "feynman-bonus.5", "feynman-ii.38.14"):
        s = pool[(pool.task == task) & (pool.r2 >= 0.99)].extrap_far_relerr
        s = s[np.isfinite(s)]
        if len(s):
            dist[task.replace("feynman-", "")] = {
                "n": int(len(s)),
                "q25": round(float(np.percentile(s, 25)), 4),
                "median": round(float(np.median(s)), 4),
                "q75": round(float(np.percentile(s, 75)), 4),
                "q90": round(float(np.percentile(s, 90)), 4),
                "q99": round(float(np.percentile(s, 99)), 4),
                "max": round(float(s.max()), 4),
                "gt10": round(float((s > 0.10).mean()), 5),
                "gt15": round(float((s > 0.15).mean()), 5),
            }
        else:
            dist[task.replace("feynman-", "")] = {"n": 0}
    res["far_distribution"] = dist
    D["D. far 分布（R²≥0.99 层）"] = dist

    # 5. 负面对照：Da 全为 0 放行的任务
    zero_rel = [r["task"] for r in per_task if r["Da"] > 0 and r["Ld"] == 0]
    res["zero_release_tasks"] = zero_rel
    D["E. 危险池放行为 0 的任务"] = zero_rel

    # 6. 阈值扫描下的 bonus.15 占比（复核报告 §3.1）
    sweep = []
    for t in (0.02, 0.05, 0.10, 0.15, 0.20, 0.50):
        S = pool[pool.r2 >= 0.99]
        Sd = S[S.extrap_far_relerr > t]
        Ld = Sd[Sd._cert]
        b15 = int((Ld.task == "feynman-bonus.15").sum())
        b5 = int((Ld.task == "feynman-bonus.5").sum())
        sweep.append({"threshold": t, "Sd": int(len(Sd)), "Ld": int(len(Ld)),
                      "bonus15": b15, "bonus5": b5, "other": int(len(Ld)) - b15 - b5,
                      "Ld/Sd": round(len(Ld) / len(Sd), 5) if len(Sd) else None})
    res["threshold_sweep"] = sweep
    D["F. 阈值扫描"] = sweep

    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(res, f, indent=2, ensure_ascii=False)

    for k, v in D.items():
        print("=" * 72)
        print(k)
        if isinstance(v, list) and v and isinstance(v[0], dict):
            for r in v:
                print("  " + "  ".join(f"{kk}={vv}" for kk, vv in r.items()))
        else:
            print("  " + json.dumps(v, ensure_ascii=False))
    print()
    print("written:", OUT)


if __name__ == "__main__":
    main()