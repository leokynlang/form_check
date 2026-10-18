"""Task-level uncertainty for the joint-release counts (read-only).

The record-level percentages in Section 4.10 treat simulation replicates as
independent observations. This audit instead treats the 17 eligible H2
equations as sampling units: it recomputes S/Sd/Ld/Da per task, bootstraps
tasks with replacement, and reports leave-one-task-out ranges. No fitting or
permutation testing is performed; inputs are the same caches used by
audit_results.py.

Run: python audit_task_intervals.py
"""
from __future__ import annotations

import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

import audit_results as ar

OUT = Path(__file__).resolve().parent
EVIDENCE = OUT / "evidence"
REPS = 10000
SEED = 20261018


def task_counts(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for task, g in df.groupby("task"):
        rows.append({
            "task": task,
            "N": int(len(g)),
            "S": int(g.screen.sum()),
            "Sd": int((g.screen & g.danger).sum()),
            "Ld": int((g.certificate & g.danger).sum()),
            "Da": int(g.danger.sum()),
        })
    return pd.DataFrame(rows).sort_values("task").reset_index(drop=True)


def pooled_ratios(counts: np.ndarray) -> dict:
    S, Sd, Ld, Da = (float(x) for x in counts)
    return {
        "Ld_over_S": Ld / S if S else None,
        "Ld_over_Sd": Ld / Sd if Sd else None,
        "Ld_over_Da": Ld / Da if Da else None,
    }


def bootstrap_ci(counts: np.ndarray, reps: int = REPS, seed: int = SEED) -> dict:
    rng = np.random.default_rng(seed)
    n = len(counts)
    draw = rng.integers(0, n, size=(reps, n))
    sums = counts[draw].sum(axis=1)
    out = {}
    for j, key in enumerate(["Ld_over_S", "Ld_over_Sd", "Ld_over_Da"]):
        num = sums[:, 2]
        den = sums[:, [0, 1, 3][j]]
        with np.errstate(divide="ignore", invalid="ignore"):
            vals = np.where(den > 0, num / den, np.nan)
        vals = vals[np.isfinite(vals)]
        out[key] = {
            "mean": float(vals.mean()),
            "ci95": np.quantile(vals, [.025, .975]).tolist(),
            "share_exact_zero": float((vals == 0).mean()),
        }
    return out


def leave_one_out(counts: np.ndarray, tasks: list[str]) -> dict:
    out = {}
    for j, key in enumerate(["Ld_over_S", "Ld_over_Sd", "Ld_over_Da"]):
        vals = []
        for i, task in enumerate(tasks):
            sub = np.delete(counts, i, axis=0).sum(axis=0)
            vals.append({"omitted": task, "value": pooled_ratios(sub)[key]})
        finite = [v["value"] for v in vals if v["value"] is not None]
        out[key] = {"min": min(finite), "max": max(finite), "values": vals}
    return out


def main() -> None:
    m = ar.load()
    valid = m[(m.sentinel_ok == True) & (m.task != "feynman-i.30.5")].copy()
    all_ex = valid[(valid.rel > 0) & np.isfinite(valid.r2)
                   & np.isfinite(valid.extrap_far_relerr)].copy()
    h1 = all_ex[all_ex.klass == "H1"].copy()
    assert h1.task.nunique() == 17

    result = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "audit_script": {"path": str(Path(__file__).resolve()),
                         "sha256": ar.sha256(__file__)},
        "inputs": {k: {"path": str(v), "sha256": ar.sha256(v), "bytes": v.stat().st_size}
                   for k, v in ar.INPUTS.items()},
        "runtime": {"python": platform.python_version(),
                    "numpy": np.__version__, "pandas": pd.__version__},
        "method": (
            "Task-level bootstrap: the 17 eligible H1 equations are resampled "
            "with replacement (B=10000); S/Sd/Ld/Da are summed over the drawn "
            "tasks and the pooled ratios recomputed. Leave-one-task-out ranges "
            "omit one equation at a time. This measures between-equation "
            "heterogeneity, not record-level binomial precision; the task "
            "bootstrap does not include threshold-estimation uncertainty."),
        "definitions": {
            "S": "screened wrong models (finite r2 >= .99)",
            "Sd": "dangerous screened wrong models (extrap_far_relerr > .10)",
            "Ld": "dangerous records passing the screen and all four valid gates",
            "Da": "all dangerous eligible H1 records",
        },
        "by_n": {},
        "pooled": {},
    }

    for n in [50, 200, 1000]:
        c = task_counts(h1[h1.n == n])
        counts = c[["S", "Sd", "Ld", "Da"]].to_numpy(dtype=float)
        result["by_n"][str(n)] = {
            "tasks": c.task.tolist(),
            "counts": c.to_dict(orient="records"),
            "pooled": pooled_ratios(counts.sum(axis=0)),
            "task_bootstrap_ci95": bootstrap_ci(counts),
            "leave_one_task_out": leave_one_out(counts, c.task.tolist()),
            "per_task": [
                {"task": r.task, "Sd": int(r.Sd), "Ld": int(r.Ld),
                 "Ld_over_Sd": (r.Ld / r.Sd) if r.Sd else None}
                for r in c.itertuples(index=False)
            ],
        }

    c = task_counts(h1)
    counts = c[["S", "Sd", "Ld", "Da"]].to_numpy(dtype=float)
    result["pooled"] = {
        "tasks": c.task.tolist(),
        "counts": c.to_dict(orient="records"),
        "pooled": pooled_ratios(counts.sum(axis=0)),
        "task_bootstrap_ci95": bootstrap_ci(counts),
        "leave_one_task_out": leave_one_out(counts, c.task.tolist()),
        "tasks_with_zero_Ld": int((c.Ld == 0).sum()),
        "tasks_with_positive_Ld": c[c.Ld > 0].task.tolist(),
        "share_of_Ld_top_two_tasks": float(
            c.Ld.nlargest(2).sum() / c.Ld.sum()) if c.Ld.sum() else None,
    }

    EVIDENCE.mkdir(exist_ok=True)
    destination = EVIDENCE / "task_intervals.json"
    destination.write_text(json.dumps(result, indent=2, ensure_ascii=False,
                                      allow_nan=False), encoding="utf-8")
    print(json.dumps({
        "output": str(destination),
        "pooled": result["pooled"]["pooled"],
        "pooled_ci95": result["pooled"]["task_bootstrap_ci95"],
        "pooled_loo": {k: [v["min"], v["max"]] for k, v in result["pooled"]["leave_one_task_out"].items()},
        "by_n": {n: {"pooled": v["pooled"], "ci95": v["task_bootstrap_ci95"]}
                 for n, v in result["by_n"].items()},
        "zero_Ld_tasks": result["pooled"]["tasks_with_zero_Ld"],
        "positive_Ld_tasks": result["pooled"]["tasks_with_positive_Ld"],
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
