"""Recount cached joint-silence evidence; never generate or fit data.

Inputs are read-only files under D:/form_check. The only output is
evidence/s7.json beside this script. Run with Python, pandas and a parquet
reader installed. No benchmark module is imported (some run on import).
"""

from __future__ import annotations

import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


SOURCE = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent / "evidence" / "s7.json"
FILES = {
    "s2_original": SOURCE / "runs/srsd_full/raw_srsd.parquet",
    "s2_svd": SOURCE / "runs/cd_svd/s2_cd_svd.csv",
    "s5_results": SOURCE / "runs/gate_cmp/results.csv",
    "s5_thresholds": SOURCE / "runs/gate_cmp/summary.json",
    "four_gate_script": SOURCE / "joint_silence.py",
    "three_gate_script": SOURCE / "joint_silence_toolbox.py",
    "s5_runner_source": SOURCE / "gate_compare.py",
    "current_wrapper_source": SOURCE / "common.py",
    "archived_four_gate_script": SOURCE / "paper/revision_loc/joint_silence_v5.py",
    "stored_four_gate_summary": SOURCE / "runs/joint_silence/joint_silence_srsd_pool.csv",
    "stored_three_gate_summary": SOURCE / "runs/joint_silence/toolbox_silence.csv",
    "stored_three_gate_h0": SOURCE / "runs/joint_silence/toolbox_silence_h0_reference.csv",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def record(mask: pd.Series) -> dict:
    count = int(mask.sum())
    return {"silent": count, "total": len(mask), "rate": count / len(mask)}


def likelihood(h0: dict, h1: dict) -> dict:
    p0, p1 = h0["rate"], h1["rate"]
    return {
        "definition": "P(silent|H0)/P(silent|H1), using the same arm and threshold regime",
        "LR_H0_over_H1": p0 / p1 if p1 else None,
        "posterior_H1_given_silent_at_equal_prior": p1 / (p0 + p1) if p0 + p1 else None,
        "prior_H1": 0.5,
        "interpretation": "Descriptive pool likelihood ratio; posterior requires the stated mixing prior.",
    }


def audit() -> dict:
    # Match original pandas parsing to preserve threshold-boundary decisions.
    base = pd.read_parquet(FILES["s2_original"], columns=[
        "task", "n", "rel", "seed", "p_fc", "p_fitted_f", "klass", "rej_fc", "rej_fitted_f",
    ])
    svd = pd.read_csv(FILES["s2_svd"], usecols=[
        "task", "n", "rel", "seed", "sentinel_ok", "p_C_svd", "p_D_svd", "rej_C_svd", "rej_D_svd",
    ])
    svd = svd[svd.sentinel_ok.eq(True) & svd.task.ne("feynman-i.30.5")]
    joined = base.merge(svd, on=["task", "n", "rel", "seed"],
                        how="inner", validate="one_to_one")
    joined = joined[joined.rel.gt(0)]
    pcols = ["p_fc", "p_fitted_f", "p_C_svd", "p_D_svd"]
    rcols = ["rej_fc", "rej_fitted_f", "rej_C_svd", "rej_D_svd"]
    four = []
    for n, group in joined.groupby("n"):
        h0 = group[group.klass.eq("H0")]
        thresholds = h0[pcols].quantile(0.05)
        entry = {
            "n": int(n), "calibration_thresholds": thresholds.to_dict(),
            "threshold_rule": "Per-gate H0-pool 5th percentile, pooled over positive noise layers at this n; reject p < threshold.",
            "counts": {},
        }
        for klass, values in group.groupby("klass"):
            clean = np.isfinite(values[pcols]).all(axis=1)
            stored_silent = clean & ~values[rcols].astype(bool).any(axis=1)
            strict_silent = clean & values[pcols].ge(0.05).all(axis=1)
            assert stored_silent.equals(strict_silent)
            entry["counts"][klass] = {
                "task_count": int(values.task.nunique()),
                "tasks": sorted(values.task.unique().tolist()),
                "rows_per_task": values.groupby("task").size().astype(int).to_dict(),
                "seed_ids_per_task_noise_cell": sorted(values.groupby(["task", "rel"]).seed.nunique().unique().tolist()),
                "rel_values": sorted(values.rel.unique().tolist()),
                "declined_any": int((~clean).sum()),
                "native": record(strict_silent),
                "native_stored_decisions": record(stored_silent),
                "inclusive_threshold_sensitivity_not_native": record(clean & values[pcols].gt(0.05).all(axis=1)),
                "calibrated": record(clean & values[pcols].ge(thresholds).all(axis=1)),
            }
        entry["likelihood"] = {
            regime: likelihood(entry["counts"]["H0"][regime], entry["counts"]["H1"][regime])
            for regime in ("native", "calibrated")
        }
        four.append(entry)

    gates = pd.read_csv(FILES["s5_results"])
    with FILES["s5_thresholds"].open(encoding="utf-8") as stream:
        summary = json.load(stream)
    three = []
    for n, group in gates[gates.bench.eq("srsd")].groupby("n"):
        thresholds = summary["cells"][f"srsd_n{n}"]["naive"]
        entry = {
            "n": int(n), "E_threshold": thresholds["d_star"], "F_threshold": thresholds["m_star"],
            "threshold_rule": "A uses stored rej_fc; E rejects D_cv>d_star; F rejects delta_aicc>m_star. E/F thresholds fitted separately on this S5 H0 pool, not a joint 5% budget.",
            "counts": {},
        }
        for klass, values in group.groupby("klass"):
            a_quiet = ~values.rej_fc.astype(bool) & values.fail_fc.isna()
            failure = values[["fail_fc", "fail_mc"]].notna().any(axis=1)
            ef_reject = values.D_cv.gt(thresholds["d_star"]) | values.delta_aicc.gt(thresholds["m_star"])
            aef_quiet = a_quiet & ~failure & ~ef_reject
            assert not (aef_quiet & ~a_quiet).any()
            entry["counts"][klass] = {
                "task_count": int(values.task.nunique()),
                "tasks": sorted(values.task.unique().tolist()),
                "A": record(a_quiet), "AEF": record(aef_quiet),
                "paired_A_silent_removed_by_EF": int((a_quiet & ~aef_quiet).sum()),
                "declined_any": int(failure.sum()),
                "nonfinite_scores_without_failure": int((~np.isfinite(values[["D_cv", "delta_aicc"]]).all(axis=1) & ~failure).sum()),
            }
        entry["likelihood"] = likelihood(entry["counts"]["H0"]["AEF"], entry["counts"]["H1"]["AEF"])
        three.append(entry)

    overlap = gates[gates.bench.eq("srsd")].merge(base, on=["task", "n", "seed"], how="inner")
    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "execution": "Read-only aggregation of cached results; no new fitting, sampling, or benchmark imports.",
        "runtime": {"python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__},
        "inputs": {key: {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)} for key, path in FILES.items()},
        "audit_script_sha256": sha256(Path(__file__).resolve()),
        "four_gate_arm": {
            "methods": ["A: permutation curvature", "B: fitted RESET F", "C: fitted HC3 SVD", "D: exog HC3 SVD"],
            "filters": ["SVD sentinel_ok == True", "task != feynman-i.30.5", "inner join original S2 and SVD by task/n/rel/seed", "rel > 0", "NaN p-value prevents silence"],
            "provenance": "Reanalysis of existing S2 fits, not independent S7 sampling. 4 positive noise levels times 500 seed IDs give 2000 records per task/n.",
            "native_threshold_audit": {
                "rule": "Reject p < .05, not p <= .05; silent requires p >= .05 and all p finite. This matches joint_silence.py:43-50 and stored A/B/C/D decisions on the eligible joined pool.",
                "per_gate": {
                    pcol: {
                        "stored_decision_column": rcol,
                        "p_equals_005": int(joined[pcol].eq(.05).sum()),
                        "stored_vs_strict_mismatches": int((joined[rcol].astype(bool) != joined[pcol].lt(.05)).sum()),
                        "stored_vs_inclusive_mismatches": int((joined[rcol].astype(bool) != joined[pcol].le(.05)).sum()),
                        "infinite_p_values": int(np.isinf(joined[pcol]).sum()),
                    }
                    for pcol, rcol in zip(pcols, rcols)
                },
                "scope": "Only the eligible S2/SVD joined pool is checked; no claim is made that all suites or historical source versions use the same convention.",
            },
            "results": four,
        },
        "three_gate_arm": {
            "methods": ["A: permutation curvature", "E: CV-quad", "F: AICc-quad"],
            "filters": ["gate_cmp bench == srsd", "klass H0 or H1", "recorded fail_fc or fail_mc prevents silence"],
            "design": {"rel": 0.03, "seeds_per_task": 200, "seed_rule": "400000 + s*97 + n_index*13", "SVD_C_D_present": False},
            "sentinel_evidence": {
                "runner_calls": "gate_compare.py:160 calls common.form_check_with_p; lines 161-162 retain p_fc/rej_fc/fail_fc only.",
                "current_wrapper_checks": "common.py:127 rejects unchecked native results; lines 145-147 compare native and wrapper curvature and record statistic_mismatch.",
                "stored_columns": list(gates.columns),
                "explicit_replay_sentinel_ok_column": "sentinel_ok" in gates.columns,
                "native_checked_column": "checked" in gates.columns,
                "all_recorded_fail_fc": gates.fail_fc.fillna("NONE").value_counts().astype(int).to_dict(),
                "all_recorded_fail_mc": gates.fail_mc.fillna("NONE").value_counts().astype(int).to_dict(),
                "i305": {
                    "rows": int(gates.task.eq("feynman-i.30.5").sum()),
                    "rows_by_n": {str(n): int(v) for n, v in gates[gates.task.eq("feynman-i.30.5")].groupby("n").size().items()},
                    "recorded_fail_fc": gates.loc[gates.task.eq("feynman-i.30.5"), "fail_fc"].fillna("NONE").value_counts().astype(int).to_dict(),
                    "recorded_fail_mc": gates.loc[gates.task.eq("feynman-i.30.5"), "fail_mc"].fillna("NONE").value_counts().astype(int).to_dict(),
                    "missing_p_fc": int(gates.loc[gates.task.eq("feynman-i.30.5"), "p_fc"].isna().sum()),
                },
                "qualification": "Current source contains a native-wrapper drift check and saved S5 rows report no failure, but CSV lacks replay sentinel_ok and raw native/recomputed statistics. No run-level code-hash manifest was located in runs/gate_cmp. This is not equivalent to the separately replayed S2/SVD sentinel evidence and does not prove historical implementation validity.",
                "scope_decision": "Retain as exploratory 19-H1 S5 A/E/F analysis. No recorded S5 i.30.5 mismatch justifies exclusion solely on the separate S2 replay outcome; keep the same S5 task policy in Table 12 and this arm.",
            },
            "results": three,
        },
        "cross_arm": {
            "same_task_n_seed_key_matches": len(overlap),
            "same_pool_four_vs_six_from_current_inputs": False,
            "reason": "S5 results have only A/E/F and no B/C/D; S2 records have A-D but no E/F on those draws. Cannot infer a paired 38.6% to 42.2% increase.",
        },
        "limitations": [
            "Reference H0 results use the same pool used to calibrate thresholds; not independent calibration validation.",
            "No shared joint false-positive budget is imposed on A versus A/E/F.",
            "Equations, not repeated datasets, are the units for between-equation generalization.",
            "Legacy revision_loc/joint_silence_v5.py omits rel>0 in the SRSD split, unlike the root script.",
            "This audit does not validate generator correctness, original model fitting, or the statistical tests.",
        ],
    }


if __name__ == "__main__":
    result = audit()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(OUTPUT)
