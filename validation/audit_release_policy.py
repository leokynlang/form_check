"""Independent record/accounting and historical-replay checks for P0-1.

Does not import the new policy calibration implementation. It checks frozen
threshold decisions and denominator arithmetic directly from saved records.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
import sys

sys.path.insert(0, str(ROOT))

from common import form_check_with_p, reset_test
from experiments.release_validation.inventory import digest
from experiments.release_validation.runner import fit_and_endpoint, load_registry
from paper.revision_loc.rerun_cd_svd import reset_test_svd


def replay(out):
    raw = pd.read_parquet(ROOT / "runs/srsd_full/raw_srsd.parquet")
    quality = pd.read_parquet(ROOT / "runs/srsd_full/quality_recompute.parquet")
    svd = pd.read_csv(ROOT / "runs/cd_svd/s2_cd_svd.csv", low_memory=False)
    keys = ["task", "n", "rel", "seed"]
    rows = raw[(raw.seed == 0) & (raw.rel == .03) & raw.n.isin([50, 200])].merge(
        quality[[*keys, "r2"]], on=keys, validate="one_to_one").merge(
            svd[[*keys, "sentinel_ok", "p_C_svd", "p_D_svd", "fail_C_svd", "fail_D_svd"]],
            on=keys, validate="one_to_one")
    registry = {r["task"]: r for r in load_registry()}
    checks = []
    for row in rows.itertuples():
        rec = registry[row.task]
        data, features = rec["fn"](int(row.base_seed), int(row.n), float(row.rel))
        features = features.split(",")
        x, y = np.log(data[features].to_numpy()), np.log(data.target.to_numpy())
        outcome = fit_and_endpoint(data, features, rec["truth"])
        a = form_check_with_p(data, "target", features, int(row.base_seed))
        b = reset_test(x, y, variant="fitted", vcov="classical")
        current = {"r2": outcome["r2"], "endpoint": outcome["legacy_endpoint_error"],
                   "A": a["p_fc"], "B": b["p"],
                   "C": reset_test_svd(x, y, variant="fitted")["p"],
                   "D": reset_test_svd(x, y, variant="exog")["p"]}
        archived = {"r2": row.r2, "endpoint": row.extrap_far_relerr,
                    "A": row.p_fc, "B": row.p_fitted_f, "C": row.p_C_svd, "D": row.p_D_svd}
        parity = {key: bool(np.isclose(current[key], archived[key], rtol=1e-6,
                                       atol=1e-10, equal_nan=True)) for key in current}
        checks.append({"task": row.task, "n": row.n, "historical_sentinel_valid": bool(row.sentinel_ok),
                       "parity": parity,
                       "absolute_differences": {key: (float(abs(current[key] - archived[key]))
                                                       if np.isfinite(current[key]) and np.isfinite(archived[key])
                                                       else None) for key in current}})
    failures = [r for r in checks if r["historical_sentinel_valid"] and not all(r["parity"].values())]
    result = {"rows": len(checks), "checked_seed": 0, "noise": .03,
              "sentinel_valid_replay_mismatches": failures,
              "sentinel_invalid_records_are_diagnostic_only": True, "checks": checks}
    out.mkdir(parents=True, exist_ok=True)
    (out / "replay.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({"replayed": len(checks), "sentinel_valid_mismatches": len(failures),
                      "output": str(out / 'replay.json')}))
    if failures:
        raise SystemExit("Replay parity failed; inspect per-metric differences before trusting reuse")


def audit(run):
    manifest = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["complete"]
    raw_path = run / "raw_policy.parquet"
    assert digest(raw_path) == manifest["raw_sha256"]
    raw = pd.read_parquet(raw_path)
    decision = pd.read_parquet(run / "decisions.parquet")
    summary = json.loads((run / "summary.json").read_text(encoding="utf-8"))
    policies = json.loads((run / "policies.json").read_text(encoding="utf-8"))
    protocol = manifest["protocol"]
    assert len(raw) == manifest["rows"] == 6592
    assert raw.base_seed.nunique() == len(raw)
    assert raw.dataset_sha256.dropna().nunique() == raw.dataset_sha256.notna().sum()
    assert raw.loc[raw.phase != "test", "klass"].eq("H0").all()
    assert not decision.duplicated(["policy", "budget", "task", "n", "rel", "replicate"]).any()
    assert (decision.released == (decision.screen & ~decision.alarm & ~decision.decline)).all()
    checked = 0
    for fitted in policies:
        if not fitted["fitted"]["feasible"]:
            continue
        n, rel = fitted["n"], fitted["rel"]
        members, refs = fitted["members"], fitted["reference"]
        cal = raw[(raw.n == n) & (raw.rel == rel) & (raw.phase == "calibration")
                  & raw.r2.ge(protocol["screen_r2"])].copy()
        scores, valid = [], np.ones(len(cal), dtype=bool)
        for member in members:
            ref_rows = raw[(raw.n == n) & (raw.rel == rel) & (raw.phase == "reference")
                           & raw.r2.ge(protocol["screen_r2"]) & raw[f"fail_{member}"].isna()
                           & np.isfinite(raw[f"score_{member}"])]
            expected_ref = np.sort(ref_rows[f"score_{member}"].to_numpy())
            assert np.array_equal(expected_ref, np.array(refs[member]))
            v = cal[f"score_{member}"].to_numpy()
            valid &= np.isfinite(v) & cal[f"fail_{member}"].isna().to_numpy()
            if not len(expected_ref):
                valid[:] = False
            scores.append(np.array([np.sum(expected_ref < value) / (len(expected_ref) + 1) for value in v]))
        joint = np.max(scores, axis=0)
        withheld = ~valid | (joint > fitted["fitted"]["threshold"])
        assert withheld.sum() <= int(np.floor(fitted["budget"] * len(cal) + 1e-10))
        d = decision[(decision.n == n) & (decision.rel == rel) & (decision.policy == fitted["policy"])
                     & (decision.budget == fitted["budget"])]
        assert len(d) == int(((raw.n == n) & (raw.rel == rel) & (raw.phase == "test")).sum())
        stratum = next(s for s in summary["strata"] if s["n"] == n and s["rel"] == rel)
        counts = stratum["budgets"][str(fitted["budget"])]["policies"][fitted["policy"]]["counts"]
        danger = np.isfinite(d.legacy_endpoint_error) & d.legacy_endpoint_error.gt(protocol["danger_threshold"])
        assert int((d.klass.eq("H1") & danger & d.released).sum()) == counts["h1_danger_released"]
        assert int((d.klass.eq("H1") & danger & d.screen).sum()) == counts["h1_danger_screened"]
        assert int((d.klass.eq("H0") & d.screen & (d.alarm | d.decline)).sum()) == counts["h0_withheld_screened"]
        assert int((d.released & ~np.isfinite(d.legacy_endpoint_error)).sum()) == counts["released_unknown_outcome"]
        checked += 1
    result = {"complete": True, "rows": len(raw), "decision_rows": len(decision),
              "checked_policy_strata": checked, "raw_sha256": digest(raw_path),
              "auditor_sha256": digest(Path(__file__)),
              "checks": ["source input hash", "phase/seed/dataset separation", "reference H0 provenance",
                         "strict joint calibration including declines", "release identity",
                         "risk denominators and missing outcome accounting"]}
    (run / "audit.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", default="runs/p01_pilot_v1")
    parser.add_argument("--replay", action="store_true")
    args = parser.parse_args()
    if args.replay:
        replay(ROOT / "runs/p01_inventory")
    else:
        audit(ROOT / args.run)
