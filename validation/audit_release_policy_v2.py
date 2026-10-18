"""Independent numerical and provenance audit for the frozen P0-1 v2 run.

This file deliberately imports no calibration, runner, or analysis implementation.
It reconstructs the reference transforms, thresholds, decisions, and all event
counts directly from the saved raw records. It writes only into runs/p01_v2.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
METHODS = tuple("ABCDEFG")
POLICY_KEYS = ["n", "rel", "calibration_round", "policy", "budget"]
TEST_KEYS = ["task", "n", "rel", "replicate"]
COUNT_NAMES = [
    "total", "released", "h0_total", "h0_screened", "h0_released",
    "h0_withheld_screened", "h0_alarm_screened", "h0_decline_screened",
    "h1_total", "h1_alarm", "h1_decline", "h1_danger_total",
    "h1_danger_screened", "h1_danger_released", "h1_safe_screened",
    "h1_safe_released", "all_danger_released", "h0_danger_released",
    "released_known_outcome", "outcome_unavailable", "released_unknown_outcome",
    "data_failures", "fit_failures",
]
RATIOS = {
    "H0_withholding_given_screen": ("h0_withheld_screened", "h0_screened"),
    "H0_alarm_given_screen": ("h0_alarm_screened", "h0_screened"),
    "H0_decline_given_screen": ("h0_decline_screened", "h0_screened"),
    "H0_retention_overall": ("h0_released", "h0_total"),
    "H1_alarm_yield": ("h1_alarm", "h1_total"),
    "H1_decline_rate": ("h1_decline", "h1_total"),
    "dangerous_H1_release_given_screen": ("h1_danger_released", "h1_danger_screened"),
    "dangerous_H1_release_overall": ("h1_danger_released", "h1_danger_total"),
    "safe_H1_retention_given_screen": ("h1_safe_released", "h1_safe_screened"),
    "release_overall": ("released", "total"),
    "dangerous_H1_share_of_known_outcome_releases": ("h1_danger_released", "released_known_outcome"),
    "all_dangerous_share_of_known_outcome_releases": ("all_danger_released", "released_known_outcome"),
    "outcome_unavailable_rate": ("outcome_unavailable", "total"),
}
PAIRS = [("AD_panel", "D"), ("AD_panel", "F"), ("AEF_panel", "F"),
         ("G_rbf_cv", "F"), ("D", "F")]


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def require(condition, description):
    if not bool(condition):
        raise AssertionError(description)


def close(a, b, description):
    require(np.allclose(a, b, rtol=1e-12, atol=0, equal_nan=True), description)


def raw_transform(frame, method):
    raw = frame[f"score_{method}"].to_numpy(dtype=float)
    failed = frame[f"fail_{method}"].fillna("").astype(str).str.strip().ne("").to_numpy()
    valid = np.isfinite(raw) & ~failed
    if method in "ABCD":
        p = -raw
        valid &= (p >= 0) & (p <= 1)
        with np.errstate(all="ignore"):
            value = -np.log10(np.maximum(p, np.finfo(np.float64).tiny))
        value[~valid] = np.nan
    else:
        value = raw.copy()
        value[~valid] = np.nan
    return value, valid


def independent_reference(frame, method):
    """Return independently computed center/scale diagnostics, without metadata."""
    transformed, valid = raw_transform(frame, method)
    values = transformed[valid]
    p = -frame[f"score_{method}"].to_numpy(dtype=float)
    p_counts = {
        "zero": int(np.sum(p == 0)) if method in "ABCD" else 0,
        "floor": int(np.sum(p == .001)) if method == "A" else 0,
        "floating_floor": int(np.sum((p >= 0) & (p <= np.finfo(np.float64).tiny)))
        if method in "ABCD" else 0,
    }
    if not len(values):
        return {"usable": False, "center": None, "scale": None, "n": 0,
                "candidates": {}, **p_counts}
    center = float(np.median(values))
    q25, q75 = np.quantile(values, [.25, .75])
    candidates = {
        "iqr_normalized": float((q75 - q25) / 1.3489795003921634),
        "mad_normalized": float(np.median(np.abs(values - center)) * 1.482602218505602),
        "population_sd": float(np.std(values, ddof=0)),
    }
    return {"usable": None, "center": center, "scale": None, "n": len(values),
            "candidates": candidates, **p_counts}


def audit_reference(frame, method, saved):
    expected = independent_reference(frame, method)
    require(int(saved["n"]) == expected["n"], f"{method}: reference valid count")
    for field, expected_name in [("reference_p_zero_count", "zero"),
                                 ("reference_p_floor_count", "floor"),
                                 ("floating_floor_count", "floating_floor")]:
        require(int(saved[field]) == expected[expected_name], f"{method}: {field}")
    require(saved["raw_transform"] == ("negative_log10_p" if method in "ABCD" else "identity"),
            f"{method}: raw transform metadata")
    if not expected["n"]:
        require(not saved["usable"], f"{method}: empty reference must be unusable")
        require(saved["reason"] == "empty_reference" and saved["scale_method"] == "unusable",
                f"{method}: empty reference metadata")
        return
    close(saved["center"], expected["center"], f"{method}: reference median")
    # Scale names are human-readable metadata; the numerical value must match
    # the fixed fallback ordering independently of that spelling.
    candidates = expected["candidates"]
    available = [(name, value) for name, value in candidates.items()
                 if np.isfinite(value) and value > 0]
    if not available:
        require(not saved["usable"], f"{method}: constant reference must be unusable")
        require(saved["reason"] == "degenerate_reference_scale" and saved["scale_method"] == "unusable",
                f"{method}: constant reference metadata")
        return
    require(saved["usable"], f"{method}: nonconstant reference must be usable")
    close(saved["scale"], available[0][1], f"{method}: reference scale fallback")
    require(saved["scale_method"] == available[0][0] and saved["reason"] is None,
            f"{method}: reference fallback metadata")
    require(np.isfinite(saved["scale"]) and saved["scale"] > 0,
            f"{method}: positive finite affine scale")


def independent_policy_score(frame, members, reference):
    valid = np.ones(len(frame), dtype=bool)
    transformed = []
    for method in members:
        score, member_valid = raw_transform(frame, method)
        ref = reference[method]
        if not ref["usable"]:
            member_valid[:] = False
            transformed.append(np.full(len(frame), np.nan))
        else:
            with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                affine = (score - ref["center"]) / ref["scale"]
            transformed.append(affine)
            member_valid &= np.isfinite(affine)
        valid &= member_valid
    joint = np.max(np.column_stack(transformed), axis=1)
    joint[~valid] = np.nan
    return joint, ~valid


def independent_threshold(score, decline, budget):
    n = len(score)
    cap = int(np.floor(float(budget) * n + 1e-10))
    failures = int(np.sum(decline))
    expected = {"n": n, "declines": failures, "allowed_withheld": cap,
                "feasible": bool(n and failures <= cap)}
    if not expected["feasible"]:
        expected["threshold"] = None
        return expected
    values = np.sort(score[~decline])
    # Budgets are below one, hence at least one non-declined observation must
    # remain. Strict exceedance makes the threshold conservative under ties.
    threshold = float(values[len(values) - (cap - failures) - 1])
    expected.update(threshold=threshold, alarms=int(np.sum(values > threshold)))
    expected["withheld_rate"] = (expected["alarms"] + failures) / n
    return expected


def independently_count(frame, danger_threshold):
    h0, h1 = frame.klass.eq("H0"), frame.klass.eq("H1")
    screen, alarm, decline, released = [frame[x].astype(bool)
                                      for x in ("screen", "alarm", "decline", "released")]
    known = np.isfinite(frame.legacy_endpoint_error)
    danger = known & frame.legacy_endpoint_error.gt(danger_threshold)
    safe = known & ~danger
    events = [
        np.ones(len(frame), dtype=bool), released, h0, h0 & screen, h0 & released,
        h0 & screen & (alarm | decline), h0 & screen & alarm, h0 & screen & decline,
        h1, h1 & alarm, h1 & decline, h1 & danger, h1 & danger & screen,
        h1 & danger & released, h1 & safe & screen, h1 & safe & released,
        danger & released, h0 & danger & released, released & known, ~known,
        released & ~known, frame.data_failure.notna(), frame.fit_failure.notna(),
    ]
    table = pd.DataFrame({name: np.asarray(event, dtype=np.int64)
                          for name, event in zip(COUNT_NAMES, events, strict=True)}, index=frame.index)
    return table.groupby(frame.task, sort=True).sum()


def compare_nullable(saved, expected, label):
    if expected is None:
        require(saved is None, label)
    else:
        close(saved, expected, label)


def independent_mean_summary(original, weighted, round_draws, all_feasible, saved):
    finite = np.isfinite(original)
    values = original[finite]
    defined = bool(all_feasible and finite.all())
    require(saved["primary_defined"] == defined, "primary estimate feasibility/definition")
    compare_nullable(saved["point_mean"], float(np.mean(original)) if defined else None, "round-equal mean")
    compare_nullable(saved["feasible_defined_round_only_sensitivity_mean"],
                     float(np.mean(values)) if len(values) else None, "secondary defined-round mean")
    sensitivity = saved["calibration_sensitivity"]
    require(sensitivity["defined_rounds"] == int(finite.sum()), "defined calibration rounds")
    for field, expected in [("sample_sd", float(np.std(values, ddof=1)) if len(values) > 1 else None),
                            ("min", float(np.min(values)) if len(values) else None),
                            ("max", float(np.max(values)) if len(values) else None),
                            ("empirical_percentile_range95", np.quantile(values, [.025, .975]).tolist()
                             if len(values) else None)]:
        compare_nullable(sensitivity[field], expected, f"calibration sensitivity {field}")
    if not defined:
        require(all(saved[name] is None for name in ("conditional_task", "round_mean", "crossed")),
                "undefined primary estimate has no main interval")
        return
    draws = {
        "conditional_task": np.mean(weighted, axis=1),
        "round_mean": np.mean(original[round_draws], axis=1),
        "crossed": np.mean(weighted[np.arange(len(weighted))[:, None], round_draws], axis=1),
    }
    for field, values in draws.items():
        finite = np.isfinite(values)
        expected = np.quantile(values[finite], [.025, .975]).tolist() if finite.any() else None
        require(saved[field]["defined_draws"] == int(finite.sum())
                and saved[field]["total_draws"] == len(values), f"{field}: defined draws")
        compare_nullable(saved[field]["ci95"], expected, f"{field}: independent interval")


def independent_ratio_arrays(cube, numerator, denominator, weights):
    num = cube[:, :, COUNT_NAMES.index(numerator)]
    den = cube[:, :, COUNT_NAMES.index(denominator)]
    original_num, original_den = num.sum(axis=1), den.sum(axis=1)
    original = np.full(original_num.shape, np.nan)
    np.divide(original_num, original_den, out=original, where=original_den > 0)
    weighted_num, weighted_den = weights @ num.T, weights @ den.T
    weighted = np.full(weighted_num.shape, np.nan)
    np.divide(weighted_num, weighted_den, out=weighted, where=weighted_den > 0)
    return original, weighted


def audit_uncertainty(run, manifest, raw, counts, policies):
    """Recreate all round summaries and three paired resampling distributions."""
    summary = read_json(run / "summary.json")
    protocol = manifest["protocol"]
    require(summary["protocol_id"] == protocol["protocol_id"]
            and summary["raw_sha256"] == manifest["raw_sha256"], "summary source identity")
    require(summary["analysis_sha256"] == digest(ROOT / "experiments/release_validation/analyze_v2.py"),
            "summary analysis source digest")
    planned_strata = {(n, rel) for n in protocol["n_grid"] for rel in protocol["noise_grid"]}
    require(len(summary["strata"]) == len(planned_strata)
            and {(s["n"], s["rel"]) for s in summary["strata"]} == planned_strata,
            "complete summary stratum inventory")
    n_rounds, repetitions = protocol["calibration_rounds"], protocol["uncertainty"]["bootstrap_repetitions"]
    lookup = {tuple(p[x] for x in POLICY_KEYS): p["fitted"] for p in policies}
    count_indices = counts.groupby(POLICY_KEYS, sort=False).indices
    round_metrics = pd.read_csv(run / "round_metrics.csv", float_precision="round_trip")
    require(not round_metrics.duplicated(POLICY_KEYS).any(), "unique round-metric keys")
    metric_indices = round_metrics.set_index(POLICY_KEYS)
    require(set(metric_indices.index) == set(lookup), "complete round-metric inventory")
    checked_metrics, checked_paired = 0, 0
    for stratum in summary["strata"]:
        n, rel = stratum["n"], stratum["rel"]
        test = raw[(raw.phase == "test") & raw.n.eq(n) & raw.rel.eq(rel)]
        task_classes = test.groupby("task", sort=True).klass.first()
        tasks = task_classes.index
        require(stratum["test_rows"] == len(test)
                and stratum["test_H0_tasks"] == int(task_classes.eq("H0").sum())
                and stratum["test_H1_tasks"] == int(task_classes.eq("H1").sum()), "summary test inventory")
        rng = np.random.default_rng(protocol["master_seed"] + int(n) + 9001)
        weights = np.zeros((repetitions, len(tasks)), dtype=np.int64)
        for klass in ("H0", "H1"):
            positions = np.flatnonzero(task_classes.to_numpy() == klass)
            selected = rng.choice(positions, size=(repetitions, len(positions)), replace=True)
            for position in positions:
                weights[:, position] = (selected == position).sum(axis=1)
        round_draws = rng.integers(0, n_rounds, size=(repetitions, n_rounds))
        require(set(stratum["budgets"]) == {str(b) for b in protocol["budgets"]},
                "complete summary budget inventory")
        for budget in protocol["budgets"]:
            result = stratum["budgets"][str(budget)]
            require(set(result["policies"]) == set(protocol["policies"]), "complete summary policy inventory")
            require(set(result["paired"]) == {f"{a}_minus_{b}" for a, b in PAIRS},
                    "complete paired comparison inventory")
            cube_by_policy, feasible_by_policy = {}, {}
            for policy in protocol["policies"]:
                saved = result["policies"][policy]
                cube = np.full((n_rounds, len(tasks), len(COUNT_NAMES)), np.nan)
                feasible = np.zeros(n_rounds, dtype=bool)
                for round_no in range(n_rounds):
                    key = (n, rel, round_no, policy, budget)
                    fitted = lookup[key]
                    row = metric_indices.loc[key]
                    feasible[round_no] = fitted["feasible"]
                    require(bool(row.feasible) == bool(fitted["feasible"]), "round-metric feasibility")
                    if not fitted["feasible"]:
                        continue
                    table = counts.iloc[count_indices[key]].set_index("task").reindex(tasks)
                    cube[round_no] = table[COUNT_NAMES].to_numpy()
                    total = table[COUNT_NAMES].sum()
                    require(np.array_equal(row[COUNT_NAMES].to_numpy(dtype=float), total.to_numpy(dtype=float)),
                            "round-metric event counts")
                    close(row.threshold, fitted["threshold"], "round-metric threshold")
                    close(row.calibration_withheld, fitted["withheld_rate"], "round-metric calibration rate")
                    for metric, (num, den) in RATIOS.items():
                        expected = float(total[num] / total[den]) if total[den] else np.nan
                        close(row[metric], expected, f"round-metric {metric}")
                cube_by_policy[policy], feasible_by_policy[policy] = cube, feasible
                require(saved["feasible_rounds"] == int(feasible.sum())
                        and saved["total_rounds"] == n_rounds
                        and saved["all_rounds_feasible"] == bool(feasible.all()), "summary round feasibility")
                for metric, (num, den) in RATIOS.items():
                    original, weighted = independent_ratio_arrays(cube, num, den, weights)
                    stats = saved["metrics"][metric]
                    independent_mean_summary(original, weighted, round_draws, bool(feasible.all()), stats)
                    contributing = np.sum(cube[:, :, COUNT_NAMES.index(den)] > 0, axis=1).astype(int).tolist()
                    require(stats["denominator_contributing_test_tasks_by_round"] == contributing,
                            "contributing denominator tasks by round")
                    if metric == "H0_withholding_given_screen":
                        require(stats["test_rounds_above_nominal"] == int(np.sum(np.isfinite(original) & (original > budget))),
                                "round count above nominal H0 budget")
                        compare_nullable(stats["test_fraction_rounds_above_nominal"],
                                         float(np.mean(original > budget))
                                         if feasible.all() and np.isfinite(original).all() else None,
                                         "fraction rounds above nominal H0 budget")
                    checked_metrics += 1
            for first, second in PAIRS:
                record = result["paired"][f"{first}_minus_{second}"]
                for metric in ("dangerous_H1_release_given_screen", "H0_retention_overall", "H0_withholding_given_screen"):
                    a, aw = independent_ratio_arrays(cube_by_policy[first], *RATIOS[metric], weights)
                    b, bw = independent_ratio_arrays(cube_by_policy[second], *RATIOS[metric], weights)
                    feasible = feasible_by_policy[first] & feasible_by_policy[second]
                    denominator_column = COUNT_NAMES.index(RATIOS[metric][1])
                    require(np.array_equal(cube_by_policy[first][feasible, :, denominator_column],
                                           cube_by_policy[second][feasible, :, denominator_column]),
                            "paired policy denominators on feasible matched rounds")
                    independent_mean_summary(a - b, aw - bw, round_draws, bool(feasible.all()), record[metric])
                    checked_paired += 1
    return {"policy_metric_summaries": checked_metrics, "paired_metric_summaries": checked_paired,
            "bootstrap_repetitions": repetitions, "independently_recreated_intervals":
            ["conditional_task", "round_mean", "crossed"]}


def audit_provenance(run, manifest, raw):
    raw_path = run / "raw_policy.parquet"
    require(digest(raw_path) == manifest["raw_sha256"], "raw data digest")
    frozen_protocol = read_json(run / "protocol_frozen.json")
    require(frozen_protocol["protocol"] == manifest["protocol"], "frozen protocol identity")
    protocol_hash = hashlib.sha256(json.dumps(manifest["protocol"], sort_keys=True).encode()).hexdigest()
    require(protocol_hash == frozen_protocol["protocol_sha256"], "frozen protocol digest")
    frozen = {name: manifest[name] for name in ("protocol", "source_sha256", "environment", "registry")}
    fingerprint = hashlib.sha256(json.dumps(frozen, sort_keys=True).encode()).hexdigest()
    require(fingerprint == manifest["fingerprint"] == frozen_protocol["fingerprint"],
            "protocol/source/environment fingerprint")
    cell_keys = ["phase", "calibration_round", "task", "n", "rel"]
    raw_cell_indices = raw.groupby(cell_keys, sort=False).indices
    cell_key_inventory = set()
    for relative, value in manifest["cells"].items():
        path = run / "cells" / relative
        require(digest(path) == value, f"cell digest: {relative}")
        cell = pd.read_parquet(path)
        identity = tuple(cell.iloc[0][name] for name in cell_keys)
        require(identity in raw_cell_indices and identity not in cell_key_inventory,
                f"cell/raw identity: {relative}")
        cell_key_inventory.add(identity)
        raw_cell = raw.iloc[raw_cell_indices[identity]]
        common = sorted(raw.columns)
        require(set(cell.columns) == set(common), f"cell/raw columns: {relative}")
        pd.testing.assert_frame_equal(
            cell[common].sort_values("replicate").reset_index(drop=True),
            raw_cell[common].sort_values("replicate").reset_index(drop=True),
            check_dtype=False, check_exact=True, obj=f"cell/raw parity {relative}")
    require(cell_key_inventory == set(raw_cell_indices), "raw/cell complete inventory parity")
    for name, value in manifest.get("source_sha256", {}).items():
        path = Path(name)
        require(digest(path if path.is_absolute() else ROOT / path) == value,
                f"frozen live source digest: {name}")
    snapshots = read_json(run / "source_snapshot.json")["files"]
    require(set(snapshots) == set(manifest["source_sha256"]), "complete source snapshot inventory")
    for name, metadata in snapshots.items():
        path = Path(metadata["snapshot"])
        require(metadata["sha256"] == manifest["source_sha256"][name],
                f"snapshot matches frozen source digest: {name}")
        require(digest(path if path.is_absolute() else run / path) == metadata["sha256"],
                f"frozen source snapshot digest: {name}")
    v1 = ROOT / "runs/p01_pilot_v1"
    snapshot = manifest["v1_recursive_snapshot"]
    current = {str(p.relative_to(v1)).replace("\\", "/"): digest(p)
               for p in v1.rglob("*") if p.is_file()}
    normalized = {str(k).replace("\\", "/"): v for k, v in snapshot.items()}
    require(current == normalized, "all v1 files preserved, no addition/removal/change")
    legacy_seeds, legacy_hashes = set(), set()
    checked_inputs = []
    for name, metadata in manifest["historical_inputs"].items():
        path = Path(name)
        if not path.is_absolute():
            path = ROOT / path
        value = metadata["sha256"] if isinstance(metadata, dict) else metadata
        require(digest(path) == value, f"historical input digest: {name}")
        if path.suffix == ".parquet":
            old = pd.read_parquet(path)
        elif path.suffix == ".csv":
            old = pd.read_csv(path, low_memory=False)
        else:
            checked_inputs.append(str(path))
            continue
        for col in ("base_seed", "seed"):
            if col in old:
                legacy_seeds.update(int(x) for x in old[col].dropna())
        if "dataset_sha256" in old:
            legacy_hashes.update(old.dataset_sha256.dropna())
        checked_inputs.append(str(path))
    old_v1 = pd.read_parquet(v1 / "raw_policy.parquet")
    legacy_seeds.update(int(x) for x in old_v1.base_seed)
    legacy_hashes.update(old_v1.dataset_sha256.dropna())
    require(not (set(raw.base_seed) & legacy_seeds), "new seeds disjoint from v1/legacy")
    require(not (set(raw.dataset_sha256.dropna()) & legacy_hashes),
            "new dataset hashes disjoint from all archived available hashes")
    return {"historical_inputs_checked": checked_inputs,
            "v1_preserved_files": len(snapshot), "legacy_seed_overlap": 0,
            "available_historical_dataset_hash_overlap": 0,
            "legacy_content_limit": "Legacy tables without dataset hashes are checked by frozen input digests and disjoint generator seeds; their dataset bytes are not exhaustively regenerated."}


def audit(run):
    run = Path(run).resolve()
    require(run == (ROOT / "runs/p01_v2").resolve(), "audit output restricted to v2 directory")
    manifest = read_json(run / "manifest.json")
    protocol = manifest["protocol"]
    require(manifest["complete"], "completed data run")
    raw = pd.read_parquet(run / "raw_policy.parquet")
    rounds = int(protocol["calibration_rounds"])
    registry = manifest["registry"]
    n_h0 = sum(r["klass"] == "H0" for r in registry)
    expected_rows = len(protocol["n_grid"]) * len(protocol["noise_grid"]) * (
        rounds * n_h0 * (protocol["replicates"]["reference"] + protocol["replicates"]["calibration"])
        + len(registry) * protocol["replicates"]["test"])
    require(len(raw) == int(manifest["rows"]) == expected_rows, "planned raw row count")
    require(not raw.duplicated(["phase", "calibration_round", *TEST_KEYS]).any(), "unique raw keys")
    require(not raw.base_seed.duplicated().any(), "all raw generator seeds unique")
    require(not raw.dataset_sha256.dropna().duplicated().any(), "all successful raw dataset hashes unique")
    require(set(raw.phase) == {"reference", "calibration", "test"}, "exact phase inventory")
    classes = {record["task"]: record["klass"] for record in registry}
    require(raw.task.map(classes).equals(raw.klass), "registered task and H0/H1 identity")
    require(set(raw.n) == set(protocol["n_grid"]) and set(raw.rel) == set(protocol["noise_grid"]),
            "exact raw n/noise inventory")
    require(raw.loc[raw.phase.ne("test"), "klass"].eq("H0").all(), "no H1 in reference/calibration")
    require(raw.loc[raw.phase.eq("test"), "calibration_round"].eq(-1).all(), "single independent test bank")
    for phase in ("reference", "calibration"):
        require(set(raw.loc[raw.phase.eq(phase), "calibration_round"]) == set(range(rounds)),
                f"{phase}: exact round inventory")
    for phase, group in raw.groupby("phase"):
        expected = set(range(protocol["replicates"][phase]))
        expected_cells = {
            (round_no, record["task"], n, rel)
            for round_no in ([-1] if phase == "test" else range(rounds))
            for record in registry if phase == "test" or record["klass"] == "H0"
            for n in protocol["n_grid"] for rel in protocol["noise_grid"]
        }
        actual_cells = set(group.groupby(["calibration_round", "task", "n", "rel"]).groups)
        require(actual_cells == expected_cells, f"{phase}: exact planned cells")
        for _, cell in group.groupby(["calibration_round", "task", "n", "rel"]):
            require(set(cell.replicate) == expected and len(cell) == len(expected), "cell replicate inventory")
    phase_numbers = {"reference": 11, "calibration": 12, "test": 13}
    task_words = {task: (int.from_bytes(hashlib.sha256(task.encode()).digest()[:4], "little"),
                         int.from_bytes(hashlib.sha256(task.encode()).digest()[4:8], "little"))
                  for task in classes}
    for row in raw.itertuples(index=False):
        entropy = [protocol["master_seed"], 2, phase_numbers[row.phase], row.calibration_round + 1,
                   *task_words[row.task], row.n, round(row.rel * 1_000_000), row.replicate]
        expected_seed = int(np.random.SeedSequence(entropy).generate_state(1, dtype=np.uint64)[0]) % (2**62) + 2**40
        require(row.base_seed == expected_seed, "independent reconstruction of every generator seed")
    provenance = audit_provenance(run, manifest, raw)
    raw["screen"] = np.isfinite(raw.r2) & raw.r2.ge(protocol["screen_r2"])
    policies = read_json(run / "policies.json")
    decision = pd.read_parquet(run / "decisions.parquet")
    counts = pd.read_csv(run / "per_task_counts.csv")
    require(not decision.duplicated([*POLICY_KEYS, "task", "replicate"]).any(), "unique decision keys")
    require(not counts.duplicated([*POLICY_KEYS, "task"]).any(), "unique task-count keys")
    d_indices = decision.groupby(POLICY_KEYS, sort=False).indices
    c_indices = counts.groupby(POLICY_KEYS, sort=False).indices
    planned_policies = len(protocol["n_grid"]) * len(protocol["noise_grid"]) * rounds * len(protocol["policies"]) * len(protocol["budgets"])
    require(len(policies) == planned_policies, "planned policy count")
    seen, checked, infeasible = set(), 0, 0
    tail_rows, tail_distinct, zero_p = 0, 0, 0
    for fitted in policies:
        key = tuple(fitted[x] for x in POLICY_KEYS)
        require(key not in seen, "unique fitted-policy keys")
        seen.add(key)
        n, rel, round_no, policy, budget = key
        require(policy in protocol["policies"] and fitted["members"] == protocol["policies"][policy], "frozen policy members")
        require(round_no in range(rounds) and budget in protocol["budgets"]
                and n in protocol["n_grid"] and rel in protocol["noise_grid"],
                "frozen policy n, noise, round and budget")
        stratum = raw[(raw.n == n) & (raw.rel == rel)]
        ref = stratum[(stratum.phase == "reference") & (stratum.calibration_round == round_no) & stratum.screen]
        cal = stratum[(stratum.phase == "calibration") & (stratum.calibration_round == round_no) & stratum.screen]
        test = stratum[stratum.phase == "test"].copy().sort_values(TEST_KEYS).reset_index(drop=True)
        for method in fitted["members"]:
            audit_reference(ref, method, fitted["reference"][method])
        score, decline = independent_policy_score(cal, fitted["members"], fitted["reference"])
        expected_fit = independent_threshold(score, decline, budget)
        require(bool(fitted["fitted"]["feasible"]) == expected_fit["feasible"], "threshold feasibility")
        for field in ("n", "declines", "allowed_withheld"):
            require(int(fitted["fitted"][field]) == expected_fit[field], f"threshold {field}")
        if not expected_fit["feasible"]:
            require(fitted["fitted"]["threshold"] is None, "infeasible threshold is undefined")
            require(fitted["fitted"]["reason"] == ("empty_screened_H0" if not len(cal)
                    else "declines_exceed_budget"), "infeasible threshold reason")
            require(key not in d_indices and key not in c_indices, "no decisions for infeasible policy")
            infeasible += 1
            continue
        for field in ("threshold", "alarms", "withheld_rate"):
            close(fitted["fitted"][field], expected_fit[field], f"threshold {field}")
        require(key in d_indices and key in c_indices, "feasible policy has decisions and counts")
        saved = decision.iloc[d_indices[key]].sort_values(TEST_KEYS).reset_index(drop=True)
        require(len(saved) == len(test), "same test bank reused by every fitted policy")
        require(saved[TEST_KEYS].equals(test[TEST_KEYS]), "test key identity")
        require(saved.klass.equals(test.klass), "test class identity")
        require(saved.dataset_sha256.fillna("").equals(test.dataset_sha256.fillna("")), "test dataset identity")
        score, decline = independent_policy_score(test, fitted["members"], fitted["reference"])
        alarm = ~decline & (score > expected_fit["threshold"])
        released = test.screen.to_numpy() & ~alarm & ~decline
        close(saved.policy_score.to_numpy(), score, "test joint affine score")
        for name, expected in (("screen", test.screen.to_numpy()), ("alarm", alarm),
                               ("decline", decline), ("released", released)):
            require(np.array_equal(saved[name].to_numpy(), expected), f"test {name}")
            test[name] = expected
        close(saved.legacy_endpoint_error, test.legacy_endpoint_error, "test endpoint identity")
        recomputed = independently_count(test, protocol["danger_threshold"])
        actual = counts.iloc[c_indices[key]].set_index("task").sort_index()
        require(actual.index.equals(recomputed.index), "task-count inventory")
        require(np.array_equal(actual[COUNT_NAMES].to_numpy(), recomputed[COUNT_NAMES].to_numpy()),
                "all per-task event counts and missing-outcome accounting")
        # Verify values above each member's reference maximum remain ordered.
        # This is an observed diagnostic, not a requirement for positive power.
        for method in fitted["members"]:
            ref_value, ref_valid = raw_transform(ref, method)
            test_value, test_valid = raw_transform(test, method)
            if ref_valid.any():
                beyond = test_valid & (test_value > np.max(ref_value[ref_valid]))
                tail_rows += int(beyond.sum())
                tail_distinct += len(np.unique(test_value[beyond]))
            if method in "ABCD":
                zero_p += int(np.sum(test_valid & (test[f"score_{method}"].to_numpy() == 0)))
        checked += 1
    require(len(d_indices) == checked and len(c_indices) == checked, "no unregistered policy outputs")
    uncertainty = audit_uncertainty(run, manifest, raw, counts, policies)
    result = {
        "complete": True, "raw_rows": len(raw), "decision_rows": len(decision),
        "checked_policy_strata_rounds": checked, "infeasible_policy_strata_rounds": infeasible,
        "calibration_rounds": rounds, "raw_sha256": digest(run / "raw_policy.parquet"),
        "auditor_sha256": digest(Path(__file__)), "provenance": provenance,
        "uncertainty_audit": uncertainty,
        "tail_diagnostics_summed_over_fitted_policy_uses": {
            "test_rows_above_member_reference_max": tail_rows,
            "distinct_transformed_tail_values": tail_distinct,
            "zero_p_uses_clipped_to_float64_tiny": zero_p,
        },
        "checks": ["raw/cell integrity", "v1 recursive preservation", "legacy input integrity",
                   "phase/round seed and dataset separation", "H0-only reference/calibration",
                   "independent robust affine transformation", "strict threshold, ties and declines",
                   "shared independent test bank", "all decisions", "all per-task accounting",
                   "round metrics and calibration sensitivity", "three paired bootstrap distributions"],
    }
    path = run / "audit.json"
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", default="runs/p01_v2")
    arguments = parser.parse_args()
    selected = Path(arguments.run)
    audit(selected if selected.is_absolute() else ROOT / selected)
