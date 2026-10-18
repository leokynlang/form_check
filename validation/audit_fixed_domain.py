"""Independent P0-2 frozen-domain, provenance, raw-outcome and policy audit.

No experimental implementation is imported. Saved expressions, arrays and raw
records are independently interpreted with NumPy/SciPy/SymPy. This audit checks
analytic-argument presence, never represents finite grid checks as a proof of
positivity on the whole compact domain.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import sympy as sp
from scipy import special

ROOT = Path(__file__).resolve().parents[1]
METHODS = tuple("ABCDEFG")
ERROR_FIELDS = ("weighted_mean_relative_error", "weighted_p95_relative_error", "max_relative_error")
COUNT_NAMES = ["total", "released", "h0_total", "h0_screened", "h0_released",
               "h0_withheld_screened", "h0_alarm_screened", "h0_decline_screened",
               "h1_total", "h1_alarm", "h1_decline", "h1_danger_total",
               "h1_danger_screened", "h1_danger_released", "h1_safe_screened",
               "h1_safe_released", "all_danger_released", "h0_danger_released",
               "released_known_outcome", "outcome_unavailable", "released_unknown_outcome",
               "data_failures", "fit_failures"]
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


def require(condition, message):
    if not bool(condition):
        raise AssertionError(message)


def close(actual, expected, message, rtol=3e-11, atol=1e-14):
    require(np.allclose(actual, expected, rtol=rtol, atol=atol, equal_nan=True), message)


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def array_digest(value):
    value = np.ascontiguousarray(value, dtype="<f8")
    header = json.dumps(list(value.shape), separators=(",", ":")).encode("ascii")
    return hashlib.sha256(header + b"\n<float64-le>\n" + value.tobytes()).hexdigest()


def tree_hashes(path):
    path = Path(path)
    return {p.relative_to(path).as_posix(): digest(p)
            for p in sorted(path.rglob("*")) if p.is_file()}


def _expression(task):
    functions = {"log1p": lambda x: sp.log(1+x), "expm1": lambda x: sp.exp(x)-1,
                 "lambertw": sp.LambertW, "abs": sp.Abs}
    return sp.sympify(task["expression"], locals=functions)


def _canonical(expr):
    """Independent exhaustive variable-renaming canonicalization."""
    expr = sp.factor_terms(expr)
    variables = sorted(expr.free_symbols, key=str)
    if variables:
        _, expr = expr.as_independent(*variables, as_Add=False)
    replacements = [sp.Symbol(f"v{j}") for j in range(len(variables))]
    return min(sp.srepr(expr.xreplace(dict(zip(variables, perm, strict=True))))
               for perm in itertools.permutations(replacements))


@lru_cache(maxsize=256)
def _truth_function(expression):
    expr = sp.sympify(expression)
    # The exact exponential difference law is evaluated stably, independently
    # of any implementation adapter that produced the frozen truth arrays.
    def difference(node):
        return isinstance(node, sp.Add) and len(node.args) == 2 and sp.Integer(-1) in node.args and any(a.func == sp.exp for a in node.args)

    def replacement(node):
        value = next(a for a in node.args if a.func == sp.exp)
        return sp.Function("expm1")(value.args[0])

    expr = expr.replace(difference, replacement)
    expr = expr.replace(lambda node: node.func == sp.Max and len(node.args) == 2,
                        lambda node: sp.Function("maximum")(*node.args))
    symbols = sorted(expr.free_symbols, key=lambda s: int(str(s)[1:]))
    return sp.lambdify(symbols, expr,
                       modules=[{"expm1": np.expm1, "maximum": np.maximum,
                                 "erf": special.erf, "besseli": special.iv,
                                 "elliptic_k": special.ellipk, "LambertW": special.lambertw},
                                "scipy", "numpy"], cse=False)


def independent_truth(task, X):
    value = _truth_function(task["expression"])(*[X[:, j] for j in range(X.shape[1])])
    value = np.asarray(value)
    if np.iscomplexobj(value):
        require(np.all(value.imag == 0), f"{task['task']}: complex independent truth")
        value = value.real
    return np.broadcast_to(value, (len(X),)).astype(float)


def independent_registry(tasks):
    names = [t["task"] for t in tasks]
    require(len(names) == len(set(names)), "Duplicate registry task identities")
    old_pool = read_json(ROOT / "runs/srsd_full/config.json")["task_pool"]
    by_unit, by_canonical = defaultdict(list), defaultdict(list)
    historical = set()
    for task in tasks:
        name, klass, split = task["task"], task["klass"], task["split"]
        require(klass in {"H0", "H1"}, f"{name}: invalid hypothesis class")
        source = task.get("source_equations", [])
        source = [source] if isinstance(source, str) else source
        hits = ({name} | set(source)) & set(old_pool)
        if hits:
            require(split == "development", f"{name}: historical alias outside development")
            historical |= hits
        unit = "design:" + task["design_group"] if klass == "H0" else "structure:" + task["structural_group"]
        by_unit[unit].append(task)
        if klass == "H1":
            by_canonical[_canonical(_expression(task))].append(task)
        p = int(task["n_features"])
        outer, train = np.asarray(task["outer_bounds"]), np.asarray(task["train_bounds"])
        require(outer.shape == train.shape == (p, 2), f"{name}: domain dimensions")
        require(np.isfinite(outer).all() and np.isfinite(train).all(), f"{name}: nonfinite bounds")
        require(np.all((outer[:, 0] > 0) & (outer[:, 0] <= train[:, 0]) & (train[:, 0] < train[:, 1]) & (train[:, 1] <= outer[:, 1])), f"{name}: positive nested bounds")
        require(len(task["coordinate_scales"]) == p, f"{name}: missing scales")
        require(all(s in {"linear", "log", "integer"} for s in task["coordinate_scales"]), f"{name}: unsupported scale")
        require(bool(task.get("validity_argument", "").strip()), f"{name}: missing analytic argument")
    for unit, members in by_unit.items():
        require(len({t["split"] for t in members}) == 1, f"Group crosses split: {unit}")
    for members in by_canonical.values():
        require(len({t["split"] for t in members}) == 1, f"Formula twins cross split: {[t['task'] for t in members]}")
    require(historical == set(old_pool), "Not all 32 historical identities retained as development")
    return {"n_tasks": len(tasks), "historical32_development": len(historical),
            "split_class_counts": {f"{split}/{klass}": int(count) for (split, klass), count in Counter((t["split"], t["klass"]) for t in tasks).items()},
            "groups": {key: {"split": members[0]["split"], "tasks": len(members)} for key, members in sorted(by_unit.items())},
            "h0_exception": "Monomial functional family is shared by design; separate frozen design groups are the resampling units.",
            "limitation": "Syntax canonicalization is not complete algebraic equivalence; manual conservative family groups remain primary."}


def _grid_entries(grid_manifest):
    entries = grid_manifest.get("tasks", grid_manifest.get("grids", grid_manifest))
    if isinstance(entries, list):
        return {e["task"]: e for e in entries}
    require(isinstance(entries, dict), "Grid manifest entries must be a mapping/list")
    return entries


def _artifact_path(run, path):
    path = Path(path)
    return path if path.is_absolute() else run / path


def audit_grids(run, tasks, grid_manifest):
    entries, grids, results = _grid_entries(grid_manifest), {}, {}
    require(set(entries) == {t["task"] for t in tasks}, "Grid manifest does not cover exactly the task registry")
    for task in tasks:
        name, entry = task["task"], entries[task["task"]]
        path_value = entry.get("file", entry.get("path", entry.get("npz_path", f"grids/{name}.npz")))
        path = _artifact_path(run, path_value)
        require(path.is_file(), f"Missing frozen grid: {path}")
        expected_file = entry.get("npz_sha256", entry.get("file_sha256", entry.get("sha256")))
        require(expected_file is not None and digest(path) == expected_file, f"{name}: NPZ file hash mismatch/missing")
        with np.load(path, allow_pickle=False) as data:
            X, weights = data["X"], data["weights"]
            truth = data["truth"]
        n, p = X.shape
        require(p == task["n_features"] and n > 0, f"{name}: grid dimensions")
        require(np.isfinite(X).all() and (X > 0).all(), f"{name}: grid inputs not finite positive")
        require(truth.shape == weights.shape == (n,), f"{name}: truth/weight dimensions")
        require(np.isfinite(truth).all() and (truth > 0).all(), f"{name}: frozen truth not finite positive")
        require(np.isfinite(weights).all() and (weights > 0).all(), f"{name}: invalid weights")
        close(weights.sum(), 1., f"{name}: weights must sum to one", rtol=0, atol=1e-12)
        require(len(np.unique(X, axis=0)) == n, f"{name}: repeated grid points")
        outer, train = np.asarray(task["outer_bounds"]), np.asarray(task["train_bounds"])
        require(((outer[:, 0] <= X) & (outer[:, 1] >= X)).all(), f"{name}: point outside valid box")
        require(np.any((train[:, 0] > X) | (train[:, 1] < X), axis=1).all(), f"{name}: interpolation point in external grid")
        for j, scale in enumerate(task["coordinate_scales"]):
            if scale == "integer":
                require(np.equal(X[:, j], np.rint(X[:, j])).all(), f"{name}: integer coordinate violation")
        labels = entry["directions"]
        require(len(labels) == n and all(label for label in labels), f"{name}: missing grid labels")
        axes = sorted({label for row in labels for label in row if label.startswith("axis_")})
        require(bool(axes), f"{name}: missing axial support")
        expected_weights = np.zeros(n)
        for label in axes:
            membership = np.asarray([label in row for row in labels])
            expected_weights[membership] += .5 / len(axes) / membership.sum()
            _, variable, side = label.split("_")
            index = int(variable[1:])-1
            require((X[membership, index] < train[index, 0]).all() if side == "lower" else (X[membership, index] > train[index, 1]).all(), f"{name}: axial direction label mismatch")
        joint = np.asarray(["joint" in row for row in labels])
        require(joint.any(), f"{name}: missing joint support")
        expected_weights[joint] += .5 / joint.sum()
        close(weights, expected_weights, f"{name}: independently reconstructed axial/joint mass", rtol=0, atol=1e-14)
        actual_truth = independent_truth(task, X)
        close(truth, actual_truth, f"{name}: independently evaluated truth", rtol=2e-11, atol=0)
        metadata = entry.get("validation", entry)
        for field, array in [("coordinates_sha256", X), ("weights_sha256", weights), ("truth_sha256", truth)]:
            require(metadata.get(field) == array_digest(array), f"{name}: {field}")
        grids[name] = (X, weights, truth)
        results[name] = {"points": n, "features": p, "npz_sha256": digest(path),
                         "truth_max_relative_difference": float(np.max(np.abs(actual_truth/truth-1))),
                         "grid_coordinates_sha256": array_digest(X)}
    return grids, results


def audit_provenance(run, manifest):
    protocol = manifest["protocol"]
    protocol_path = run / "protocol_frozen.json"
    require(protocol_path.exists() and read_json(protocol_path) == protocol, "Frozen protocol JSON differs from manifest")
    require(digest(protocol_path) == manifest["protocol_sha256"], "Frozen protocol digest mismatch")
    fingerprint = {key: value for key, value in manifest.items() if key != "fingerprint"}
    require(hashlib.sha256(json.dumps(fingerprint, sort_keys=True).encode()).hexdigest() == manifest["fingerprint"], "Root manifest fingerprint mismatch")
    require(manifest["freeze_before_study_sampling"] is True, "Missing before-sampling freeze declaration")
    source = manifest.get("source_sha256", manifest.get("source_hashes"))
    require(isinstance(source, dict) and source, "Missing source hash inventory")
    source_count = 0
    for name, value in source.items():
        expected = value.get("sha256") if isinstance(value, dict) else value
        require(digest(_artifact_path(ROOT, name)) == expected, f"Frozen source changed: {name}")
        source_count += 1
    snapshot_path = run / "source_snapshot.json"
    require(snapshot_path.exists(), "Missing frozen source snapshot inventory")
    snapshots = read_json(snapshot_path)["files"]
    require(set(snapshots) == set(source), "Source snapshot inventory differs from frozen source inventory")
    for name, record in snapshots.items():
        require(digest(run / record["snapshot"]) == record["sha256"], f"Source snapshot changed: {name}")
        expected_source = source[name].get("sha256") if isinstance(source[name], dict) else source[name]
        require(record["sha256"] == expected_source, f"Source snapshot/source digest disagreement: {name}")
    preservation_path = run / "preservation_before.json"
    preservation = read_json(preservation_path) if preservation_path.exists() else manifest.get("preservation_before", manifest.get("preservation"))
    require(isinstance(preservation, dict), "Missing frozen v1/v2 preservation inventory")
    preserved_files = {}
    for name, record in preservation.items():
        values = record["files"] if isinstance(record, dict) and "files" in record else record
        if isinstance(values, dict) and all(isinstance(v, str) and len(v) == 64 for v in values.values()):
            path = _artifact_path(ROOT, name)
            require(tree_hashes(path) == values, f"Preserved directory changed: {name}")
            preserved_files[name] = len(values)
    require(len(preserved_files) >= 2, "Both v1 and v2 must have complete tree-hash preservation records")
    for field, path in [("task_registry_sha256", run / "task_registry.json"),
                        ("grid_manifest_sha256", run / "grid_manifest.json")]:
        require(manifest.get(field) == digest(path), f"Missing/mismatched frozen {field}")
    return {"source_files": source_count, "snapshot_files": len(snapshots), "preserved_files": preserved_files,
            "protocol_id": protocol["protocol_id"]}


def audit_raw(stage_dir, tasks, protocol, grids, grid_results):
    path = stage_dir / "raw.parquet"
    raw = pd.read_parquet(path)
    registry = {t["task"]: t for t in tasks}
    rep_name = "rep" if "rep" in raw else "replicate"
    keys = ["phase", "calibration_round", "task", "n", "rel", rep_name]
    require(not raw.duplicated(keys).any(), "Duplicate raw experiment identities")
    seed_name = "seed" if "seed" in raw else "base_seed"
    require(raw[seed_name].nunique() == len(raw), "Seed collision across frozen experimental identities")
    require(set(raw.task) <= set(registry), "Unknown task in raw outcomes")
    require(not raw.split.eq("reserve").any(), "Reserved task has noisy data or gate outcomes")
    require(raw.phase.isin(["reference", "calibration", "test"]).all(), "Unknown raw phase")
    for phase in ("reference", "calibration"):
        part = raw[raw.phase.eq(phase)]
        require(part.klass.eq("H0").all() and part.split.eq("calibration").all(), f"{phase}: wrong design pool")
        require(set(part.calibration_round) == set(range(protocol["calibration_rounds"])), f"{phase}: missing calibration rounds")
    test = raw[raw.phase.eq("test")]
    require(test.split.isin(["development", "confirmation", "stress"]).all(), "Test phase contains forbidden split")
    require(test.calibration_round.eq(-1).all(), "Test bank must be shared rather than regenerated across calibration rounds")
    stage_grid = protocol["stages"][stage_dir.name]
    require(set(raw.n) == set(stage_grid["n_grid"]) and set(raw.rel) == set(stage_grid["noise_grid"]), "Raw n/noise grid mismatch")
    planned = {}
    for task in tasks:
        for phase in ("reference", "calibration", "test"):
            allowed = task["split"] == "calibration" if phase != "test" else task["split"] in {"development", "confirmation", "stress"}
            if not allowed:
                continue
            for rnd in (range(protocol["calibration_rounds"]) if phase != "test" else [-1]):
                for n in stage_grid["n_grid"]:
                    for rel in stage_grid["noise_grid"]:
                        if stage_dir.name == "stage2" and n in protocol["stages"]["stage1"]["n_grid"] and rel in protocol["stages"]["stage1"]["noise_grid"]:
                            continue
                        planned[(phase, rnd, task["task"], n, rel)] = protocol[f"{phase}_replicates"]
    observed = raw.groupby(keys[:-1], sort=False).size().to_dict()
    require(observed == planned, "Raw cell inventory/counts differ from the complete frozen plan")
    for key, part in raw.groupby(keys[:-1], sort=False):
        require(set(part[rep_name]) == set(range(planned[key])), f"Noncontiguous/reused replicates in cell {key}")
    phase_codes = {"reference": 21, "calibration": 22, "test": 23}
    task_words = {name: [int(v) for v in np.frombuffer(hashlib.sha256(name.encode()).digest()[:8], dtype="<u4")] for name in registry}
    for phase, rnd, name, n, rel, rep, actual_seed in raw[[*keys, seed_name]].itertuples(index=False, name=None):
        entropy = [protocol["master_seed"], 3, phase_codes[phase], rnd+1, *task_words[name], n, round(rel*1_000_000), rep]
        expected_seed = int(np.random.SeedSequence(entropy).generate_state(1, dtype=np.uint64)[0]) % 2**62 + 2**40
        require(actual_seed == expected_seed, f"Independent seed mismatch: {(phase, rnd, name, n, rel, rep)}")
    require(not raw.dataset_sha256.dropna().duplicated().any(), "Dataset hashes reused inside P0-2")
    history_overlap = {}
    for old in ("p01_pilot_v1", "p01_v2"):
        prior = pd.read_parquet(ROOT / f"runs/{old}/raw_policy.parquet", columns=["base_seed", "dataset_sha256"])
        overlap = {"seed_overlap": len(set(prior.base_seed) & set(raw[seed_name])),
                   "dataset_hash_overlap": len(set(prior.dataset_sha256.dropna()) & set(raw.dataset_sha256.dropna()))}
        require(overlap == {"seed_overlap": 0, "dataset_hash_overlap": 0}, f"Old data reused from {old}")
        history_overlap[old] = overlap
    legacy = set(pd.read_parquet(ROOT / "runs/srsd_full/raw_srsd.parquet", columns=["base_seed"]).base_seed)
    legacy.update(pd.read_csv(ROOT / "runs/gate_cmp/results.csv", usecols=["seed"]).seed)
    require(not (legacy & set(raw[seed_name])), "P0-2 seed overlap with original S2/S5 historical experiments")
    expected_rows, audited_outcomes = 0, 0
    max_differences = dict.fromkeys(ERROR_FIELDS, 0.)
    for name, frame in raw.groupby("task", sort=True):
        task = registry[name]
        for field in ("klass", "split", "structural_group", "source_kind", "regime_kind"):
            require(frame[field].eq(task[field]).all(), f"{name}: raw {field} differs from frozen registry")
        unit = "design:" + task["design_group"] if task["klass"] == "H0" else "structure:" + task["structural_group"]
        require(frame.group_unit.eq(unit).all(), f"{name}: raw bootstrap unit mismatch")
        require(frame.grid_sha256.eq(grid_results[name]["grid_coordinates_sha256"]).all(), f"{name}: raw frozen grid hash")
        X, weights, truth = grids[name]
        successful = frame[frame.fit_failure.isna() & frame.data_failure.isna() & frame.coef_json.notna()]
        missing_coefficient = frame.fit_failure.isna() & frame.data_failure.isna() & frame.coef_json.isna()
        require(frame.loc[missing_coefficient, "outcome_failure"].notna().all(), f"{name}: unexplained absent fitted coefficients")
        if len(successful):
            coefficients = np.asarray([json.loads(v) for v in successful.coef_json], dtype=float)
            require(coefficients.shape == (len(successful), X.shape[1] + 1), f"{name}: coefficient shape")
            require(np.isfinite(coefficients).all(), f"{name}: nonfinite stored coefficients")
            design = np.column_stack([np.ones(len(X)), np.log(X)])
            with np.errstate(over="ignore", invalid="ignore"):
                errors = np.abs(np.expm1(coefficients @ design.T - np.log(truth)[None, :]))
            finite = np.isfinite(errors).all(axis=1)
            predicted = {
                "weighted_mean_relative_error": np.full(len(successful), np.nan),
                "weighted_p95_relative_error": np.full(len(successful), np.nan),
                "max_relative_error": np.full(len(successful), np.nan),
            }
            if finite.any():
                valid_errors = errors[finite]
                order = np.argsort(valid_errors, axis=1, kind="stable")
                cumulative = np.cumsum(weights[order], axis=1)
                p95_positions = np.minimum(np.sum(cumulative < .95, axis=1), len(weights)-1)
                predicted[ERROR_FIELDS[0]][finite] = valid_errors @ weights
                predicted[ERROR_FIELDS[1]][finite] = np.take_along_axis(valid_errors, order, axis=1)[np.arange(finite.sum()), p95_positions]
                predicted[ERROR_FIELDS[2]][finite] = valid_errors.max(axis=1)
            for field, expected in predicted.items():
                actual = successful[field].to_numpy(dtype=float)
                close(actual, expected, f"{name}: independently recomputed {field}", rtol=3e-9, atol=2e-13)
                if finite.any():
                    max_differences[field] = max(max_differences[field], float(np.max(np.abs(actual[finite]-expected[finite]))))
            audited_outcomes += len(successful)
    expected_rows = sum(planned.values())
    require(len(raw) == expected_rows, f"Raw count {len(raw)} differs from planned count {expected_rows}")
    expected_tasks = {t["task"] for t in tasks if t["split"] != "reserve"}
    require(set(raw.task) == expected_tasks, "Raw task inventory differs from complete planned nonreserve set")
    stage_manifest = read_json(stage_dir / "manifest.json")
    require(stage_manifest["complete"] and stage_manifest["raw_rows"] == len(raw), "Stage manifest not complete or wrong row count")
    require(stage_manifest["expected_cells"] == len(planned) and stage_manifest["expected_datasets"] == expected_rows, "Stage manifest plan counts mismatch")
    require(stage_manifest["raw_sha256"] == digest(path), "Stage raw digest mismatch")
    for field in ("data_failures", "fit_failures", "outcome_failures"):
        require(stage_manifest[field] == int(raw[field[:-1]].notna().sum()), f"Stage failure count: {field}")
    return raw, {"rows": len(raw), "test_rows": len(test), "independently_recomputed_grid_outcomes": audited_outcomes,
                 "maximum_absolute_error_differences": max_differences, "raw_sha256": digest(path),
                 "phase_rows": raw.phase.value_counts().to_dict(), "no_reserved_outcomes": True,
                 "unique_seeds": int(raw[seed_name].nunique()), "cells": len(planned),
                 "v1_v2_data_disjointness": history_overlap}


def raw_score(frame, method):
    value = frame[f"score_{method}"].to_numpy(dtype=float)
    valid = np.isfinite(value) & frame[f"fail_{method}"].fillna("").astype(str).str.strip().eq("").to_numpy()
    if method in "ABCD":
        p = -value
        valid &= (p >= 0) & (p <= 1)
        with np.errstate(divide="ignore", invalid="ignore"):
            value = -np.log10(np.maximum(p, np.finfo(float).tiny))
    value = value.copy()
    value[~valid] = np.nan
    return value, valid


def independent_reference(frame, method):
    value, valid = raw_score(frame, method)
    selected = value[valid]
    center = float(np.median(selected)) if len(selected) else None
    scale, scale_method = None, "unusable"
    if len(selected):
        quartiles = np.quantile(selected, [.25, .75])
        choices = [("iqr_normalized", float(np.diff(quartiles)[0]/1.3489795003921634)),
                   ("mad_normalized", float(np.median(np.abs(selected-center))*1.482602218505602)),
                   ("population_sd", float(np.std(selected, ddof=0)))]
        for name, candidate in choices:
            if np.isfinite(candidate) and candidate > 0:
                scale, scale_method = candidate, name
                break
    return {"usable": scale is not None, "center": center, "scale": scale,
            "scale_method": scale_method, "n": len(selected)}


def independent_score(frame, members, reference):
    valid = np.ones(len(frame), dtype=bool)
    columns = []
    for method in members:
        value, okay = raw_score(frame, method)
        ref = reference[method]
        if ref["usable"]:
            with np.errstate(invalid="ignore", over="ignore"):
                value = (value-ref["center"])/ref["scale"]
            okay &= np.isfinite(value)
        else:
            okay[:] = False
            value[:] = np.nan
        valid &= okay
        columns.append(value)
    joint = np.max(np.column_stack(columns), axis=1)
    joint[~valid] = np.nan
    return joint, ~valid


def independent_threshold(score, decline, budget):
    n, failures = len(score), int(np.sum(decline))
    cap = int(np.floor(budget*n + 1e-10))
    feasible = bool(n and failures <= cap)
    result = {"n": n, "declines": failures, "allowed_withheld": cap, "feasible": feasible, "threshold": None}
    if feasible:
        ordered = np.sort(score[~decline])
        result["threshold"] = float(ordered[len(ordered)-(cap-failures)-1])
        result["alarms"] = int(np.sum(ordered > result["threshold"]))
        result["withheld_rate"] = (result["alarms"]+failures)/n
    return result


def nullable(actual, expected, description):
    if expected is None:
        require(actual is None or (isinstance(actual, (float, np.floating)) and np.isnan(actual)), description)
    else:
        close(actual, expected, description)


def independent_pool(frame, pool):
    h0 = frame.klass.eq("H0") & frame.split.eq("confirmation")
    h1 = frame.klass.eq("H1") & frame.split.eq("confirmation")
    if pool == "confirmation_all":
        return h0 | h1
    if pool == "confirmation_physical_sources":
        return h0 | (h1 & frame.source_kind.eq("new_source_positive_restriction"))
    if pool == "confirmation_designed_controls":
        return h0 | (h1 & frame.source_kind.eq("synthetic_positive_mechanism"))
    if pool == "development":
        return frame.split.eq("development")
    if pool == "stress":
        return h0 | (frame.klass.eq("H1") & frame.split.eq("stress"))
    raise AssertionError(f"Unregistered pool {pool}")


def independent_counts(frame, alarm, decline, released, outcome, threshold):
    h0, h1 = frame.klass.eq("H0").to_numpy(), frame.klass.eq("H1").to_numpy()
    screen = frame.screen.to_numpy()
    known = np.isfinite(frame[outcome].to_numpy(dtype=float))
    danger = known & (frame[outcome].to_numpy(dtype=float) > threshold)
    safe = known & ~danger
    events = [np.ones(len(frame), dtype=bool), released, h0, h0 & screen, h0 & released,
              h0 & screen & (alarm | decline), h0 & screen & alarm, h0 & screen & decline,
              h1, h1 & alarm, h1 & decline, h1 & danger, h1 & danger & screen,
              h1 & danger & released, h1 & safe & screen, h1 & safe & released,
              danger & released, h0 & danger & released, released & known, ~known,
              released & ~known, frame.data_failure.notna().to_numpy(), frame.fit_failure.notna().to_numpy()]
    return pd.DataFrame({name: value.astype(np.int64) for name, value in zip(COUNT_NAMES, events, strict=True)}, index=frame.index)


def divide(num, den):
    return np.divide(num, den, out=np.full(np.broadcast_shapes(np.shape(num), np.shape(den)), np.nan), where=np.asarray(den) > 0)


def independent_draws(classes, rounds, repetitions, seed):
    rng = np.random.default_rng(seed)
    weights = np.zeros((repetitions, len(classes)), dtype=int)
    for klass in ("H0", "H1"):
        ids = np.flatnonzero(np.asarray(classes) == klass)
        sampled = rng.choice(ids, size=(repetitions, len(ids)), replace=True)
        for position in ids:
            weights[:, position] = np.count_nonzero(sampled == position, axis=1)
    return weights, rng.integers(0, rounds, size=(repetitions, rounds))


def audit_stats(saved, original, weighted, round_draws, feasible, description):
    finite = np.isfinite(original)
    values = original[finite]
    defined = bool(np.all(feasible) and finite.all())
    require(saved["primary_defined"] == defined, description+": definition")
    nullable(saved["point_mean"], float(original.mean()) if defined else None, description+": point")
    nullable(saved["feasible_defined_round_only_sensitivity_mean"], float(values.mean()) if len(values) else None, description+": sensitivity point")
    sensitivity = saved["calibration_sensitivity"]
    require(sensitivity["defined_rounds"] == int(finite.sum()), description+": defined rounds")
    for name, expected in [("sample_sd", float(np.std(values, ddof=1)) if len(values)>1 else None),
                           ("min", float(values.min()) if len(values) else None),
                           ("max", float(values.max()) if len(values) else None),
                           ("empirical_percentile_range95", np.quantile(values, [.025,.975]) if len(values) else None)]:
        nullable(sensitivity[name], expected, description+": "+name)
    if not defined:
        require(all(saved[field] is None for field in ("conditional_group", "round_mean", "crossed")), description+": undefined intervals")
        return
    distributions = {"conditional_group": weighted.mean(axis=1), "round_mean": original[round_draws].mean(axis=1),
                     "crossed": np.take_along_axis(weighted, round_draws, axis=1).mean(axis=1)}
    for field, draws in distributions.items():
        good = np.isfinite(draws)
        interval = saved[field]
        require(interval["defined_draws"] == int(good.sum()) and interval["total_draws"] == len(draws), description+": "+field+" draw counts")
        nullable(interval["ci95"], np.quantile(draws[good], [.025,.975]) if good.any() else None, description+": "+field+" interval")


def ratio_arrays(cube, numerator, denominator, weights):
    num, den = cube[:, :, COUNT_NAMES.index(numerator)], cube[:, :, COUNT_NAMES.index(denominator)]
    return divide(num.sum(axis=1), den.sum(axis=1)), divide(weights @ num.T, weights @ den.T)


def audit_analysis(stage_dir, raw, protocol):
    summary = read_json(stage_dir / "summary.json")
    require(summary["raw_sha256"] == digest(stage_dir / "raw.parquet"), "Summary raw digest")
    require(summary["analysis_sha256"] == digest(ROOT / "experiments/fixed_domain/analyze.py"), "Summary analysis source")
    require(summary["protocol_id"] == protocol["protocol_id"] and summary["stage"] == stage_dir.name, "Summary protocol/stage")
    policies = read_json(stage_dir / "policies.json")
    policy_keys = ["n", "rel", "budget", "policy", "calibration_round"]
    policy_lookup = {tuple(p[key] for key in policy_keys): p for p in policies}
    require(len(policy_lookup) == len(policies), "Repeated fitted policy identities")
    decisions = pd.read_parquet(stage_dir / "decisions.parquet")
    decision_keys = [*policy_keys, "task", "replicate"]
    require(not decisions.duplicated(decision_keys).any(), "Repeated saved decisions")
    decision_indices = decisions.groupby(policy_keys, sort=False).indices
    saved_counts = pd.read_parquet(stage_dir / "group_counts.parquet")
    count_keys = [*policy_keys, "pool", "outcome", "danger_threshold", "group_unit"]
    require(not saved_counts.duplicated(count_keys).any(), "Repeated group counts")
    counts_lookup = saved_counts.set_index(count_keys)
    raw = raw.copy()
    raw["screen"] = np.isfinite(raw.r2) & raw.r2.ge(protocol["screen_r2"])
    specs = [("weighted_mean_relative_error", float(t)) for t in protocol["danger_thresholds"]]
    specs += [(field, float(protocol["primary_danger_threshold"])) for field in protocol["outcome"]["secondary"] if field in raw]
    specs = list(dict.fromkeys(specs))
    primary = ("weighted_mean_relative_error", float(protocol["primary_danger_threshold"]))
    counts_checked, decisions_checked, policies_checked, stats_checked, paired_checked = 0, 0, 0, 0, 0
    strata = {(s["n"], s["rel"]): s for s in summary["strata"]}
    require(set(strata) == set(raw.groupby(["n", "rel"]).groups), "Summary stratum inventory")
    for (n, rel), frame in raw.groupby(["n", "rel"], sort=True):
        test = frame[frame.phase.eq("test")].sort_values(["task", "replicate"]).copy()
        rounds = int(protocol["calibration_rounds"])
        references, calibration = {}, {}
        for rnd in range(rounds):
            bank = frame[frame.phase.eq("reference") & frame.calibration_round.eq(rnd) & frame.screen]
            references[rnd] = {method: independent_reference(bank, method) for method in protocol["methods"]}
            calibration[rnd] = frame[frame.phase.eq("calibration") & frame.calibration_round.eq(rnd) & frame.screen]
        pool_info = {}
        for pool in protocol["evaluation_pools"]:
            selected = test[independent_pool(test, pool)]
            units = sorted(selected.group_unit.unique())
            classes = [selected.loc[selected.group_unit.eq(unit), "klass"].iloc[0] for unit in units]
            seed = protocol["master_seed"] + int(n) + round(float(rel)*1e6) + int.from_bytes(hashlib.sha256(pool.encode()).digest()[:4], "little")
            weights, round_draws = independent_draws(classes, rounds, int(protocol["group_inference"]["bootstrap_repetitions"]), seed)
            saved = strata[(n, rel)]["pools"][pool]
            require(saved["group_units"] == units and saved["test_groups"] == len(units), f"{pool}: group inventory")
            require(saved["test_tasks"] == selected.task.nunique() and saved["test_rows"] == len(selected), f"{pool}: pool population")
            pool_info[pool] = (units, weights, round_draws, independent_pool(test, pool).to_numpy())
        for budget in protocol["budgets"]:
            cubes, feasible_maps = {}, {}
            for policy, members in protocol["policies"].items():
                collections = {pool: [] for pool in pool_info}
                feasible = []
                for rnd in range(rounds):
                    key = (n, rel, float(budget), policy, rnd)
                    saved_policy = policy_lookup[key]
                    require(saved_policy["members"] == members, "Frozen policy membership")
                    for method in members:
                        for field, expected in references[rnd][method].items():
                            actual = saved_policy["reference"][method][field]
                            if isinstance(expected, (bool, str)):
                                require(actual == expected, f"{key}: reference {method}/{field}")
                            else:
                                nullable(actual, expected, f"{key}: reference {method}/{field}")
                    cal_score, cal_decline = independent_score(calibration[rnd], members, references[rnd])
                    fitted = independent_threshold(cal_score, cal_decline, float(budget))
                    for field, expected in fitted.items():
                        nullable(saved_policy["fitted"][field], expected, f"{key}: threshold {field}")
                    policies_checked += 1
                    feasible.append(fitted["feasible"])
                    if fitted["feasible"]:
                        score, decline = independent_score(test, members, references[rnd])
                        alarm = ~decline & (score > fitted["threshold"])
                        released = test.screen.to_numpy() & ~alarm & ~decline
                        saved = decisions.iloc[decision_indices[key]].sort_values(["task", "replicate"])
                        require(saved.task.tolist() == test.task.tolist() and saved.replicate.tolist() == test.replicate.tolist(), f"{key}: decision population")
                        require(saved.dataset_sha256.fillna("").tolist() == test.dataset_sha256.fillna("").tolist(), f"{key}: decision data identity")
                        close(saved.policy_score, score, f"{key}: policy scores")
                        for field, expected in [("screen", test.screen), ("alarm", alarm), ("decline", decline), ("released", released)]:
                            require(np.array_equal(saved[field].to_numpy(), np.asarray(expected)), f"{key}: {field}")
                        decisions_checked += len(saved)
                    else:
                        require(key not in decision_indices, f"{key}: infeasible policy has saved decisions")
                    for spec in specs:
                        events = independent_counts(test, alarm, decline, released, *spec) if fitted["feasible"] else None
                        for pool, (units, _, _, mask) in pool_info.items():
                            if fitted["feasible"]:
                                grouped = events.loc[mask].groupby(test.loc[mask, "group_unit"], sort=True).sum().reindex(units, fill_value=0)
                                values = grouped[COUNT_NAMES].to_numpy(dtype=float)
                            else:
                                values = np.full((len(units), len(COUNT_NAMES)), np.nan)
                            saved = counts_lookup.loc[[(*key, pool, spec[0], spec[1], unit) for unit in units]]
                            close(saved[COUNT_NAMES], values, f"{key}/{pool}/{spec}: independent group event counts", rtol=0, atol=0)
                            require(saved.feasible.eq(fitted["feasible"]).all(), f"{key}: group count feasibility")
                            counts_checked += len(units)
                            if spec == primary:
                                collections[pool].append(values)
                cubes[policy] = {pool: np.stack(values) for pool, values in collections.items()}
                feasible_maps[policy] = np.asarray(feasible, dtype=bool)
            for pool, (_, weights, round_draws, _) in pool_info.items():
                saved_budget = strata[(n, rel)]["pools"][pool]["budgets"][str(float(budget))]
                for policy in protocol["policies"]:
                    cube = cubes[policy][pool]
                    feasible = feasible_maps[policy]
                    saved = saved_budget["policies"][policy]
                    require(saved["feasible_rounds"] == int(feasible.sum()) and saved["total_rounds"] == rounds, "Summary feasibility")
                    for metric, (num, den) in RATIOS.items():
                        original, weighted = ratio_arrays(cube, num, den, weights)
                        stats = saved["metrics"][metric]
                        audit_stats(stats, original, weighted, round_draws, feasible, f"{n}/{pool}/{budget}/{policy}/{metric}")
                        denominator = cube[:, :, COUNT_NAMES.index(den)]
                        require(stats["denominator_contributing_test_groups_by_round"] == np.sum(denominator > 0, axis=1).tolist(), "Contributing groups")
                        if metric == "H0_withholding_given_screen":
                            require(stats["test_rounds_above_nominal"] == int(np.sum(np.isfinite(original) & (original > budget))), "Independent-test threshold transfer")
                            nullable(stats["test_fraction_rounds_above_nominal"], float(np.mean(original > budget)) if feasible.all() and np.isfinite(original).all() else None, "Independent-test nominal exceedance fraction")
                        stats_checked += 1
                    numerator, denominator = RATIOS["dangerous_H1_release_given_screen"]
                    group_ratios = divide(cube[:, :, COUNT_NAMES.index(numerator)], cube[:, :, COUNT_NAMES.index(denominator)])
                    defined = np.isfinite(group_ratios)
                    original = divide(np.nansum(group_ratios, axis=1), defined.sum(axis=1))
                    weighted = divide(weights @ np.nan_to_num(group_ratios, nan=0).T, weights @ defined.astype(float).T)
                    equal = saved["equal_group_dangerous_release"]
                    audit_stats(equal, original, weighted, round_draws, feasible, "Equal-group risk")
                    require(equal["defined_groups_by_round"] == defined.sum(axis=1).tolist() and equal["omitted_groups_by_round"] == (~defined).sum(axis=1).tolist(), "Equal-group omitted counts")
                for label, saved_pair in saved_budget["paired"].items():
                    first, second = label.split("_minus_")
                    both = feasible_maps[first] & feasible_maps[second]
                    for metric, saved_stats in saved_pair.items():
                        num, den = RATIOS[metric]
                        first_original, first_weighted = ratio_arrays(cubes[first][pool], num, den, weights)
                        second_original, second_weighted = ratio_arrays(cubes[second][pool], num, den, weights)
                        audit_stats(saved_stats, first_original-second_original, first_weighted-second_weighted, round_draws, both, f"Paired {label}/{metric}")
                        paired_checked += 1
    require(policies_checked == len(policies) and decisions_checked == len(decisions) and counts_checked == len(saved_counts), "Unreviewed extra policy/decision/group-count records")
    return {"status": "passed", "policy_configurations": policies_checked, "decisions": decisions_checked,
            "independent_group_count_rows": counts_checked, "primary_metric_uncertainty_summaries": stats_checked,
            "paired_uncertainty_summaries": paired_checked,
            "reconstructed": ["H0 reference transforms/scales", "budgeted thresholds", "all decisions", "all outcomes/thresholds/pools group event counts", "primary points and calibration sensitivity", "conditional-group/round/crossed intervals", "equal-group risks", "paired points and intervals"],
            "not_independently_recomputed": ["Individual A-G test-statistic numerical implementations from full sampled data", "secondary CSV flattening/formatting; underlying secondary group event counts verified"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", default="runs/p02_fixed_domain")
    parser.add_argument("--stage", default="stage1")
    parser.add_argument("--freeze-only", action="store_true")
    parser.add_argument("--raw-only", action="store_true")
    parser.add_argument("--reuse-raw-audit", action="store_true")
    args = parser.parse_args()
    require(not (args.freeze_only and (args.raw_only or args.reuse_raw_audit)), "Freeze-only cannot be combined with raw audit modes")
    require(not (args.raw_only and args.reuse_raw_audit), "Raw-only cannot reuse an earlier raw audit")
    run = _artifact_path(ROOT, args.run)
    tasks = read_json(run / "task_registry.json")
    tasks = tasks["tasks"] if isinstance(tasks, dict) else tasks
    manifest = read_json(run / "manifest.json")
    result = {"audit": "independent_p02_fixed_domain", "audit_source_sha256": digest(Path(__file__)),
              "imports_experimental_implementation": False, "raw_only": bool(args.raw_only),
              "frozen_protocol_fingerprint": manifest["fingerprint"],
              "registry": independent_registry(tasks), "provenance": audit_provenance(run, manifest)}
    grids, grid_results = audit_grids(run, tasks, read_json(run / "grid_manifest.json"))
    result["frozen_grids"] = {"tasks": len(grid_results), "points": sum(v["points"] for v in grid_results.values()),
                              "per_task": grid_results}
    stage_dir = run / args.stage
    if not args.freeze_only:
        require(read_json(stage_dir / "manifest.json")["fingerprint"] == manifest["fingerprint"], "Stage belongs to another frozen protocol")
        audited_files = ("raw.parquet",) if args.raw_only else (
            "raw.parquet", "summary.json", "policies.json", "decisions.parquet", "group_counts.parquet",
            "reference_scales.csv", "primary_metrics.csv", "sensitivity_metrics.csv", "round_metrics.csv")
        audited_hashes = {name: digest(stage_dir / name) for name in audited_files}
        if args.reuse_raw_audit:
            raw_audit_path = stage_dir / "independent_raw_audit.json"
            certificate = read_json(raw_audit_path)
            require(certificate["passed"] and certificate["raw_only"], "Saved raw-only audit is incomplete")
            require(certificate["audit_source_sha256"] == digest(Path(__file__)), "Raw audit used another audit implementation")
            require(certificate["frozen_protocol_fingerprint"] == manifest["fingerprint"], "Raw audit used another frozen protocol")
            require(certificate["raw"]["raw_sha256"] == audited_hashes["raw.parquet"], "Raw data changed since independent raw audit")
            raw = pd.read_parquet(stage_dir / "raw.parquet")
            result["raw"] = certificate["raw"]
            result["raw_audit_reused_sha256"] = digest(raw_audit_path)
        else:
            raw, result["raw"] = audit_raw(stage_dir, tasks, manifest["protocol"], grids, grid_results)
        result["analysis"] = ({"status": "not_requested_raw_only"} if args.raw_only else
                              audit_analysis(stage_dir, raw, manifest["protocol"]))
        require(all(digest(stage_dir / name) == expected for name, expected in audited_hashes.items()), "An audited raw/analysis artifact changed during the audit")
        result["audited_artifact_sha256"] = audited_hashes
    result["passed"] = True
    output = (run / "freeze_independent_audit.json" if args.freeze_only else
              stage_dir / ("independent_raw_audit.json" if args.raw_only else "independent_audit.json"))
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"passed": True, "output": str(output), "grids": len(grid_results), "raw_rows": result.get("raw", {}).get("rows")}))


if __name__ == "__main__":
    main()
