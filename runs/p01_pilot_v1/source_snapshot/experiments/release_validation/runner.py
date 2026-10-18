"""Run a frozen shared-dataset P0-1 pilot, with resumable atomic cell files."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

import srsd_data
from experiments.release_validation.gates import METHODS, run_gates
from experiments.release_validation.inventory import ROOT, digest

PACKAGE = Path(__file__).resolve().parent
PHASES = {"reference": 1, "calibration": 2, "test": 3}


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    temporary.replace(path)


def sample_seed(master, phase, task, n, rel, replicate):
    task_id = int.from_bytes(hashlib.sha256(task.encode()).digest()[:4], "little")
    seq = np.random.SeedSequence([master, PHASES[phase], task_id, n,
                                  round(rel * 1_000_000), replicate])
    # Keep parquet signed-int64 safe and outside all historical seed namespaces.
    return int(seq.generate_state(1, dtype=np.uint64)[0]) % (2**62) + 2**40


def fit_and_endpoint(data, features, truth):
    x = np.log(data[features].to_numpy())
    y = np.log(data.target.to_numpy())
    design = np.column_stack([np.ones(len(x)), x])
    coef, _, rank, _ = np.linalg.lstsq(design, y, rcond=None)
    result = {"r2": np.nan, "fit_failure": None,
              "legacy_endpoint_error": np.nan, "outcome_failure": None,
              "coef_json": json.dumps(coef.tolist())}
    if rank != design.shape[1]:
        result["fit_failure"] = "rank_deficient"
        return result
    tss = float(np.sum((y - y.mean()) ** 2))
    if not np.isfinite(tss) or tss <= 0:
        result["fit_failure"] = "degenerate_response"
        return result
    result["r2"] = 1 - float(np.sum((y - design @ coef) ** 2)) / tss
    point = np.median(x, axis=0)
    point[0] = x[:, 0].max() + np.log(3.0)
    try:
        with np.errstate(all="ignore"):
            target = float(truth(np.exp(point)))
            if not np.isfinite(target) or target <= 0:
                result["outcome_failure"] = "nonpositive_or_nonfinite_truth_at_legacy_endpoint"
                return result
            error = float(abs(np.expm1(np.r_[1.0, point] @ coef - np.log(target))))
        if not np.isfinite(error):
            result["outcome_failure"] = "nonfinite_prediction_at_legacy_endpoint"
        else:
            result["legacy_endpoint_error"] = error
    except (ValueError, OverflowError, FloatingPointError, ZeroDivisionError) as exc:
        result["outcome_failure"] = f"endpoint:{type(exc).__name__}:{exc}"
    return result


def source_record():
    paths = list(PACKAGE.glob("*.py")) + [
        ROOT / "common.py", ROOT / "gate_compare.py", ROOT / "srsd_data.py",
        ROOT / "_bootstrap.py", ROOT / "dimreflex/discover.py", ROOT / "dimreflex/_ols.py",
        ROOT / "paper/revision_loc/rerun_cd_svd.py", ROOT / "runs/srsd_full/config.json",
    ]
    upstream = Path(srsd_data._UPSTREAM)
    paths += sorted((upstream / "datasets").rglob("*.py"))
    return {str(path): digest(path) for path in paths}


def load_registry():
    # Use exactly the historical registered task set without reclassifying or
    # screening equations based on the outcomes of this pilot.
    config = json.loads((ROOT / "runs/srsd_full/config.json").read_text(encoding="utf-8"))
    classes = srsd_data.load_feynman()
    by_name = {cls()._eq_name: cls for cls in classes.values()}
    registry = []
    for name, meta in sorted(config["task_pool"].items()):
        cls = by_name[name]
        inst = cls()
        if inst.get_var_count() != meta["n_vars"]:
            raise ValueError(f"Source variable count drift: {name}")
        sampling = [[float(obj.min_value), float(obj.max_value), type(obj).__name__,
                     bool(getattr(obj, "uses_positive", True)),
                     bool(getattr(obj, "uses_negative", True))] for obj in inst.sampling_objs]
        if sampling != meta["sampling"]:
            raise ValueError(f"Source sampling drift: {name}")
        registry.append({"task": name, "klass": meta["klass"], "cls": cls,
                         "truth": inst.eq_func, "formula": str(inst.sympy_eq),
                         "sampling": sampling, "fn": srsd_data.make_dataset_fn(cls)})
    return registry


def assert_dataset_isolation(frame):
    keys = ["phase", "task", "n", "rel", "replicate"]
    if frame.duplicated(keys).any() or frame.base_seed.duplicated().any():
        raise AssertionError("Repeated run key or cross-phase/task seed")
    hashes = frame.dataset_sha256.dropna()
    if hashes.duplicated().any():
        raise AssertionError("A dataset was reused across run records")
    if not frame.loc[frame.phase.ne("test"), "klass"].eq("H0").all():
        raise AssertionError("H1 leaked into reference/calibration")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", default=str(PACKAGE / "protocol_pilot.json"))
    parser.add_argument("--out", default="runs/p01_pilot_v1")
    args = parser.parse_args()
    protocol = json.loads(Path(args.protocol).read_text(encoding="utf-8"))
    # Fail closed if a future protocol asks for unsupported knobs.
    if (protocol["permutations"] != 999 or protocol["cv_folds"] != 5
            or protocol["rbf"] != {"landmarks": 32, "bandwidth_squared": "number_of_features",
                                    "ridge_per_observation": 0.01}
            or protocol["methods"] != list(METHODS)):
        raise ValueError("Protocol requires an adapter implementation update")
    out = ROOT / args.out
    # Only create new dedicated run directories under runs, never an R2 path.
    resolved = out.resolve()
    if resolved.parent != (ROOT / "runs").resolve() or not resolved.name.startswith("p01_"):
        raise ValueError("Output must be a dedicated runs/p01_* directory")
    registry = load_registry()
    source = source_record()
    frozen = {"protocol": protocol, "source_sha256": source,
              "environment": {"python": platform.python_version(), "numpy": np.__version__,
                              "pandas": pd.__version__, "scipy": scipy.__version__},
              "registry": [{k: rec[k] for k in ("task", "klass", "formula", "sampling")}
                           for rec in registry]}
    fingerprint = hashlib.sha256(json.dumps(frozen, sort_keys=True).encode()).hexdigest()
    out.mkdir(parents=True, exist_ok=True)
    cells_dir = out / "cells"
    cells_dir.mkdir(exist_ok=True)
    manifest_path = out / "manifest.json"
    if manifest_path.exists():
        old = json.loads(manifest_path.read_text(encoding="utf-8"))
        if old["fingerprint"] != fingerprint:
            raise ValueError("Resume rejected: protocol/source/environment changed; use a new output directory")
    manifest = {**frozen, "fingerprint": fingerprint, "complete": False}
    write_json(manifest_path, manifest)
    planned = []
    for phase in PHASES:
        for task_no, rec in enumerate(registry):
            if phase != "test" and rec["klass"] != "H0":
                continue
            for n in protocol["n_grid"]:
                for rel in protocol["noise_grid"]:
                    planned.append((phase, task_no, rec, n, rel))
    started = time.perf_counter()
    all_frames = []
    for cell_no, (phase, task_no, rec, n, rel) in enumerate(planned, 1):
        path = cells_dir / f"{phase}_{task_no:03d}_n{n}_r{round(rel*1e6)}.parquet"
        expected = protocol["replicates"][phase]
        if path.exists():
            frame = pd.read_parquet(path)
            expected_replicates = set(range(expected))
            if (len(frame) != expected or set(frame.replicate) != expected_replicates
                    or not frame.phase.eq(phase).all() or not frame.task.eq(rec["task"]).all()
                    or not frame.n.eq(n).all() or not frame.rel.eq(rel).all()):
                raise ValueError(f"Checkpoint schema/count mismatch: {path}")
        else:
            rows = []
            for replicate in range(expected):
                seed = sample_seed(protocol["master_seed"], phase, rec["task"], n, rel, replicate)
                row = {"phase": phase, "task": rec["task"], "klass": rec["klass"],
                       "n": n, "rel": rel, "replicate": replicate, "base_seed": seed,
                       "data_failure": None, "dataset_sha256": None, "r2": np.nan,
                       "fit_failure": None, "outcome_failure": None, "legacy_endpoint_error": np.nan}
                for method in METHODS:
                    row.update({f"score_{method}": np.nan, f"fail_{method}": "data_unavailable"})
                try:
                    data, features = rec["fn"](seed, n, rel)
                    features = features.split(",")
                except (ValueError, OverflowError, FloatingPointError) as exc:
                    row["data_failure"] = f"generation:{type(exc).__name__}:{exc}"
                    rows.append(row)
                    continue
                row.update(run_gates(data, features, seed))
                row.update(fit_and_endpoint(data, features, rec["truth"]))
                rows.append(row)
            frame = pd.DataFrame(rows)
            temporary = path.with_suffix(".parquet.tmp")
            frame.to_parquet(temporary, index=False)
            temporary.replace(path)
        all_frames.append(frame)
        progress = {"completed_cells": cell_no, "total_cells": len(planned),
                    "last_cell": {"phase": phase, "task": rec["task"], "n": n, "rel": rel},
                    "elapsed_seconds_this_invocation": time.perf_counter() - started}
        write_json(out / "progress.json", progress)
        print(f"[{cell_no}/{len(planned)}] {phase} {rec['task']} n={n} rows={len(frame)}", flush=True)
    frame = pd.concat(all_frames, ignore_index=True)
    assert_dataset_isolation(frame)
    destination = out / "raw_policy.parquet"
    temporary = destination.with_suffix(".parquet.tmp")
    frame.to_parquet(temporary, index=False)
    temporary.replace(destination)
    historical_seeds = set(pd.read_parquet(ROOT / "runs/srsd_full/raw_srsd.parquet",
                                            columns=["base_seed"]).base_seed)
    historical_seeds.update(pd.read_csv(ROOT / "runs/gate_cmp/results.csv", usecols=["seed"]).seed)
    if historical_seeds & set(frame.base_seed):
        raise AssertionError("New sample seed overlaps a historical run")
    manifest.update({"complete": True, "rows": len(frame), "raw_sha256": digest(destination),
                     "cells": {str(p.relative_to(out)): digest(p) for p in sorted(cells_dir.glob("*.parquet"))},
                     "phase_rows": frame.phase.value_counts().to_dict(),
                     "historical_seed_overlap": 0,
                     "elapsed_seconds_this_invocation": time.perf_counter() - started})
    write_json(manifest_path, manifest)
    print(json.dumps({"complete": True, "rows": len(frame), "output": str(out)}), flush=True)


if __name__ == "__main__":
    main()
