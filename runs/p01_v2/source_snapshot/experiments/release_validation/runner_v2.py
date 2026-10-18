"""Frozen v2 new data generation with independent repeated calibration rounds."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

import srsd_data
from experiments.release_validation.gates import METHODS, run_gates
from experiments.release_validation.inventory import INPUTS, ROOT, digest
from experiments.release_validation.runner import (
    fit_and_endpoint,
    load_registry,
    write_json,
)

PACKAGE = Path(__file__).resolve().parent
V1 = ROOT / "runs/p01_pilot_v1"
PHASES = {"reference": 11, "calibration": 12, "test": 13}
_WORKER_REGISTRY = None


def sample_seed(master, phase, calibration_round, task, n, rel, replicate):
    task_id = np.frombuffer(hashlib.sha256(task.encode()).digest()[:8], dtype="<u4")
    entropy = [master, 2, PHASES[phase], calibration_round + 1,
               *(int(x) for x in task_id), n, round(rel * 1_000_000), replicate]
    return int(np.random.SeedSequence(entropy).generate_state(1, dtype=np.uint64)[0]) % (2**62) + 2**40


def tree_hashes(root):
    return {p.relative_to(root).as_posix(): digest(p) for p in sorted(root.rglob("*")) if p.is_file()}


def source_hashes():
    paths = [*PACKAGE.glob("*.py"), PACKAGE / "protocol_v2.json",
             ROOT / "common.py", ROOT / "gate_compare.py", ROOT / "srsd_data.py",
             ROOT / "_bootstrap.py", ROOT / "dimreflex/discover.py", ROOT / "dimreflex/_ols.py",
             ROOT / "paper/revision_loc/rerun_cd_svd.py", ROOT / "runs/srsd_full/config.json"]
    paths += sorted((Path(srsd_data._UPSTREAM) / "datasets").rglob("*.py"))
    return {str(p): digest(p) for p in paths}


def freeze_sources(out, source):
    saved = {}
    upstream = Path(srsd_data._UPSTREAM)
    for name, expected in source.items():
        path = Path(name)
        if digest(path) != expected:
            raise ValueError(f"Source changed during freeze: {path}")
        relative = (path.relative_to(ROOT) if path.is_relative_to(ROOT)
                    else Path("_upstream") / path.relative_to(upstream))
        target = out / "source_snapshot" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and digest(target) != expected:
            raise ValueError(f"Frozen source snapshot changed: {target}")
        if not target.exists():
            shutil.copyfile(path, target)
        saved[name] = {"snapshot": target.relative_to(out).as_posix(), "sha256": expected}
    license_path = upstream / "LICENSE"
    if license_path.exists():
        target = out / "source_snapshot/_upstream/LICENSE"
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            shutil.copyfile(license_path, target)
    write_json(out / "source_snapshot.json", {"files": saved})


def initialize_worker():
    # Each spawned process owns this cache; no inter-process mutable state.
    _WORKER_REGISTRY = {r["task"]: r for r in load_registry()}
    compute_cell.registry = _WORKER_REGISTRY


def compute_cell(job):
    phase, calibration_round, task, n, rel, repeats, master = job
    rec = compute_cell.registry[task]
    rows = []
    for replicate in range(repeats):
        seed = sample_seed(master, phase, calibration_round, task, n, rel, replicate)
        row = {"phase": phase, "calibration_round": calibration_round, "task": task,
               "klass": rec["klass"], "n": n, "rel": rel, "replicate": replicate,
               "base_seed": seed, "data_failure": None, "dataset_sha256": None,
               "r2": np.nan, "fit_failure": None, "outcome_failure": None,
               "legacy_endpoint_error": np.nan}
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
    return pd.DataFrame(rows)


def jobs_for(protocol, registry):
    jobs = []
    for phase in ("reference", "calibration", "test"):
        rounds = range(protocol["calibration_rounds"]) if phase != "test" else [-1]
        for calibration_round in rounds:
            for rec in registry:
                if phase != "test" and rec["klass"] != "H0":
                    continue
                for n in protocol["n_grid"]:
                    for rel in protocol["noise_grid"]:
                        jobs.append((phase, calibration_round, rec["task"], n, rel,
                                     protocol["replicates"][phase], protocol["master_seed"]))
    return jobs


def cell_name(job):
    phase, calibration_round, task, n, rel, _, _ = job
    return f"{phase}_r{calibration_round:02d}_{task}_n{n}_noise{round(rel*1e6)}.parquet"


def verify_cell(frame, job):
    phase, calibration_round, task, n, rel, count, master = job
    if (len(frame) != count or set(frame.replicate) != set(range(count))
            or not frame.phase.eq(phase).all() or not frame.task.eq(task).all()
            or not frame.n.eq(n).all() or not frame.rel.eq(rel).all()
            or not frame.calibration_round.eq(calibration_round).all()):
        raise ValueError(f"Cell key/count mismatch: {cell_name(job)}")
    expected = [sample_seed(master, phase, calibration_round, task, n, rel, int(rep)) for rep in frame.replicate]
    if frame.base_seed.tolist() != expected:
        raise ValueError(f"Cell seed mismatch: {cell_name(job)}")


def verify_new_records(frame, protocol):
    keys = ["phase", "calibration_round", "task", "n", "rel", "replicate"]
    if len(frame) != protocol["expected_datasets"] or frame.duplicated(keys).any():
        raise AssertionError("Unexpected v2 cardinality or duplicate record keys")
    if frame.base_seed.duplicated().any() or frame.dataset_sha256.dropna().duplicated().any():
        raise AssertionError("A sample was reused within v2")
    if not frame.loc[frame.phase != "test", "klass"].eq("H0").all():
        raise AssertionError("H1 in reference/calibration")
    if not frame.loc[frame.phase == "test", "calibration_round"].eq(-1).all():
        raise AssertionError("Test bank must not be attached to a fitted round")
    v1 = pd.read_parquet(V1 / "raw_policy.parquet", columns=["base_seed", "dataset_sha256"])
    old_seeds = set(v1.base_seed)
    old_seeds.update(pd.read_parquet(ROOT / INPUTS[0], columns=["base_seed"]).base_seed)
    old_seeds.update(pd.read_csv(ROOT / INPUTS[3], usecols=["seed"]).seed)
    if old_seeds & set(frame.base_seed):
        raise AssertionError("v2 seed overlaps v1 or a legacy run")
    if set(v1.dataset_sha256.dropna()) & set(frame.dataset_sha256.dropna()):
        raise AssertionError("v2 dataset overlaps v1")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", default=str(PACKAGE / "protocol_v2.json"))
    parser.add_argument("--out", default="runs/p01_v2")
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    protocol = json.loads(Path(args.protocol).read_text(encoding="utf-8"))
    if (protocol["protocol_id"] != "p01-independent-recalibration-v2"
            or protocol["permutations"] != 999 or protocol["cv_folds"] != 5
            or protocol["rbf"] != {"landmarks": 32, "bandwidth_squared": "number_of_features",
                                    "ridge_per_observation": 0.01}
            or protocol["methods"] != list(METHODS)):
        raise ValueError("Unsupported protocol; do not silently change a frozen design")
    out = (ROOT / args.out).resolve()
    if out.parent != (ROOT / "runs").resolve() or not out.name.startswith("p01_v2"):
        raise ValueError("v2 output must be a dedicated runs/p01_v2* directory")
    registry = load_registry()
    source = source_hashes()
    frozen = {"protocol": protocol, "source_sha256": source,
              "environment": {"python": platform.python_version(), "numpy": np.__version__,
                              "pandas": pd.__version__, "scipy": scipy.__version__,
                              "worker_blas_threads": 1},
              "registry": [{k: r[k] for k in ("task", "klass", "formula", "sampling")} for r in registry]}
    fingerprint = hashlib.sha256(json.dumps(frozen, sort_keys=True).encode()).hexdigest()
    jobs = jobs_for(protocol, registry)
    if sum(job[5] for job in jobs) != protocol["expected_datasets"]:
        raise ValueError("Frozen expected dataset count does not match planned jobs")
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = out / "manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest["fingerprint"] != fingerprint:
            raise ValueError("Source/protocol drift; frozen run may not be overwritten")
        if tree_hashes(V1) != manifest["v1_recursive_snapshot"]:
            raise ValueError("v1 changed since v2 freeze")
        if manifest["complete"]:
            if digest(out / "raw_policy.parquet") != manifest["raw_sha256"]:
                raise ValueError("Completed v2 raw hash mismatch")
            print("Existing complete frozen v2 verified; no generation or overwrite.")
            return
    else:
        if (out / "cells").exists():
            raise ValueError("Refuse orphan cells without a frozen manifest")
        manifest = {**frozen, "fingerprint": fingerprint, "complete": False,
                    "v1_recursive_snapshot": tree_hashes(V1),
                    "historical_inputs": {str(ROOT / name): {"sha256": digest(ROOT / name)} for name in INPUTS}}
        freeze_sources(out, source)
        write_json(manifest_path, manifest)
        write_json(out / "protocol_frozen.json", {"protocol": protocol,
                                                  "protocol_sha256": hashlib.sha256(json.dumps(protocol, sort_keys=True).encode()).hexdigest(),
                                                  "fingerprint": fingerprint,
                                                  "generation_started_after_freeze": True})
    cells = out / "cells"
    cells.mkdir(exist_ok=True)
    pending = []
    for job in jobs:
        path = cells / cell_name(job)
        if path.exists():
            verify_cell(pd.read_parquet(path), job)
        else:
            pending.append(job)
    # Spawned Windows workers inherit these before importing NumPy. No shared
    # mutable RNG, no nested BLAS oversubscription; job identity determines data.
    for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS"):
        os.environ[variable] = "1"
    started, done = time.perf_counter(), len(jobs) - len(pending)
    with ProcessPoolExecutor(max_workers=args.workers, initializer=initialize_worker) as pool:
        futures = {pool.submit(compute_cell, job): job for job in pending}
        for future in as_completed(futures):
            job = futures[future]
            frame = future.result()
            verify_cell(frame, job)
            path = cells / cell_name(job)
            temp = path.with_suffix(".parquet.tmp")
            frame.to_parquet(temp, index=False)
            temp.replace(path)
            done += 1
            if done % 10 == 0 or done == len(jobs):
                progress = {"completed_cells": done, "total_cells": len(jobs),
                            "last_cell": cell_name(job), "elapsed_seconds": time.perf_counter() - started}
                write_json(out / "progress.json", progress)
                print(f"[{done}/{len(jobs)}] {cell_name(job)}", flush=True)
    frame = pd.concat([pd.read_parquet(cells / cell_name(job)) for job in jobs], ignore_index=True)
    verify_new_records(frame, protocol)
    if tree_hashes(V1) != manifest["v1_recursive_snapshot"]:
        raise AssertionError("v1 preservation failed")
    target = out / "raw_policy.parquet"
    temp = target.with_suffix(".parquet.tmp")
    frame.to_parquet(temp, index=False)
    temp.replace(target)
    manifest.update({"complete": True, "rows": len(frame), "raw_sha256": digest(target),
                     "phase_rows": frame.phase.value_counts().to_dict(),
                     "cells": {cell_name(job): digest(cells / cell_name(job)) for job in jobs},
                     "v1_preserved": True, "v1_seed_overlap": 0, "v1_dataset_hash_overlap": 0,
                     "historical_seed_overlap": 0, "workers": args.workers,
                     "elapsed_seconds_this_invocation": time.perf_counter() - started})
    write_json(manifest_path, manifest)
    print(json.dumps({"complete": True, "rows": len(frame), "out": str(out)}), flush=True)


if __name__ == "__main__":
    main()
