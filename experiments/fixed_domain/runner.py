"""Fresh P0-2 draws and resumable cells under a pre-frozen physical domain."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

from experiments.fixed_domain.freeze import verify_frozen
from experiments.fixed_domain.grids import error_metrics
from experiments.fixed_domain.grouping import group_unit
from experiments.fixed_domain.tasks import evaluate, sample
from experiments.release_validation.gates import METHODS, run_gates
from experiments.release_validation.inventory import INPUTS, ROOT, digest
from experiments.release_validation.runner import write_json

PHASES = {"reference": 21, "calibration": 22, "test": 23}


def sample_seed(master, phase, calibration_round, task, n, rel, replicate):
    task_id = np.frombuffer(hashlib.sha256(task.encode()).digest()[:8], dtype="<u4")
    entropy = [master, 3, PHASES[phase], calibration_round + 1, *(int(v) for v in task_id),
               n, round(rel * 1_000_000), replicate]
    return int(np.random.SeedSequence(entropy).generate_state(1, dtype=np.uint64)[0]) % (2**62) + 2**40


def jobs_for(protocol, registry, stage):
    spec = protocol["stages"][stage]
    jobs = []
    for phase in PHASES:
        for rnd in (range(protocol["calibration_rounds"]) if phase != "test" else [-1]):
            for task in registry:
                if phase != "test" and task["split"] != "calibration":
                    continue
                if phase == "test" and task["split"] not in {"development", "confirmation", "stress"}:
                    continue
                for n in spec["n_grid"]:
                    for rel in spec["noise_grid"]:
                        if stage == "stage2" and n in protocol["stages"]["stage1"]["n_grid"] and rel in protocol["stages"]["stage1"]["noise_grid"]:
                            continue
                        repeats = protocol[f"{phase}_replicates"] if phase != "calibration" else protocol["calibration_replicates"]
                        jobs.append((phase, rnd, task["task"], n, rel, repeats, protocol["master_seed"]))
    return jobs


def cell_name(job):
    phase, rnd, task, n, rel, _, _ = job
    return f"{phase}_r{rnd:02d}_{task}_n{n}_noise{round(rel*1e6)}.parquet"


def initialize_worker(run):
    run = Path(run)
    compute_cell.registry = {t["task"]: t for t in json.loads((run / "task_registry.json").read_text(encoding="utf-8"))}
    manifest = json.loads((run / "grid_manifest.json").read_text(encoding="utf-8"))
    compute_cell.grids = {}
    for name, record in manifest.items():
        with np.load(run / record["file"], allow_pickle=False) as grid:
            compute_cell.grids[name] = {k: grid[k] for k in ("X", "weights", "truth")}
        compute_cell.grids[name]["sha256"] = record["coordinates_sha256"]


def fit_outcomes(data, features, task, grid):
    x = np.log(data[features].to_numpy())
    y = np.log(data.target.to_numpy())
    design = np.column_stack([np.ones(len(x)), x])
    coef, _, rank, _ = np.linalg.lstsq(design, y, rcond=None)
    result = {"r2": np.nan, "coef_json": json.dumps(coef.tolist()), "fit_failure": None,
              "outcome_failure": None, "legacy_endpoint_error": np.nan,
              "legacy_endpoint_failure": None}
    if rank != design.shape[1]:
        result["fit_failure"] = "rank_deficient"
        return result
    tss = float(np.sum((y - y.mean()) ** 2))
    if not np.isfinite(tss) or tss <= 0:
        result["fit_failure"] = "degenerate_response"
        return result
    result["r2"] = 1 - float(np.sum((y - design @ coef) ** 2)) / tss
    pred_log = np.column_stack([np.ones(len(grid["X"])), np.log(grid["X"])]) @ coef
    result.update(error_metrics(pred_log, grid["truth"], grid["weights"]))
    point = np.median(x, axis=0)
    point[0] = x[:, 0].max() + np.log(3.0)
    external = np.exp(point)
    bounds = np.asarray(task["outer_bounds"])
    if (external < bounds[:, 0]).any() or (external > bounds[:, 1]).any():
        result["legacy_endpoint_failure"] = "outside_registered_physical_domain"
    elif any(scale == "integer" and not float(external[j]).is_integer()
             for j, scale in enumerate(task["coordinate_scales"])):
        result["legacy_endpoint_failure"] = "outside_registered_integer_support"
    else:
        truth = evaluate(task, external[None])[0]
        with np.errstate(over="ignore"):
            error = float(abs(np.expm1(np.r_[1., point] @ coef - np.log(truth))))
        if np.isfinite(error):
            result["legacy_endpoint_error"] = error
        else:
            result["legacy_endpoint_failure"] = "nonfinite_prediction"
    return result


def compute_cell(job):
    phase, rnd, name, n, rel, repeats, master = job
    task, grid = compute_cell.registry[name], compute_cell.grids[name]
    rows = []
    for replicate in range(repeats):
        seed = sample_seed(master, phase, rnd, name, n, rel, replicate)
        row = {k: task.get(k) for k in ("task", "klass", "split", "structural_group", "design_group", "source_kind", "regime_kind")}
        row.update({"phase": phase, "calibration_round": rnd, "n": n, "rel": rel,
                    "replicate": replicate, "base_seed": seed, "group_unit": group_unit(task),
                    "data_failure": None, "fit_failure": None, "outcome_failure": None,
                    "dataset_sha256": None, "grid_sha256": grid["sha256"], "r2": np.nan,
                    "coef_json": None, "weighted_mean_relative_error": np.nan,
                    "weighted_p95_relative_error": np.nan, "max_relative_error": np.nan,
                    "legacy_endpoint_error": np.nan, "legacy_endpoint_failure": None})
        for method in METHODS:
            row.update({f"score_{method}": np.nan, f"fail_{method}": "data_unavailable"})
        try:
            data, features = sample(task, seed, n, rel)
        except (ValueError, OverflowError, FloatingPointError) as exc:
            row["data_failure"] = f"generation:{type(exc).__name__}:{exc}"
            rows.append(row)
            continue
        row.update(run_gates(data, features, seed))
        # A nonfinite fixed-grid prediction is an explicit outcome failure;
        # it can never be removed and reaveraged over the remaining points.
        try:
            row.update(fit_outcomes(data, features, task, grid))
        except (ValueError, OverflowError, FloatingPointError) as exc:
            row["outcome_failure"] = f"fixed_grid:{type(exc).__name__}:{exc}"
        rows.append(row)
    return pd.DataFrame(rows)


def verify_cell(frame, job):
    phase, rnd, name, n, rel, count, master = job
    if (len(frame) != count or set(frame.replicate) != set(range(count))
            or not frame.phase.eq(phase).all() or not frame.task.eq(name).all()
            or not frame.n.eq(n).all() or not frame.rel.eq(rel).all()
            or not frame.calibration_round.eq(rnd).all()):
        raise ValueError(f"Cell key/count mismatch: {cell_name(job)}")
    expected = [sample_seed(master, phase, rnd, name, n, rel, int(rep)) for rep in frame.replicate]
    if frame.base_seed.tolist() != expected:
        raise ValueError(f"Cell seed mismatch: {cell_name(job)}")


def verify_records(raw, jobs):
    keys = ["phase", "calibration_round", "task", "n", "rel", "replicate"]
    if len(raw) != sum(j[5] for j in jobs) or raw.duplicated(keys).any():
        raise AssertionError("Unexpected dataset cardinality or duplicate record")
    if raw.base_seed.duplicated().any() or raw.dataset_sha256.dropna().duplicated().any():
        raise AssertionError("Data reused within study")
    if not raw.loc[raw.phase.ne("test"), "split"].eq("calibration").all():
        raise AssertionError("Reference/calibration role leakage")
    for prior in ("p01_pilot_v1", "p01_v2"):
        old = pd.read_parquet(ROOT / f"runs/{prior}/raw_policy.parquet", columns=["base_seed", "dataset_sha256"])
        if set(old.base_seed) & set(raw.base_seed) or set(old.dataset_sha256.dropna()) & set(raw.dataset_sha256.dropna()):
            raise AssertionError(f"Dataset/seed overlap with {prior}")
    old_seeds = set(pd.read_parquet(ROOT / INPUTS[0], columns=["base_seed"]).base_seed)
    old_seeds.update(pd.read_csv(ROOT / INPUTS[3], usecols=["seed"]).seed)
    if old_seeds & set(raw.base_seed):
        raise AssertionError("Legacy seed overlap")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", default="runs/p02_fixed_domain")
    parser.add_argument("--stage", choices=["stage1", "stage2"], default="stage1")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    run = (ROOT / args.run).resolve()
    frozen = verify_frozen(run)
    protocol = frozen["protocol"]
    registry = json.loads((run / "task_registry.json").read_text(encoding="utf-8"))
    jobs = jobs_for(protocol, registry, args.stage)
    out = run / args.stage
    out.mkdir(exist_ok=True)
    manifest_path = out / "manifest.json"
    if manifest_path.exists():
        existing = json.loads(manifest_path.read_text(encoding="utf-8"))
        if existing["fingerprint"] != frozen["fingerprint"]:
            raise ValueError("Stage belongs to another protocol")
        if existing["complete"]:
            if digest(out / "raw.parquet") != existing["raw_sha256"]:
                raise ValueError("Completed raw data changed")
            print("Existing completed stage verified; no redraw or overwrite.")
            return
    cells = out / "cells"
    cells.mkdir(exist_ok=True)
    pending = []
    for job in jobs:
        path = cells / cell_name(job)
        if path.exists():
            verify_cell(pd.read_parquet(path), job)
        else:
            pending.append(job)
    started = time.perf_counter()
    stage_manifest = {"fingerprint": frozen["fingerprint"], "stage": args.stage, "complete": False,
                      "expected_cells": len(jobs), "expected_datasets": sum(j[5] for j in jobs)}
    write_json(manifest_path, stage_manifest)
    for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS"):
        os.environ[variable] = "1"
    done = len(jobs) - len(pending)
    with ProcessPoolExecutor(max_workers=args.workers, initializer=initialize_worker, initargs=(str(run),)) as pool:
        futures = {pool.submit(compute_cell, job): job for job in pending}
        for future in as_completed(futures):
            job = futures[future]
            frame = future.result()
            verify_cell(frame, job)
            path = cells / cell_name(job)
            temp = path.with_suffix(".tmp.parquet")
            frame.to_parquet(temp, index=False)
            temp.replace(path)
            done += 1
            if done % 10 == 0 or done == len(jobs):
                progress = {"completed_cells": done, "total_cells": len(jobs),
                            "elapsed_seconds": round(time.perf_counter() - started, 2)}
                write_json(out / "progress.json", progress)
                print(json.dumps(progress), flush=True)
    raw = pd.concat([pd.read_parquet(cells / cell_name(job)) for job in jobs], ignore_index=True)
    verify_records(raw, jobs)
    verify_frozen(run)
    temp = out / "raw.tmp.parquet"
    raw.to_parquet(temp, index=False)
    temp.replace(out / "raw.parquet")
    stage_manifest.update({"complete": True, "raw_sha256": digest(out / "raw.parquet"),
                           "raw_rows": len(raw), "elapsed_seconds": time.perf_counter() - started,
                           "data_failures": int(raw.data_failure.notna().sum()),
                           "fit_failures": int(raw.fit_failure.notna().sum()),
                           "outcome_failures": int(raw.outcome_failure.notna().sum())})
    write_json(manifest_path, stage_manifest)
    print(json.dumps(stage_manifest), flush=True)


if __name__ == "__main__":
    main()
