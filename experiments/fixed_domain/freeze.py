"""Freeze domain geometry, structural allocation and all executable sources."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

from experiments.fixed_domain.grids import build_grid, validate_grid
from experiments.fixed_domain.grouping import assert_group_isolation
from experiments.fixed_domain.tasks import build_registry, evaluate
from experiments.release_validation.inventory import ROOT, digest
from experiments.release_validation.runner import write_json
from experiments.release_validation.runner_v2 import (
    freeze_sources,
    tree_hashes,
)
from experiments.release_validation.runner_v2 import (
    source_hashes as prior_sources,
)

PACKAGE = Path(__file__).resolve().parent


def verify_frozen(run):
    run = Path(run)
    manifest = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    for name, expected in manifest["source_hashes"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"Frozen source drift: {name}")
    for relative, hashes in manifest["preservation"].items():
        if tree_hashes(ROOT / relative) != hashes:
            raise ValueError(f"Historical results changed: {relative}")
    for name, key in (("task_registry.json", "task_registry_sha256"),
                      ("grid_manifest.json", "grid_manifest_sha256"),
                      ("protocol_frozen.json", "protocol_sha256")):
        if digest(run / name) != manifest[key]:
            raise ValueError(f"Frozen artifact drift: {name}")
    for task, record in json.loads((run / "grid_manifest.json").read_text(encoding="utf-8")).items():
        if digest(run / record["file"]) != record["file_sha256"]:
            raise ValueError(f"Frozen grid drift: {task}")
    return manifest


def freeze(run):
    run = Path(run).resolve()
    if run.parent != (ROOT / "runs").resolve() or not run.name.startswith("p02_fixed_domain"):
        raise ValueError("Use a dedicated runs/p02_fixed_domain* output")
    if (run / "manifest.json").exists():
        return verify_frozen(run)
    if run.exists() and any(run.iterdir()):
        raise ValueError("Refuse an existing nonempty directory without a frozen manifest")
    protocol = json.loads((PACKAGE / "protocol.json").read_text(encoding="utf-8"))
    registry = build_registry()
    groups = assert_group_isolation(registry)
    sources = prior_sources()
    sources.update({str(path): digest(path) for path in PACKAGE.glob("*.py")})
    sources[str(PACKAGE / "protocol.json")] = digest(PACKAGE / "protocol.json")
    for name in ("runs/p02_inventory/inventory.json", "runs/p02_inventory/source_domain_review.json"):
        sources[str(ROOT / name)] = digest(ROOT / name)
    for required in ("runner.py", "analyze.py"):
        if str(PACKAGE / required) not in sources:
            raise ValueError(f"Missing analysis/execution implementation before freeze: {required}")
    preservation = {name: tree_hashes(ROOT / name) for name in protocol["preservation"]}
    run.mkdir(parents=True, exist_ok=True)
    (run / "grids").mkdir()
    write_json(run / "protocol_frozen.json", protocol)
    write_json(run / "task_registry.json", registry)
    write_json(run / "group_isolation.json", groups)
    grids = {}
    for task in registry:
        grid = build_grid(task, seed=protocol["grid"]["seed"])
        truth = evaluate(task, grid["X"])
        valid = validate_grid(task, grid, truth)
        relative = f"grids/{task['task']}.npz"
        np.savez_compressed(run / relative, X=grid["X"], weights=grid["weights"], truth=truth)
        grids[task["task"]] = {**valid, "file": relative, "file_sha256": digest(run / relative),
                               "directions": grid["directions"], "metadata": grid["metadata"]}
    write_json(run / "grid_manifest.json", grids)
    freeze_sources(run, sources)
    manifest = {"protocol": protocol, "source_hashes": sources, "preservation": preservation,
                "task_registry_sha256": digest(run / "task_registry.json"),
                "grid_manifest_sha256": digest(run / "grid_manifest.json"),
                "protocol_sha256": digest(run / "protocol_frozen.json"),
                "frozen_at_utc": datetime.now(UTC).isoformat(),
                "freeze_before_study_sampling": True,
                "preflight_scope": "Truth/domain and deterministic software regression checks with debug seeds; no study draws, no outcome-driven selection",
                "environment": {"python": platform.python_version(), "numpy": np.__version__,
                                "pandas": pd.__version__, "scipy": scipy.__version__, "worker_blas_threads": 1}}
    manifest["fingerprint"] = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest()
    write_json(run / "manifest.json", manifest)
    verify_frozen(run)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", default="runs/p02_fixed_domain")
    args = parser.parse_args()
    manifest = freeze(ROOT / args.run)
    print(f"Frozen {manifest['fingerprint']}; prior v1/v2 trees verified.")


if __name__ == "__main__":
    main()
