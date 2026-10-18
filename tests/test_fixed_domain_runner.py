"""Execution identity and physically restricted endpoint regression tests."""
import json

import numpy as np

from experiments.fixed_domain.grids import build_grid
from experiments.fixed_domain.runner import fit_outcomes, jobs_for, sample_seed
from experiments.fixed_domain.tasks import build_registry, evaluate, sample
from experiments.release_validation.inventory import ROOT


def test_registered_stage_cardinality_and_no_reserve_draws():
    protocol = json.loads((ROOT / "experiments/fixed_domain/protocol.json").read_text())
    registry = build_registry()
    jobs = jobs_for(protocol, registry, "stage1")
    assert len(jobs) == 660
    assert sum(j[5] for j in jobs) == 41040
    splits = {r["task"]: r["split"] for r in registry}
    assert all(splits[j[2]] != "reserve" for j in jobs)
    assert all(splits[j[2]] == "calibration" for j in jobs if j[0] != "test")
    stage2 = jobs_for(protocol, registry, "stage2")
    assert not {j[:5] for j in stage2} & {j[:5] for j in jobs}


def test_seed_namespace_is_disjoint_by_role_round_and_task():
    seeds = {sample_seed(2026101803, phase, rnd, task, 50, .03, rep)
             for phase in ("reference", "calibration", "test")
             for rnd in (0, 1) for task in ("a", "b") for rep in range(3)}
    assert len(seeds) == 36
    assert min(seeds) > 2**40
    assert max(seeds) < 2**63


def test_fixed_grid_exact_monomial_and_invalid_legacy_endpoint():
    task = next(t for t in build_registry() if t["split"] == "calibration" and t["n_features"] == 2)
    grid = build_grid(task)
    grid["truth"] = evaluate(task, grid["X"])
    data, features = sample(task, 202610180312345, 200, 0)
    result = fit_outcomes(data, features, task, grid)
    assert result["r2"] > .99999999
    assert result["max_relative_error"] < 1e-12
    assert np.isfinite(result["weighted_mean_relative_error"])


def test_legacy_endpoint_never_calls_invalid_physical_truth():
    task = next(t for t in build_registry() if t["family"] == "ellipk")
    grid = build_grid(task)
    grid["truth"] = evaluate(task, grid["X"])
    data, features = sample(task, 202610180399991, 200, .03)
    # The first-coordinate 3x endpoint exceeds the predeclared upper m=.7.
    result = fit_outcomes(data, features, task, grid)
    assert result["legacy_endpoint_failure"] == "outside_registered_physical_domain"
    assert np.isnan(result["legacy_endpoint_error"])
    assert np.isfinite(result["weighted_mean_relative_error"])
