"""Domain/truth preflight tests, independent of study samples and diagnostics.

Seeds below are unit-test/debug seeds only. These tests never fit a model,
compute a gate statistic, calibrate a threshold, or generate study outcomes.
"""
from __future__ import annotations

import itertools
import json
from collections import Counter

import numpy as np
import pandas as pd
import pytest

from experiments.fixed_domain.grids import build_grid, validate_grid
from experiments.fixed_domain.tasks import (
    _source_instances,
    build_registry,
    evaluate,
    sample,
)

# Fixed preflight-only seeds, outside the study seed-generation API.
PREFLIGHT_SEED = 89012345
REGISTRY = build_registry()
SOURCE_TASKS = [task for task in REGISTRY if task["family"] == "upstream_source"]
NEW_NULLS = [task for task in REGISTRY if task["family"] == "monomial_null"]
INTEGER_TASKS = [task for task in REGISTRY if "integer" in task["coordinate_scales"]]


def _corners(task):
    return np.asarray(list(itertools.product(*task["outer_bounds"])), dtype=float)


def _domain_anchor(task):
    anchor = []
    for (lo, hi), scale in zip(task["outer_bounds"], task["coordinate_scales"], strict=False):
        value = np.sqrt(lo * hi) if scale == "log" else (lo + hi) / 2
        anchor.append(float(round(value)) if scale == "integer" else float(value))
    return np.asarray(anchor)


def _high_precision_source_truth(task, points):
    """Evaluate the registered mathematical expression, without float exp−1."""
    expression = _source_instances()[task["params"]["source_id"]].sympy_eq
    symbols = sorted(expression.free_symbols, key=lambda value: int(str(value)[1:]))
    return np.asarray([
        float(expression.evalf(65, subs=dict(zip(symbols, point, strict=False))))
        for point in points
    ])


def test_registry_is_deterministic_json_with_registered_counts_and_roles():
    encoded = json.dumps(REGISTRY, sort_keys=True, allow_nan=False)
    assert json.loads(encoded) == REGISTRY
    assert encoded == json.dumps(build_registry(), sort_keys=True, allow_nan=False)
    assert len(REGISTRY) == len({task["task"] for task in REGISTRY}) == 110
    assert Counter(task["split"] for task in REGISTRY) == {
        "development": 36, "calibration": 12, "confirmation": 50,
        "reserve": 8, "stress": 4,
    }
    assert Counter(task["klass"] for task in REGISTRY) == {"H0": 37, "H1": 73}
    assert len(SOURCE_TASKS) == 38
    assert len(NEW_NULLS) == 24
    assert all(task["klass"] == "H0" for task in REGISTRY if task["split"] == "calibration")
    assert all(task["split"] == "reserve" for task in REGISTRY
               if task["family"] in {"asinh", "lambertw"})
    assert all(task["split"] == "stress" for task in REGISTRY
               if task["family"] == "threshold")


def test_original_source_identities_and_expressions_are_preserved():
    for task in SOURCE_TASKS:
        source_id = task["params"]["source_id"]
        original = _source_instances()[source_id]
        assert task["source_equations"] == [source_id]
        assert task["task"] == source_id
        assert task["expression"] == str(original.sympy_eq)
        assert task["n_features"] == original.get_var_count()
        if task["source_kind"] == "reused_source_restricted_domain":
            assert task["split"] == "development"


@pytest.mark.parametrize("task", REGISTRY, ids=lambda task: task["task"])
def test_complete_outer_corners_and_fixed_grid_have_positive_truth(task):
    outer = np.asarray(task["outer_bounds"], dtype=float)
    train = np.asarray(task["train_bounds"], dtype=float)
    assert outer.shape == train.shape == (task["n_features"], 2)
    assert (outer[:, 0] > 0).all()
    assert (outer[:, 0] <= train[:, 0]).all()
    assert (train[:, 0] < train[:, 1]).all()
    assert (train[:, 1] <= outer[:, 1]).all()
    assert isinstance(task["validity_argument"], str) and task["validity_argument"].strip()
    corner_truth = evaluate(task, _corners(task))
    assert np.isfinite(corner_truth).all() and (corner_truth > 0).all()
    grid = build_grid(task)
    facts = validate_grid(task, grid, lambda points: evaluate(task, points))
    assert facts["validated"] and facts["inside_outer_domain"] and facts["outside_train"]
    assert facts["finite_positive_truth"] and facts["analyticproofpresent"]
    assert facts["n_points"] == len(grid["X"])
    # Rebuilding geometry cannot depend on samples, truth magnitude or noise.
    again = build_grid(task)
    np.testing.assert_array_equal(grid["X"], again["X"])
    np.testing.assert_array_equal(grid["weights"], again["weights"])


@pytest.mark.parametrize("task", SOURCE_TASKS, ids=lambda task: task["task"])
def test_source_truth_matches_high_precision_mathematical_expression(task):
    corners = _corners(task)
    points = np.vstack([corners[0], _domain_anchor(task), corners[-1]])
    actual = evaluate(task, points)
    expected = _high_precision_source_truth(task, points)
    np.testing.assert_allclose(actual, expected, rtol=3e-13, atol=0)


def test_i41_stable_expm1_corrects_cancellation_without_changing_expression():
    task = next(task for task in SOURCE_TASKS if task["task"] == "feynman-i.41.16")
    # These coordinates lie in the fixed outer box; Planck arguments are tiny.
    points = np.asarray([[.13, 790.], [7.9, 31.7], [.63, 158.5]])
    bounds = np.asarray(task["outer_bounds"])
    assert ((points >= bounds[:, 0]) & (points <= bounds[:, 1])).all()
    expected = _high_precision_source_truth(task, points)
    stable = evaluate(task, points)
    np.testing.assert_allclose(stable, expected, rtol=3e-14, atol=0)
    source = _source_instances()[task["params"]["source_id"]]
    unstable = np.asarray([source.eq_func(point) for point in points])
    # Deliberately do not assert legacy numeric parity for cancellation-prone
    # exp(argument)-1. The expression itself remains exactly the source law.
    assert np.max(np.abs(unstable - expected) / expected) > .001
    assert task["expression"] == str(source.sympy_eq)


@pytest.mark.parametrize("task", NEW_NULLS, ids=lambda task: task["task"])
def test_new_null_truth_is_exact_registered_log_affine_law(task):
    points = np.vstack([_corners(task), _domain_anchor(task)])
    parameters = task["params"]
    exponents = np.asarray(parameters["exponents"])
    np.testing.assert_array_equal(
        exponents,
        np.asarray(parameters["exponent_numerators"]) / parameters["exponent_denominator"],
    )
    expected_log = np.log(parameters["coefficient"]) + np.log(points) @ exponents
    np.testing.assert_allclose(np.log(evaluate(task, points)), expected_log, rtol=0, atol=2e-15)


def test_new_null_design_groups_and_exponents_are_distinct_across_roles():
    calibration = [task for task in NEW_NULLS if task["split"] == "calibration"]
    confirmation = [task for task in NEW_NULLS if task["split"] == "confirmation"]
    assert Counter(task["n_features"] for task in calibration) == {1: 3, 2: 3, 3: 3, 4: 3}
    assert Counter(task["n_features"] for task in confirmation) == {1: 3, 2: 3, 3: 3, 4: 3}
    assert {task["design_group"] for task in calibration}.isdisjoint(
        task["design_group"] for task in confirmation)
    assert {tuple(task["params"]["exponents"]) for task in calibration}.isdisjoint(
        tuple(task["params"]["exponents"]) for task in confirmation)
    assert all(task["structural_group"] == "H0_null_family" for task in NEW_NULLS)
    assert {task["params"]["design_distribution"] for task in calibration} == {"log_uniform"}
    assert {task["params"]["design_distribution"] for task in confirmation} == {"log_beta"}


@pytest.mark.parametrize("task", INTEGER_TASKS, ids=lambda task: task["task"])
def test_integer_sampling_uses_inclusive_support_and_grid_stays_integral(task):
    data, features = sample(task, PREFLIGHT_SEED, 4096, .03)
    grid = build_grid(task)
    assert "inclusive" in task["integer_support"]
    for j, scale in enumerate(task["coordinate_scales"]):
        if scale != "integer":
            continue
        lo, hi = task["train_bounds"][j]
        values = data[features[j]].to_numpy()
        assert lo == int(lo) and hi == int(hi)
        np.testing.assert_array_equal(values, np.rint(values))
        assert values.min() == lo and values.max() == hi
        np.testing.assert_array_equal(grid["X"][:, j], np.rint(grid["X"][:, j]))


@pytest.mark.parametrize("task", REGISTRY, ids=lambda task: task["task"])
def test_debug_sampling_is_deterministic_positive_and_inside_train_box(task):
    data, features = sample(task, PREFLIGHT_SEED, 13, .03)
    repeated, repeated_features = sample(task, PREFLIGHT_SEED, 13, .03)
    pd.testing.assert_frame_equal(data, repeated, check_exact=True)
    assert features == repeated_features == [f"x{j + 1}" for j in range(task["n_features"])]
    values = data[features].to_numpy()
    bounds = np.asarray(task["train_bounds"])
    assert ((values >= bounds[:, 0]) & (values <= bounds[:, 1])).all()
    assert np.isfinite(data.to_numpy()).all() and (data.to_numpy() > 0).all()
    different, _ = sample(task, PREFLIGHT_SEED + 1, 13, .03)
    assert not data.equals(different)
    noiseless, _ = sample(task, PREFLIGHT_SEED, 13, 0.)
    np.testing.assert_array_equal(noiseless[features].to_numpy(), values)
    np.testing.assert_allclose(
        noiseless.target.to_numpy(), evaluate(task, values), rtol=3e-14, atol=0,
    )


def test_source_domain_constraints_are_strict_over_complete_boxes():
    by_name = {task["task"]: task for task in SOURCE_TASKS}
    for task in SOURCE_TASKS:
        corners = _corners(task)
        if task["structural_group"] == "lorentz_radical":
            assert np.max(corners[:, 1]) < 299792458
        if task["task"] in {"feynman-i.39.11", "feynman-i.43.43"}:
            assert np.min(corners[:, 0]) >= 1.05
        if task["task"] in {"feynman-ii.11.27", "feynman-ii.11.28"}:
            assert np.max(corners[:, 0] * corners[:, 1]) < 3
    assert by_name["feynman-i.37.4"]["outer_bounds"][2] == [.1, 1.4]
    assert by_name["feynman-i.26.2"]["outer_bounds"][0][0] > 0
    log_box = np.asarray(by_name["feynman-i.44.4"]["outer_bounds"])
    assert log_box[2, 0] / log_box[3, 1] >= 3 - 1e-14


def test_evaluate_fails_on_nonpositive_or_malformed_coordinates():
    task = REGISTRY[0]
    with pytest.raises(ValueError, match="n-by"):
        evaluate(task, np.ones(task["n_features"]))
    with pytest.raises(ValueError, match="n-by"):
        evaluate(task, np.ones((3, task["n_features"] + 1)))
    points = np.tile(_domain_anchor(task), (3, 1))
    points[1, 0] = 0
    with pytest.raises(ValueError, match="strictly positive"):
        evaluate(task, points)
