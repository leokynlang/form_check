"""Meaningful invariants for fixed geometry and the registered region outcome."""
import copy

import numpy as np
import pytest

from experiments.fixed_domain.grids import build_grid, error_metrics, validate_grid


def example_task(scales=("log", "linear")):
    dimensions = len(scales)
    return {
        "task": "grid-test-only",
        "n_features": dimensions,
        "outer_bounds": [[1., 16.] for _ in scales],
        "train_bounds": [[3., 9.] for _ in scales],
        "coordinate_scales": list(scales),
        "validity_argument": "The test truth is a positive monomial on a positive compact box.",
    }


def test_grid_is_external_inside_outer_positive_and_component_balanced():
    task = example_task()
    before = copy.deepcopy(task)
    grid = build_grid(task)
    points = grid["X"]
    report = validate_grid(task, grid, lambda x: np.prod(x, axis=1))
    assert task == before
    assert report["validated"] and report["outside_train"]
    assert np.all(points >= 1) and np.all(points <= 16)
    assert np.any((points < 3) | (points > 9), axis=1).all()
    assert np.isclose(grid["weights"].sum(), 1)
    assert len(np.unique(points, axis=0)) == len(points)
    components = grid["metadata"]["axial_components"]
    assert len(components) == 4
    assert all(component["mass"] == .125 for component in components.values())
    assert all(component["unique_points"] == 8 for component in components.values())
    assert grid["metadata"]["joint_total_mass"] == .5
    np.testing.assert_allclose(grid["metadata"]["anchor"], [np.sqrt(27), 6])
    assert any("axis_x1_lower" in labels for labels in grid["directions"])


def test_geometry_and_hashes_are_deterministic_without_mutating_sampling_rng():
    task = example_task()
    np.random.seed(7821)
    state = np.random.get_state()
    first = build_grid(task)
    second = build_grid(task)
    after = np.random.get_state()
    assert state[0] == after[0] and state[2:] == after[2:]
    np.testing.assert_array_equal(state[1], after[1])
    np.testing.assert_array_equal(first["X"], second["X"])
    np.testing.assert_array_equal(first["weights"], second["weights"])
    one = validate_grid(task, first, lambda x: x[:, 0] ** .7)
    two = validate_grid(task, second, lambda x: x[:, 0] ** .7)
    assert one == two
    alternative = build_grid(task, seed=123)
    assert one["coordinates_sha256"] != validate_grid(task, alternative, lambda x: x[:, 0] ** .7)["coordinates_sha256"]


def test_integer_support_deduplicates_and_merges_component_mass():
    task = example_task(("integer",))
    task["outer_bounds"] = [[1, 10]]
    task["train_bounds"] = [[3, 7]]
    grid = build_grid(task)
    np.testing.assert_array_equal(grid["X"][:, 0], [1, 2, 8, 9, 10])
    assert np.equal(grid["X"], np.rint(grid["X"])).all()
    assert len(grid["metadata"]["axial_components"]) == 2
    assert grid["metadata"]["axial_components"]["axis_x1_lower"]["unique_points"] == 2
    # Each direction receives .25 and the joint unique support receives .5.
    np.testing.assert_allclose(grid["weights"], [.225, .225, 1/12 + .1, 1/12 + .1, 1/12 + .1])
    assert all("joint" in labels for labels in grid["directions"])
    validate_grid(task, grid, lambda x: x[:, 0])


def test_grid_rejects_truth_failure_without_changing_geometry():
    task = example_task()
    grid = build_grid(task)
    points = grid["X"].copy()
    weights = grid["weights"].copy()
    truth = np.ones(len(points))
    truth[5] = 0
    with pytest.raises(ValueError, match="truth must be finite and strictly positive"):
        validate_grid(task, grid, truth)
    np.testing.assert_array_equal(points, grid["X"])
    np.testing.assert_array_equal(weights, grid["weights"])


def test_validation_requires_analytic_review_argument_and_entire_grid_validity():
    task = example_task()
    grid = build_grid(task)
    task["validity_argument"] = ""
    with pytest.raises(ValueError, match="analytic validity argument"):
        validate_grid(task, grid, np.ones(len(grid["X"])))
    task["validity_argument"] = "positive test monomial"
    grid["X"][0] = [6, 6]
    with pytest.raises(ValueError, match="inside or on the training box"):
        validate_grid(task, grid, np.ones(len(grid["X"])))


def test_weighted_region_metrics_use_weighted_not_unweighted_p95():
    # A large-error point has only 1% mass, so the weighted p95 remains .2.
    truth = np.array([1., 2., 1e-100])
    errors = np.array([.1, .2, 10.])
    predictions = np.log(truth) + np.log1p(errors)
    metrics = error_metrics(predictions, truth, [.5, .49, .01])
    assert metrics["weighted_mean_relative_error"] == pytest.approx(.248)
    assert metrics["weighted_p95_relative_error"] == pytest.approx(.2)
    assert metrics["max_relative_error"] == pytest.approx(10)
    perfect = error_metrics(np.log(truth), truth, [.5, .49, .01])
    assert perfect == {"weighted_mean_relative_error": 0, "weighted_p95_relative_error": 0, "max_relative_error": 0}


@pytest.mark.parametrize("predictions, truth, weights, message", [
    ([0., np.nan], [1., 1.], [.5, .5], "log prediction"),
    ([0., 1e4], [1., 1.], [.5, .5], "nonfinite_relative_error"),
    ([0., 0.], [1., 0.], [.5, .5], "truth value"),
    ([0., 0.], [1., 1.], [.5, 0.], "strictly positive"),
    ([0., 0.], [1., 1.], [.2, .2], "sum to one"),
])
def test_metrics_never_drop_invalid_points(predictions, truth, weights, message):
    with pytest.raises(ValueError, match=message):
        error_metrics(predictions, truth, weights)


def test_non_nested_domains_and_empty_discrete_external_support_fail():
    task = example_task(("integer",))
    task["train_bounds"] = [[.5, 17]]
    with pytest.raises(ValueError, match="nested box"):
        build_grid(task)
    task["outer_bounds"] = [[1, 2]]
    task["train_bounds"] = [[1, 2]]
    with pytest.raises(ValueError, match="No representable axial"):
        build_grid(task)
