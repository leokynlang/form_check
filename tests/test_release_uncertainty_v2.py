"""Regression checks for the frozen repeated-calibration uncertainty estimand."""
from __future__ import annotations

import numpy as np
import pytest

from experiments.release_validation.analyze import RATIOS
from experiments.release_validation.analyze_v2 import (
    aggregate_policy,
    bootstrap_design,
    mean_uncertainty,
    paired_aggregate,
    ratio_arrays,
)


def _event_cube():
    """Three calibration rounds on two H0 and two H1 test tasks."""
    columns = sorted({column for pair in RATIOS.values() for column in pair})
    cube = np.zeros((3, 4, len(columns)), dtype=float)
    for round_id in range(3):
        for task in range(4):
            is_h0 = task < 2
            dangerous = task == 2
            withheld = 1 + round_id if is_h0 else (8 - round_id if dangerous else 2)
            released = 10 - withheld
            counts = {
                "total": 10, "released": released,
                "h0_total": 10 if is_h0 else 0,
                "h0_screened": 10 if is_h0 else 0,
                "h0_released": released if is_h0 else 0,
                "h0_withheld_screened": withheld if is_h0 else 0,
                "h0_alarm_screened": withheld if is_h0 else 0,
                "h0_decline_screened": 0,
                "h1_total": 0 if is_h0 else 10,
                "h1_alarm": 0 if is_h0 else withheld,
                "h1_decline": 0,
                "h1_danger_total": 10 if dangerous else 0,
                "h1_danger_screened": 10 if dangerous else 0,
                "h1_danger_released": released if dangerous else 0,
                "h1_safe_screened": 10 if task == 3 else 0,
                "h1_safe_released": released if task == 3 else 0,
                "all_danger_released": released if dangerous else 0,
                "released_known_outcome": released,
                "outcome_unavailable": 0,
            }
            for column, value in counts.items():
                cube[round_id, task, columns.index(column)] = value
    return cube, columns


def test_bootstrap_preserves_13_h0_and_19_h1_tasks_in_every_draw():
    # Interleave classes to catch accidental assumptions that strata are slices.
    classes = np.array(["H1", "H0"] * 13 + ["H1"] * 6)
    weights, rounds = bootstrap_design(classes, 20, 400, 20261018)
    assert weights.shape == (400, 32) and rounds.shape == (400, 20)
    assert np.issubdtype(weights.dtype, np.integer) and np.all(weights >= 0)
    assert np.array_equal(weights[:, classes == "H0"].sum(axis=1), np.full(400, 13))
    assert np.array_equal(weights[:, classes == "H1"].sum(axis=1), np.full(400, 19))
    assert np.all((rounds >= 0) & (rounds < 20))
    again_weights, again_rounds = bootstrap_design(classes, 20, 400, 20261018)
    assert np.array_equal(weights, again_weights) and np.array_equal(rounds, again_rounds)
    new_weights, new_rounds = bootstrap_design(classes, 20, 400, 20261019)
    assert not np.array_equal(weights, new_weights)
    assert not np.array_equal(rounds, new_rounds)


def test_one_task_weight_vector_is_shared_across_all_rounds_and_policies():
    cube = np.array([[[0., 10.], [10., 10.]],
                     [[10., 10.], [0., 10.]]])
    weights = np.array([[2, 0], [0, 2], [1, 1]])
    original, weighted = ratio_arrays(cube, ["num", "den"], "num", "den", weights)
    assert np.array_equal(original, [.5, .5])
    assert np.array_equal(weighted, [[0., 1.], [1., 0.], [.5, .5]])
    # A second policy exactly complements the first on every task and round.
    complement = cube.copy()
    complement[:, :, 0] = complement[:, :, 1] - complement[:, :, 0]
    _, complement_weighted = ratio_arrays(complement, ["num", "den"], "num", "den", weights)
    assert np.array_equal(weighted + complement_weighted, np.ones_like(weighted))
    # The anticorrelated task effects cancel when both fitted rounds are kept.
    result = mean_uncertainty(original, weighted, np.array([[0, 1]] * 3), True)
    assert result["conditional_task"]["ci95"] == [.5, .5]


def test_primary_mean_of_ratios_differs_from_pooled_count_ratio():
    # Round 1 has 0/1 dangerous releases and round 2 has 100/100.
    # Equal averaging over calibration policies targets 1/2, not 100/101.
    cube = np.array([[[0., 1.]], [[100., 100.]]])
    original, weighted = ratio_arrays(cube, ["released", "screened"], "released", "screened", np.array([[1], [1]]))
    result = mean_uncertainty(original, weighted, np.array([[0, 1], [1, 0]]), True)
    assert result["primary_defined"] and result["point_mean"] == .5
    pooled = cube[:, :, 0].sum() / cube[:, :, 1].sum()
    assert pooled == pytest.approx(100 / 101)
    assert abs(result["point_mean"] - pooled) > .49


def test_crossed_intervals_match_explicitly_enumerated_task_and_round_draws():
    original = np.array([.2, .8, .5])
    weighted = np.array([[.1, .7, .4], [.3, .9, .5], [.2, .8, .6], [.4, 1., .7]])
    round_draws = np.array([[0, 0, 2], [1, 1, 1], [2, 0, 1], [2, 2, 0]])
    result = mean_uncertainty(original, weighted, round_draws, True)
    # Each draw keeps its own shared task weights, then averages the selected
    # calibration policies. Values are written out to avoid reusing the code's
    # indexing operation as the test oracle.
    crossed = np.array([(.1 + .1 + .4) / 3, .9, (.6 + .2 + .8) / 3, (.7 + .7 + .4) / 3])
    conditional = np.array([(.1 + .7 + .4) / 3, (.3 + .9 + .5) / 3,
                            (.2 + .8 + .6) / 3, (.4 + 1. + .7) / 3])
    round_only = np.array([(.2 + .2 + .5) / 3, .8, (.5 + .2 + .8) / 3, (.5 + .5 + .2) / 3])
    for name, samples in (("crossed", crossed), ("conditional_task", conditional), ("round_mean", round_only)):
        assert result[name]["defined_draws"] == result[name]["total_draws"] == 4
        assert np.allclose(result[name]["ci95"], np.quantile(samples, [.025, .975], method="linear"))
    sensitivity = result["calibration_sensitivity"]
    assert result["point_mean"] == .5
    assert sensitivity["sample_sd"] == pytest.approx(.3)
    assert sensitivity["min"] == .2 and sensitivity["max"] == .8
    assert "NOT a confidence interval" in sensitivity["interpretation"]


def test_between_task_uncertainty_does_not_vanish_with_many_calibration_rounds():
    # The same uncertain task identity applies to every fitted policy; treating
    # task resampling as independent per round would incorrectly average away
    # this persistent source of variation.
    weights = np.array([[2, 0], [1, 1], [0, 2]] * 100)
    for rounds in (2, 20, 200):
        cube = np.tile(np.array([[[0., 10.], [10., 10.]]]), (rounds, 1, 1))
        original, weighted = ratio_arrays(cube, ["num", "den"], "num", "den", weights)
        draws = np.random.default_rng(81).integers(0, rounds, size=(len(weights), rounds))
        result = mean_uncertainty(original, weighted, draws, True)
        assert result["point_mean"] == .5
        assert result["conditional_task"]["ci95"] == [0., 1.]
        assert result["crossed"]["ci95"] == [0., 1.]
        assert result["round_mean"]["ci95"] == [.5, .5]


def test_identical_policies_have_zero_paired_effect_and_zero_width_intervals():
    cube, columns = _event_cube()
    weights, draws = bootstrap_design(["H0", "H0", "H1", "H1"], 3, 400, 901)
    effects = paired_aggregate(cube, cube.copy(), columns, weights, draws, np.ones(3, dtype=bool))
    assert len(effects) == 3
    for result in effects.values():
        assert result["primary_defined"] and result["point_mean"] == 0
        for kind in ("conditional_task", "round_mean", "crossed"):
            assert result[kind]["ci95"] == [0., 0.]
        assert result["calibration_sensitivity"]["sample_sd"] == 0


def test_paired_policies_reject_different_test_task_denominators():
    first, columns = _event_cube()
    second = first.copy()
    second[0, 2, columns.index("h1_danger_screened")] += 1
    weights, draws = bootstrap_design(["H0", "H0", "H1", "H1"], 3, 20, 82)
    with pytest.raises(ValueError, match="denominators differ"):
        paired_aggregate(first, second, columns, weights, draws, np.ones(3, dtype=bool))


@pytest.mark.parametrize("original,feasible", [([.2, .8], False), ([.2, np.nan], True),
                                               ([.2, np.nan], False)])
def test_infeasible_or_undefined_round_suppresses_primary_and_keeps_sensitivity(original, feasible):
    original = np.asarray(original)
    weighted = np.tile(original, (3, 1))
    result = mean_uncertainty(original, weighted, np.array([[0, 1]] * 3), feasible)
    assert not result["primary_defined"] and result["point_mean"] is None
    assert result["conditional_task"] is result["round_mean"] is result["crossed"] is None
    finite = original[np.isfinite(original)]
    assert result["feasible_defined_round_only_sensitivity_mean"] == np.mean(finite)
    assert result["calibration_sensitivity"]["defined_rounds"] == len(finite)


def test_aggregate_marks_every_primary_metric_undefined_when_one_policy_round_infeasible():
    cube, columns = _event_cube()
    cube[1, :, :] = np.nan
    weights, draws = bootstrap_design(["H0", "H0", "H1", "H1"], 3, 100, 601)
    result = aggregate_policy(cube, columns, weights, draws, np.array([True, False, True]), .05)
    assert result["feasible_rounds"] == 2 and result["total_rounds"] == 3
    assert not result["all_rounds_feasible"]
    for metric in result["metrics"].values():
        assert metric["point_mean"] is None and not metric["primary_defined"]
        assert metric["conditional_task"] is metric["round_mean"] is metric["crossed"] is None
        assert metric["calibration_sensitivity"]["defined_rounds"] == 2
    assert result["metrics"]["H0_withholding_given_screen"]["test_fraction_rounds_above_nominal"] is None


def test_zero_denominator_draws_remain_undefined_and_are_counted():
    # A task bootstrap can miss the only task contributing a risk denominator.
    # Such draws must not be reported as zero risk or as all defined draws.
    cube = np.array([[[1., 2.], [0., 0.]], [[2., 2.], [0., 0.]]])
    weights = np.array([[0, 2], [1, 1], [2, 0]])
    original, weighted = ratio_arrays(cube, ["num", "den"], "num", "den", weights)
    assert np.isnan(weighted[0]).all()
    result = mean_uncertainty(original, weighted, np.array([[0, 1]] * 3), True)
    assert result["point_mean"] == .75
    for kind in ("conditional_task", "crossed"):
        assert result[kind]["defined_draws"] == 2
        assert result[kind]["total_draws"] == 3
        assert result[kind]["ci95"] == [.75, .75]
