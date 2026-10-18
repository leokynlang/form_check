"""V2 scientific invariants: unbounded tails, honest budgets, and data isolation."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from experiments.release_validation import calibration as v1
from experiments.release_validation.calibration_v2 import (
    apply_policy,
    calibrate,
    fit_reference,
    policy_scores,
    transform_raw,
)
from experiments.release_validation.runner import sample_seed as sample_seed_v1
from experiments.release_validation.runner_v2 import sample_seed


def _frame(values, phase="reference", klass="H0"):
    """Construct gate records, preserving the stored convention score_A-D=-p."""
    n = len(next(iter(values.values())))
    frame = pd.DataFrame({"klass": [klass] * n, "phase": [phase] * n})
    for method, scores in values.items():
        frame[f"score_{method}"] = np.asarray(scores, dtype=float)
        frame[f"fail_{method}"] = [None] * n
    return frame


@pytest.mark.parametrize("method", list("ABCD"))
def test_p_transform_preserves_evidence_order_and_records_numeric_floor(method):
    tiny = np.finfo(float).tiny
    p = np.array([1., .5, .01, 1e-20, tiny, tiny / 2, 0.])
    transformed = transform_raw(-p, method)
    assert np.all(np.isfinite(transformed))
    assert transformed[0] == 0
    assert np.all(np.diff(transformed) >= 0)
    assert np.all(np.diff(transformed[:5]) > 0)
    # Floating-point zeros and subnormal p values saturate only at the explicit
    # numeric floor, never at the finite empirical reference maximum.
    assert transformed[-1] == transformed[-2] == transformed[-3]
    assert transformed[-1] > 300
    reference = fit_reference(_frame({method: -p}), [method])[method]
    assert reference["reference_p_zero_count"] == 1
    assert reference["floating_floor_count"] >= 2


def test_permutation_resolution_floor_is_distinct_from_floating_point_underflow():
    reference = fit_reference(_frame({"A": [-.9, -.1, -.001, -.001, 0.]}), ["A"])["A"]
    assert reference["reference_p_floor_count"] == 2
    assert reference["reference_p_zero_count"] == 1
    assert reference["floating_floor_count"] == 1


@pytest.mark.parametrize("method", list("ABCD"))
def test_invalid_p_values_decline_instead_of_becoming_strong_or_quiet_evidence(method):
    reference = fit_reference(_frame({method: [-.8, -.5, -.2, -.01]}), [method])
    invalid = np.array([.01, -1.01, np.inf, -np.inf, np.nan])
    scores, decline = policy_scores(_frame({method: invalid}, phase="test"), [method], reference)
    assert decline.all()
    assert np.isnan(scores).all()
    # p=0 is a separately reported numerical floor, not a missing gate outcome.
    scores, decline = policy_scores(_frame({method: [0.]}, phase="test"), [method], reference)
    assert not decline.any() and np.isfinite(scores).all()


@pytest.mark.parametrize("method", list("EFG"))
def test_comparison_score_tail_remains_strictly_ordered_beyond_reference(method):
    reference = fit_reference(_frame({method: np.arange(-5., 6.)}), [method])
    extreme = _frame({method: [10., 100., 1000.]}, phase="test")
    scores, decline = policy_scores(extreme, [method], reference)
    assert not decline.any()
    assert np.all(np.diff(scores) > 0)
    assert scores[-1] > 10 * scores[0]


@pytest.mark.parametrize("method", list("ABCD"))
def test_p_evidence_beyond_smallest_reference_p_remains_ordered(method):
    reference = fit_reference(_frame({method: [-.9, -.5, -.1, -.01]}), [method])
    frame = _frame({method: [-1e-3, -1e-10, -1e-100]}, phase="test")
    scores, decline = policy_scores(frame, [method], reference)
    assert not decline.any() and np.all(np.diff(scores) > 0)


def test_v1_cdf_all_quiet_pathology_is_rejected_in_v2():
    reference_frame = _frame({"E": np.arange(10.)})
    calibration_frame = _frame({"E": np.arange(100., 200.)}, phase="calibration")
    test_frame = _frame({"E": [201., 1000.]}, phase="test")
    old_reference = v1.fit_reference(reference_frame, ["E"])
    old_fit = v1.calibrate(calibration_frame, ["E"], old_reference, .01)
    _, old_alarm, _ = v1.apply_policy(test_frame, ["E"], old_reference, old_fit)
    assert not old_alarm.any()  # All beyond-reference values tie under v1.
    reference = fit_reference(reference_frame, ["E"])
    fitted = calibrate(calibration_frame, ["E"], reference, .01)
    _, calibration_alarm, calibration_decline = apply_policy(
        calibration_frame, ["E"], reference, fitted
    )
    scores, alarm, decline = apply_policy(test_frame, ["E"], reference, fitted)
    assert fitted["feasible"]
    assert int(np.sum(calibration_alarm | calibration_decline)) == 1
    assert alarm.all() and not decline.any()
    assert scores[1] > scores[0] > fitted["threshold"]


@pytest.mark.parametrize("method", list("EFG"))
def test_positive_affine_rescaling_does_not_change_policy_decisions(method):
    reference_frame = _frame({method: [-8., -2., -1., 0., 1., 2., 9.]})
    calibration_frame = _frame({method: np.linspace(-3., 4., 50)}, phase="calibration")
    test_frame = _frame({method: [-4., 0., 3., 5., 50.]}, phase="test")
    reference = fit_reference(reference_frame, [method])
    fitted = calibrate(calibration_frame, [method], reference, .1)
    original = apply_policy(test_frame, [method], reference, fitted)
    for frame in (reference_frame, calibration_frame, test_frame):
        frame[f"score_{method}"] = 17. * frame[f"score_{method}"] + 113.
    rescaled_reference = fit_reference(reference_frame, [method])
    rescaled_fit = calibrate(calibration_frame, [method], rescaled_reference, .1)
    rescaled = apply_policy(test_frame, [method], rescaled_reference, rescaled_fit)
    assert np.allclose(original[0], rescaled[0], rtol=1e-12, atol=1e-12)
    assert np.array_equal(original[1], rescaled[1])
    assert np.array_equal(original[2], rescaled[2])


@pytest.mark.parametrize("phase,klass", [("test", "H0"), ("calibration", "H0"),
                                          ("reference", "H1")])
def test_reference_rejects_test_calibration_and_alternative_records(phase, klass):
    with pytest.raises(ValueError):
        fit_reference(_frame({"E": [0., 1., 2.]}, phase=phase, klass=klass), ["E"])


@pytest.mark.parametrize("phase,klass", [("test", "H0"), ("reference", "H0"),
                                          ("calibration", "H1")])
def test_threshold_rejects_test_reference_and_alternative_records(phase, klass):
    reference = fit_reference(_frame({"E": [0., 1., 2.]}), ["E"])
    with pytest.raises(ValueError):
        calibrate(_frame({"E": [0., 1., 2.]}, phase=phase, klass=klass), ["E"], reference, .1)


@pytest.mark.parametrize("values", [[], [3.] * 20])
def test_absent_or_constant_reference_declines_and_cannot_claim_feasible_budget(values):
    reference = fit_reference(_frame({"E": values}), ["E"])
    assert not reference["E"]["usable"]
    calibration_frame = _frame({"E": np.arange(100.)}, phase="calibration")
    scores, decline = policy_scores(calibration_frame, ["E"], reference)
    assert decline.all() and np.isnan(scores).all()
    fitted = calibrate(calibration_frame, ["E"], reference, .1)
    assert not fitted["feasible"]
    with pytest.raises(ValueError):
        apply_policy(calibration_frame, ["E"], reference, fitted)


def test_discrete_nonconstant_reference_uses_nonzero_fallback_scale():
    # Both IQR and MAD vanish, but a genuinely nonconstant reference must remain
    # usable via the preregistered standard-deviation fallback.
    reference = fit_reference(_frame({"G": [0.] * 19 + [1.]}), ["G"])
    assert reference["G"]["usable"] and reference["G"]["scale"] > 0
    scores, decline = policy_scores(_frame({"G": [0., 1., 10.]}, phase="test"), ["G"], reference)
    assert not decline.any() and np.all(np.diff(scores) > 0)


def test_reference_excludes_failed_or_invalid_values_from_scale_estimation():
    clean = _frame({"D": [-.9, -.6, -.3, -.01]})
    contaminated = _frame({"D": [-.9, -.6, -.3, -.01, 1., -10., -1e-100, np.nan]})
    contaminated.loc[6, "fail_D"] = "declared_gate_failure"
    clean_reference = fit_reference(clean, ["D"])
    reference = fit_reference(contaminated, ["D"])
    assert reference["D"]["n"] == 4
    assert reference["D"]["center"] == clean_reference["D"]["center"]
    assert reference["D"]["scale"] == clean_reference["D"]["scale"]


def test_joint_calibration_prices_declines_and_uses_panel_maximum():
    reference = fit_reference(_frame({"E": np.arange(50.), "F": np.arange(50.)}), ["E", "F"])
    frame = _frame({"E": np.arange(100.), "F": np.arange(99., -1., -1.)}, phase="calibration")
    frame.loc[:4, "fail_F"] = "unavailable"
    score_e, _ = policy_scores(frame, ["E"], reference)
    score_f, _ = policy_scores(frame, ["F"], reference)
    panel_score, panel_decline = policy_scores(frame, ["E", "F"], reference)
    assert panel_decline.sum() == 5
    assert np.allclose(panel_score[5:], np.maximum(score_e[5:], score_f[5:]))
    fitted = calibrate(frame, ["E", "F"], reference, .1)
    _, alarm, decline = apply_policy(frame, ["E", "F"], reference, fitted)
    assert fitted["feasible"] and not np.any(alarm & decline)
    assert int(np.sum(alarm | decline)) <= 10
    frame.loc[:10, "fail_F"] = "unavailable"
    assert not calibrate(frame, ["E", "F"], reference, .1)["feasible"]


def test_negative_standardized_calibration_scores_do_not_use_nonnegative_sentinel():
    reference = fit_reference(_frame({"E": np.arange(100., 200.)}), ["E"])
    calibration_frame = _frame({"E": np.arange(10.)}, phase="calibration")
    fitted = calibrate(calibration_frame, ["E"], reference, .1)
    _, alarm, decline = apply_policy(calibration_frame, ["E"], reference, fitted)
    assert fitted["threshold"] < -1
    assert alarm.sum() == 1 and not decline.any()
    _, test_alarm, test_decline = apply_policy(_frame({"E": [0., 1000.]}, phase="test"), ["E"], reference, fitted)
    assert test_alarm.tolist() == [False, True] and not test_decline.any()


def test_gate_failure_cannot_be_overridden_by_finite_quiet_score():
    reference = fit_reference(_frame({"E": [0., 1., 2., 3.], "G": [0., 1., 2., 3.]}), ["E", "G"])
    frame = _frame({"E": [0., 0.], "G": [0., np.inf]}, phase="test")
    frame.loc[0, "fail_G"] = "numerical_failure"
    scores, decline = policy_scores(frame, ["E", "G"], reference)
    assert decline.all() and np.isnan(scores).all()


def test_finite_raw_score_with_nonfinite_standardization_declines():
    reference = fit_reference(_frame({"E": np.array([-1., 0., 1.]) * 1e-12}), ["E"])
    scores, decline = policy_scores(_frame({"E": [0., np.finfo(float).max]}, phase="test"), ["E"], reference)
    assert not decline[0] and np.isfinite(scores[0])
    assert decline[1] and np.isnan(scores[1])


def test_independent_calibration_rounds_cannot_be_pooled():
    reference_frame = _frame({"E": np.arange(10.)})
    reference_frame["calibration_round"] = [0] * 5 + [1] * 5
    with pytest.raises(ValueError):
        fit_reference(reference_frame, ["E"])
    reference_frame["calibration_round"] = 0
    reference = fit_reference(reference_frame, ["E"])
    frame = _frame({"E": np.arange(10.)}, phase="calibration")
    frame["calibration_round"] = [0] * 5 + [1] * 5
    with pytest.raises(ValueError):
        calibrate(frame, ["E"], reference, .1)


def test_repeated_calibration_seed_namespaces_are_disjoint_and_reproducible():
    seeds = [sample_seed(20261018, phase, round_id, task, n, .03, replicate)
             for phase in ("reference", "calibration")
             for round_id in range(5) for task in ("a", "b")
             for n in (50, 200) for replicate in range(4)]
    seeds += [sample_seed(20261018, "test", -1, task, n, .03, replicate)
              for task in ("a", "b") for n in (50, 200) for replicate in range(4)]
    assert len(set(seeds)) == len(seeds)
    assert min(seeds) > 2**40 and max(seeds) < 2**63
    old_seeds = {sample_seed_v1(20261018, phase, task, n, .03, replicate)
                 for phase in ("reference", "calibration", "test")
                 for task in ("a", "b") for n in (50, 200) for replicate in range(4)}
    assert not old_seeds.intersection(seeds)
    assert sample_seed(20261018, "reference", 0, "a", 50, .03, 0) == sample_seed(
        20261018, "reference", 0, "a", 50, .03, 0
    )
