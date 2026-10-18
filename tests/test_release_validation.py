"""Scientific invariants: leakage, budget accounting, parity and missing outcomes."""
import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm

from common import form_check_with_p, reset_test
from experiments.release_validation.analyze import per_task_counts
from experiments.release_validation.calibration import (
    apply_policy,
    calibrate,
    choose_threshold,
    fit_reference,
    policy_scores,
)
from experiments.release_validation.gates import rbf_cv_score, run_gates
from experiments.release_validation.runner import fit_and_endpoint, sample_seed
from gate_compare import gate_model_compare
from scenarios import h0_multivar


def test_budget_prices_declines_and_ties():
    fitted = choose_threshold([0.1, .2, .3, .4, .4, .5, .6, .8, .9, np.nan],
                              [False] * 9 + [True], .2)
    assert fitted["feasible"] and fitted["alarms"] + fitted["declines"] <= 2
    assert fitted["threshold"] == .8
    assert not choose_threshold([np.nan, np.nan, .2], [True, True, False], .5)["feasible"]
    assert not choose_threshold([], [], .05)["feasible"]
    assert choose_threshold([.9] * 100, [False] * 100, .05)["alarms"] == 0


def test_reference_and_threshold_reject_wrong_phase_or_H1():
    frame = pd.DataFrame({"klass": ["H1"], "phase": ["reference"], "score_A": [.5], "fail_A": [None]})
    with pytest.raises(ValueError):
        fit_reference(frame, ["A"])
    frame["klass"], frame["phase"] = "H0", "test"
    with pytest.raises(ValueError):
        calibrate(frame, ["A"], {"A": [.1]}, .05)


def test_joint_budget_is_calibrated_as_panel_and_missing_member_withholds():
    frame = pd.DataFrame({"klass": ["H0"] * 10, "phase": ["calibration"] * 10,
                          "score_A": np.linspace(0, 1, 10), "score_D": np.linspace(1, 0, 10),
                          "fail_A": [None] * 10, "fail_D": [None] * 10})
    ref = {"A": np.linspace(0, 1, 20).tolist(), "D": np.linspace(0, 1, 20).tolist()}
    fitted = calibrate(frame, ["A", "D"], ref, .2)
    _, alarm, decline = apply_policy(frame, ["A", "D"], ref, fitted)
    assert np.sum(alarm | decline) <= 2
    frame.loc[0, "fail_D"] = "sentinel_failure"
    frame.loc[0, "score_D"] = .01  # Even a finite quiet score cannot override a failure.
    _, alarm, decline = apply_policy(frame, ["A", "D"], ref, fitted)
    assert decline[0] and not alarm[0]
    _, missing = policy_scores(frame, ["D"], {"D": []})
    assert missing.all()


def test_seed_namespaces_are_disjoint_and_deterministic():
    seeds = [sample_seed(20261018, phase, task, n, .03, rep)
             for phase in ("reference", "calibration", "test")
             for task in ("a", "b") for n in (50, 200) for rep in range(10)]
    assert len(set(seeds)) == len(seeds)
    assert min(seeds) > 2**40 and max(seeds) < 2**63
    assert sample_seed(1, "test", "a", 50, .03, 0) == sample_seed(1, "test", "a", 50, .03, 0)


def test_v1_empirical_cdf_ceiling_exposes_low_budget_power_loss():
    # Document the v1 pilot pathology: all observations beyond the reference
    # maximum tie. A strict threshold at that ceiling cannot reject any of them.
    frame = pd.DataFrame({"klass": ["H0"] * 100, "phase": ["calibration"] * 100,
                          "score_A": [1000.] * 100, "fail_A": [None] * 100})
    reference = {"A": [0., 1., 2.]}
    fitted = calibrate(frame, ["A"], reference, .01)
    scores, alarm, decline = apply_policy(frame, ["A"], reference, fitted)
    assert fitted["threshold"] == 3 / 4
    assert np.all(scores == 3 / 4) and not alarm.any() and not decline.any()


def test_adapter_matches_registered_A_B_E_F_and_independent_HC3():
    data, _ = h0_multivar(123, 200, .03)
    x, y = np.log(data[["x1", "x2"]].to_numpy()), np.log(data.target.to_numpy())
    got = run_gates(data, ["x1", "x2"], 123)
    a = form_check_with_p(data, "target", ["x1", "x2"], 123)
    b = reset_test(x, y, variant="fitted", vcov="classical")
    ef = gate_model_compare(x, y, 123)
    assert got["score_A"] == -a["p_fc"]
    assert got["score_B"] == -b["p"]
    assert got["score_E"] == ef["D_cv"] and got["score_F"] == ef["delta_aicc"]
    base = np.column_stack([np.ones(len(y)), x])
    yhat = base @ np.linalg.lstsq(base, y, rcond=None)[0]
    for method, added in [("C", np.column_stack([yhat**2, yhat**3])),
                           ("D", np.column_stack([x[:, j]**k for j in range(2) for k in (2, 3)]))]:
        design = np.column_stack([base, added])
        model = sm.OLS(y, design).fit(cov_type="HC3")
        restriction = np.zeros((added.shape[1], design.shape[1]))
        restriction[:, -added.shape[1]:] = np.eye(added.shape[1])
        p = float(model.wald_test(restriction, use_f=True, scalar=True).pvalue)
        assert np.isclose(-got[f"score_{method}"], p, atol=1e-9, rtol=1e-7)
        assert got[f"fail_{method}"] is None


def test_rbf_training_only_preprocessing_is_affine_invariant_and_reproducible():
    rng = np.random.default_rng(93)
    x = rng.normal(size=(100, 2))
    y = x[:, 0] + .2 * x[:, 1]**2 + rng.normal(scale=.05, size=100)
    score = rbf_cv_score(x, y, 17)
    assert score == rbf_cv_score(x, y, 17)
    assert np.isclose(score, rbf_cv_score(x * [3, 7] + [30, -17], y, 17), atol=1e-10)


def test_endpoint_matches_explicit_power_prediction():
    data, _ = h0_multivar(51, 50, .03)
    result = fit_and_endpoint(data, ["x1", "x2"], lambda x: 1.2*x[0]**.5*x[1]**.9)
    assert result["fit_failure"] is None and result["outcome_failure"] is None
    assert np.isfinite(result["r2"]) and result["legacy_endpoint_error"] >= 0
    invalid = fit_and_endpoint(data, ["x1", "x2"], lambda x: -1)
    assert np.isnan(invalid["legacy_endpoint_error"]) and invalid["outcome_failure"]


def test_unknown_outcome_is_not_counted_safe_or_removed_from_operations():
    frame = pd.DataFrame({"task": ["a", "a"], "klass": ["H1", "H1"],
                          "screen": [True, True], "alarm": [False, False],
                          "decline": [False, False], "released": [True, True],
                          "legacy_endpoint_error": [.2, np.nan], "data_failure": [None, None],
                          "fit_failure": [None, None]})
    frame.attrs["danger_threshold"] = .1
    counts = per_task_counts(frame).loc["a"]
    assert counts["released"] == 2
    assert counts["released_known_outcome"] == 1
    assert counts["h1_danger_released"] == 1
    assert counts["h1_safe_screened"] == 0
    assert counts["released_unknown_outcome"] == 1
