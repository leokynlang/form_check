"""Group-aware inference and independent-calibration integration checks."""
import json

import numpy as np
import pandas as pd
import pytest

from experiments.fixed_domain import analyze
from experiments.release_validation.analyze import RATIOS
from experiments.release_validation.analyze_v2 import bootstrap_design
from experiments.release_validation.inventory import digest


def test_pool_subsetting_precedes_grouping_for_shared_physical_synthetic_family():
    frame = pd.DataFrame({
        "task": ["null", "physical", "designed", "stress", "old"],
        "klass": ["H0", "H1", "H1", "H1", "H1"],
        "split": ["confirmation", "confirmation", "confirmation", "stress", "development"],
        "source_kind": ["synthetic_null_design", "new_source_positive_restriction",
                        "synthetic_positive_mechanism", "synthetic_positive_mechanism", "old"],
        "group_unit": ["design:null", "structure:shared", "structure:shared", "structure:stress", "structure:old"],
    })
    mask = analyze.pool_mask(frame, "confirmation_physical_sources")
    assert frame.loc[mask, "task"].tolist() == ["null", "physical"]
    counts = pd.DataFrame({"total": [3, 4, 11, 8, 9]}, index=frame.task)
    grouped = analyze.group_counts(counts, frame.set_index("task"),
                                   ["null", "physical"], ["design:null", "structure:shared"])
    assert grouped.loc["structure:shared", "total"] == 4
    assert frame.loc[analyze.pool_mask(frame, "stress"), "task"].tolist() == ["null", "stress"]
    with pytest.raises(ValueError, match="Unknown registered"):
        analyze.pool_mask(frame, "picked_after_results")


def test_equal_group_mean_differs_from_pooled_and_does_not_fill_undefined_groups():
    columns = ["h1_danger_released", "h1_danger_screened"]
    # One large group releases10/100; one small group releases1/1; a third
    # group has no dangerous screened cases. Equal-group risk is.55, not11/101.
    cube = np.array([[[10., 100.], [1., 1.], [0., 0.]],
                     [[20., 100.], [1., 1.], [0., 0.]]])
    weights, draws = bootstrap_design(["H1", "H1", "H1"], 2, 300, 321)
    result = analyze.equal_group_risk(cube, columns, weights, draws, [True, True])
    assert result["point_mean"] == pytest.approx(.575)
    assert result["defined_groups_by_round"] == [2, 2]
    assert result["omitted_groups_by_round"] == [1, 1]
    assert 0 < result["crossed"]["defined_draws"] < 300
    assert "conditional_group" in result and "conditional_task" not in result
    assert result["point_mean"] != pytest.approx((11/101 + 21/101)/2)
    infeasible = analyze.equal_group_risk(cube, columns, weights, draws, [True, False])
    assert infeasible["point_mean"] is None
    assert infeasible["crossed"] is None


def test_group_bootstrap_preserves_family_cluster_and_paired_round_weights():
    classes = ["H0", "H0", "H1", "H1", "H1"]
    weights, rounds = bootstrap_design(classes, 3, 200, 799)
    assert np.all(weights[:, :2].sum(axis=1) == 2)
    assert np.all(weights[:, 2:].sum(axis=1) == 3)
    assert rounds.shape == (200, 3)
    num = np.array([[0, 0, 20, 2, 0], [0, 0, 15, 3, 0], [0, 0, 10, 4, 0]], float)
    den = np.tile([0, 0, 100, 5, 0], (3, 1))
    # Identical policies must have identically paired bootstrap differences.
    original = analyze.divide(num.sum(axis=1), den.sum(axis=1))
    weighted = analyze.divide(weights @ num.T, weights @ den.T)
    result = analyze.mean_uncertainty(original-original, weighted-weighted, rounds, True)
    assert result["point_mean"] == 0
    assert result["crossed"]["ci95"] == [0., 0.]


def test_unavailable_grid_outcome_remains_operational_record():
    frame = pd.DataFrame({
        "task": ["h1"] * 3, "klass": ["H1"] * 3,
        "screen": [True] * 3, "alarm": [False] * 3,
        "decline": [False] * 3, "released": [True] * 3,
        "data_failure": [None] * 3, "fit_failure": [None] * 3,
        "weighted_mean_relative_error": [.2, .05, np.nan],
    })
    counts = analyze.outcome_task_counts(frame, "weighted_mean_relative_error", .1).loc["h1"]
    assert counts.total == 3 and counts.released == 3
    assert counts.outcome_unavailable == 1 and counts.released_unknown_outcome == 1
    assert counts.released_known_outcome == 2 and counts.h1_danger_screened == 1


def _small_run(tmp_path):
    run = tmp_path / "p02_fixed_domain_test"
    stage = run / "stage1"
    stage.mkdir(parents=True)
    registry = [
        {"task": "cal", "klass": "H0", "split": "calibration", "design_group": "calibration",
         "structural_group": "null", "source_kind": "synthetic_null_design"},
        {"task": "null", "klass": "H0", "split": "confirmation", "design_group": "confirmation",
         "structural_group": "null", "source_kind": "synthetic_null_design"},
        {"task": "physical", "klass": "H1", "split": "confirmation", "structural_group": "shared",
         "source_kind": "new_source_positive_restriction"},
        {"task": "designed", "klass": "H1", "split": "confirmation", "structural_group": "shared",
         "source_kind": "synthetic_positive_mechanism"},
        {"task": "old", "klass": "H1", "split": "development", "structural_group": "old",
         "source_kind": "reused_source_restricted_domain"},
        {"task": "stress", "klass": "H1", "split": "stress", "structural_group": "stress",
         "source_kind": "synthetic_positive_mechanism"},
    ]
    (run / "task_registry.json").write_text(json.dumps(registry), encoding="utf-8")
    protocol = {
        "protocol_id": "test-only", "stages": {"stage1": {}}, "master_seed": 786,
        "calibration_rounds": 2, "screen_r2": .99, "budgets": [.05], "primary_budget": .05,
        "danger_thresholds": [.02, .05, .1, .2], "primary_danger_threshold": .1,
        "methods": ["A", "B", "C", "D", "E", "F", "G"],
        "policies": {"D": ["D"], "F": ["F"], "AD_panel": ["A", "B", "C", "D"],
                     "AEF_panel": ["A", "E", "F"], "G_rbf_cv": ["G"]},
        "outcome": {"secondary": ["weighted_p95_relative_error", "max_relative_error", "legacy_endpoint_error"]},
        "group_inference": {"bootstrap_repetitions": 30},
        "evaluation_pools": list(analyze.POOLS),
    }
    rng = np.random.default_rng(177)
    rows = []
    for phase in ["reference", "calibration", "test"]:
        for r in ([-1] if phase == "test" else [0, 1]):
            for task in (registry[1:] if phase == "test" else [registry[0]]):
                for replicate in range(20):
                    row = {"phase": phase, "calibration_round": r, "task": task["task"],
                           "klass": task["klass"], "n": 50, "rel": .03, "replicate": replicate,
                           "r2": .999, "data_failure": None, "fit_failure": None,
                           "dataset_sha256": f"{phase}-{r}-{task['task']}-{replicate}",
                           "weighted_mean_relative_error": .2 if task["klass"] == "H1" else .01,
                           "weighted_p95_relative_error": .3, "max_relative_error": .6,
                           "legacy_endpoint_error": .15}
                    for method in protocol["methods"]:
                        row[f"score_{method}"] = -float(rng.uniform(.001, 1)) if method in "ABCD" else float(rng.normal())
                        row[f"fail_{method}"] = None
                    rows.append(row)
    frame = pd.DataFrame(rows)
    frame.to_parquet(stage / "raw.parquet", index=False)
    (stage / "manifest.json").write_text(json.dumps({"complete": True, "raw_rows": len(frame),
                                                      "raw_sha256": digest(stage / "raw.parquet")}), encoding="utf-8")
    (run / "manifest.json").write_text(json.dumps({"protocol": protocol,
                                                     "task_registry_sha256": digest(run / "task_registry.json")}), encoding="utf-8")
    return run


def test_small_end_to_end_uses_group_fields_and_saves_auditable_outputs(tmp_path, monkeypatch):
    run = _small_run(tmp_path)
    checks = []
    monkeypatch.setattr(analyze, "_verify_frozen", lambda p: checks.append(str(p)))
    summary = analyze.analyze_run(run)
    assert len(checks) == 2
    pool = summary["strata"][0]["pools"]["confirmation_all"]
    assert pool["H0_tasks"] == 1 and pool["H1_tasks"] == 2
    assert pool["H0_design_groups"] == 1 and pool["H1_structure_groups"] == 1
    record = pool["budgets"]["0.05"]["policies"]["F"]
    risk = record["metrics"]["dangerous_H1_release_given_screen"]
    assert risk["denominator_contributing_test_groups_by_round"] == [1, 1]
    assert "conditional_group" in risk
    assert "conditional_task" not in json.dumps(summary)
    groups = pd.read_parquet(run / "stage1/group_counts.parquet")
    physical = groups.loc[(groups.pool == "confirmation_physical_sources") &
                          (groups.group_unit == "structure:shared")]
    assert physical.total.eq(20).all()  # never40 from accidentally included designed task
    decisions = pd.read_parquet(run / "stage1/decisions.parquet")
    assert len(decisions) == 5 * 20 * 5 * 2  # unique test bank, five policies, two calibrations
    assert not decisions.duplicated(["task", "replicate", "n", "rel", "policy", "budget", "calibration_round"]).any()
    assert set(pd.read_csv(run / "stage1/primary_metrics.csv").metric) == set(RATIOS)
    assert (run / "stage1/report.md").exists()
    assert (run / "stage1/policies.json").exists()


def test_analysis_rejects_changed_raw_data_before_policy_fitting(tmp_path, monkeypatch):
    run = _small_run(tmp_path)
    monkeypatch.setattr(analyze, "_verify_frozen", lambda _: None)
    path = run / "stage1/raw.parquet"
    raw = pd.read_parquet(path)
    raw.loc[0, "r2"] = .1
    raw.to_parquet(path, index=False)
    with pytest.raises(ValueError, match="changed raw data"):
        analyze.analyze_run(run)


def test_infeasible_policy_is_preserved_and_never_given_primary_risk(tmp_path, monkeypatch):
    run = _small_run(tmp_path)
    monkeypatch.setattr(analyze, "_verify_frozen", lambda _: None)
    path = run / "stage1/raw.parquet"
    raw = pd.read_parquet(path)
    failing = raw.phase.eq("calibration")
    raw.loc[failing, "fail_G"] = "synthetic_decline"
    raw.loc[failing, "score_G"] = np.nan
    raw.to_parquet(path, index=False)
    execution_path = run / "stage1/manifest.json"
    execution = json.loads(execution_path.read_text(encoding="utf-8"))
    execution["raw_sha256"] = digest(path)
    execution_path.write_text(json.dumps(execution), encoding="utf-8")
    result = analyze.analyze_run(run)
    policy = result["strata"][0]["pools"]["confirmation_all"]["budgets"]["0.05"]["policies"]["G_rbf_cv"]
    assert policy["feasible_rounds"] == 0 and not policy["all_rounds_feasible"]
    assert policy["metrics"]["dangerous_H1_release_given_screen"]["point_mean"] is None
    assert policy["metrics"]["dangerous_H1_release_given_screen"]["crossed"] is None
    summary_text = (run / "stage1/summary.json").read_text(encoding="utf-8")
    assert "NaN" not in summary_text and "Infinity" not in summary_text
    groups = pd.read_parquet(run / "stage1/group_counts.parquet")
    assert not groups.loc[groups.policy.eq("G_rbf_cv"), "feasible"].any()
    assert groups.loc[groups.policy.eq("G_rbf_cv"), "h1_danger_released"].isna().all()
