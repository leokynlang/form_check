"""Guard against counting the two registered stages as repeated estimates."""
import copy

import numpy as np
import pandas as pd
import pytest

from experiments.reporting.p02_complete import combine_summaries, coverage


def _summary(strata, rows):
    return {"protocol_id": "same", "primary_outcome": "weighted_mean_relative_error",
            "primary_danger_threshold": .1, "primary_budget": .05,
            "calibration_rounds": 10, "group_inference": {"resampling": "groups"},
            "point_estimator": "mean of round ratios", "scope": "conditional",
            "uncertainty_names": {}, "raw_sha256": "unchanged", "raw_rows": rows, "strata": strata}


def test_combination_requires_exact_disjoint_grid_and_retains_each_result():
    strata = [{"n": n, "rel": r, "pools": {"arbitrary_saved_stat": n+r}}
              for n in (50, 200, 1000) for r in (.005, .01, .03, .05)]
    early = [s for s in strata if s["n"] in (50, 200) and s["rel"] == .03]
    later = [s for s in strata if s not in early]
    result = combine_summaries([("stage1", _summary(early, 41040)), ("stage2", _summary(later, 205200))])
    assert result["strata"] == strata
    assert result["raw_rows"] == 246240
    assert result["source_stages"]["stage1"]["raw_rows"] == 41040
    assert "pooled_effect" not in result
    with pytest.raises(ValueError, match="12 disjoint"):
        combine_summaries([("stage1", _summary(early, 41040)), ("stage2", _summary([*later, early[0]], 205200))])


def test_combination_rejects_changed_inference_rules():
    a = _summary([], 0)
    b = copy.deepcopy(a)
    b["primary_danger_threshold"] = .2
    with pytest.raises(ValueError, match="identical estimands"):
        combine_summaries([("stage1", a), ("stage2", b)])


def test_physical_coverage_separates_danger_from_screen_and_secondary():
    frame = pd.DataFrame({"n": [50]*5, "rel": [.03]*5, "split": ["confirmation"]*5,
                          "klass": ["H0", "H0", "H1", "H1", "H1"],
                          "source_kind": ["null", "null", "new_source_positive_restriction", "new_source_positive_restriction", "synthetic_positive_mechanism"],
                          "group_unit": ["design:one", "design:two", "structure:g", "structure:l", "structure:s"],
                          "r2": [.999, .95, .95, .999, .999],
                          "weighted_mean_relative_error": [.01, .01, .4, .03, .2],
                          "weighted_p95_relative_error": [.01, .01, .6, .15, .3]})
    result = coverage(frame).iloc[0]
    assert result.confirmation_h0_screen_rate == .5
    assert result.confirmation_dangerous_screened == 1
    assert result.confirmation_danger_structure_groups == 1
    assert result.physical_dangerous_screened == 0
    assert result.physical_danger_structure_groups == 0
    assert result.physical_p95_dangerous_screened == 1
    assert result.physical_p95_danger_structure_groups == 1
    assert np.isfinite(result.physical_total)
