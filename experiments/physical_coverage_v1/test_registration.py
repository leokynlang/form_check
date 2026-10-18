"""Source-only registration invariants; no samples, fits or gate outcomes."""
from __future__ import annotations

import json
import shutil

import pytest

from . import registration as reg


@pytest.fixture(scope="module")
def metadata():
    return reg.build_registration()


def test_complete_census_keeps_math_and_execution_obligations_separate(metadata):
    protocol, rows, _, _, summary = metadata
    assert len(rows) == len({row["eq_name"] for row in rows}) == 120
    assert summary["scientific_status"] == {
        "eligible_analytic": 86, "exclusion_unresolved": 31,
        "excluded_nonpositive_target": 2, "excluded_no_positive_input_support": 1,
    }
    assert all(row["analytic_math_status"] == row["scientific_status"] for row in rows)
    assert all(row["physical_source_contract_review_status"].startswith("pending_") for row in rows)
    assert all(row["inclusion_reason"] if row["inclusion"] else row["exclusion_reason"] for row in rows)
    assert all("positive_laurent_sum_tag" in row for row in rows)
    numerical = [row for row in rows if "stable_log_truth_adapter_pending" in row["execution_gates"]]
    assert len(numerical) == 2 and all(row["inclusion"] for row in numerical)
    high_dimensional = [row for row in rows if row["inclusion"] and row["n_vars"] > 4]
    assert len(high_dimensional) == 8
    assert all("diagnostic_and_grid_adapter_pending" in row["execution_gates"] for row in high_dimensional)
    type_pending = [row for row in rows if "source_type_review_pending" in row["execution_gates"]]
    assert [row["eq_name"] for row in type_pending] == ["feynman-i.39.22"]
    assert type_pending[0]["inclusion"]  # no silent discretization or math exclusion
    assert protocol["registration_status"] == "pre_execution_registration_incomplete"
    assert summary["study_samples_generated"] == summary["study_fits_generated"] == summary["diagnostic_outcomes_generated"] == 0


def test_disputed_fine_combination_is_not_a_primary_new_mechanism(metadata):
    _, rows, _, _, summary = metadata
    by_name = {row["eq_name"]: row for row in rows}
    disputed = by_name["feynman-bonus.16"]
    assert disputed["inclusion"] and disputed["fine_exposure"] == "prospective_candidate"
    assert disputed["exposure"] == "development_conservative"
    assert disputed["structural_assignment_status"] == "candidate_pending_structural_review"
    assert disputed["source_exposure_level"] == "metadata_only"
    assert disputed["group_exposure_level"] == "prior_study_outcomes"
    assert not disputed["positive_laurent_sum_tag"]["matches"]
    assert summary["eligible_H1_prospective_groups"] == 0
    assert summary["eligible_H1_disputed_fine_structure_candidates"] == 1
    assert sum(row["source_exposure_level"] == "prior_study_outcomes" for row in rows) == 38


def test_saved_registration_is_read_only_and_rejects_artifact_or_input_drift(tmp_path, monkeypatch):
    source = reg.ROOT / "runs/p03_physical_protocol"
    manifest = json.loads((source / "review_manifest.json").read_text(encoding="utf-8"))
    destination = tmp_path / "runs/p03_physical_protocol"
    destination.mkdir(parents=True)
    names = [*manifest["metadata_artifact_sha256"], "review_manifest.json"]
    for name in names:
        shutil.copy2(source / name, destination / name)
    monkeypatch.setattr(reg, "ROOT", tmp_path)
    before = {name: (reg.digest(destination / name), (destination / name).stat().st_mtime_ns) for name in names}
    assert reg.write_registration(destination) == manifest["summary"]
    after = {name: (reg.digest(destination / name), (destination / name).stat().st_mtime_ns) for name in names}
    assert before == after
    artifact = destination / "candidates.json"
    artifact.write_bytes(artifact.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="Metadata artifact drift"):
        reg.write_registration(destination)
    assert artifact.read_bytes().endswith(b"\n\n")
    shutil.copy2(source / "candidates.json", artifact)
    first_input = next(iter(manifest["source_and_input_sha256"]))
    manifest["source_and_input_sha256"][first_input] = "0" * 64
    (destination / "review_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="Registration input drift"):
        reg.write_registration(destination)
