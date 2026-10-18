"""Audit the read-only 12-stratum P0-2 assembly and historical preservation.

The two stages are portions of one frozen protocol. Their statistical estimates
stay in separate (n, noise) strata; this script does not turn execution stages,
calibration rounds, parameter windows or repeated rows into new physical groups.
No experimental implementation is imported and no single-gate statistic is
recomputed. Independent per-stage audits certify the underlying calculations.
"""
from __future__ import annotations

import argparse
import ast
import json
import math
import re
from collections import Counter
from pathlib import Path

import mpmath as mp
import numpy as np
import pandas as pd

try:
    from .audit_fixed_domain import ROOT, digest, read_json, require, tree_hashes
except ImportError:
    from audit_fixed_domain import ROOT, digest, read_json, require, tree_hashes


def assert_same(actual, expected, label):
    require(actual == expected, label)


def audit_preservation(run):
    before = read_json(run / "stage2_preservation_before.json")
    actual_stage1 = tree_hashes(run / "stage1")
    assert_same(actual_stage1, before["stage1_files"], "Stage1 file tree changed during stage2/reporting")
    paper = ROOT / "paper/research_R2"
    assert_same(tree_hashes(paper), before["paper_R2_files"], "Original R2 file tree changed")
    original = read_json(run / "manifest.json")
    result = {"stage1_files": len(actual_stage1), "paper_R2_files": len(before["paper_R2_files"]), "p01": {}}
    for name, expected in original["preservation"].items():
        assert_same(tree_hashes(ROOT / name), expected, f"Historical P0-1 tree changed: {name}")
        result["p01"][name] = len(expected)
    return result


def audit_stage_certificates(run):
    root = read_json(run / "manifest.json")
    result = {}
    for stage in ("stage1", "stage2"):
        directory = run / stage
        audit = read_json(directory / "independent_audit.json")
        execution = read_json(directory / "manifest.json")
        summary = read_json(directory / "summary.json")
        require(audit["passed"] and audit["analysis"]["status"] == "passed", f"Missing successful independent {stage} audit")
        require(execution["complete"] and execution["fingerprint"] == root["fingerprint"], f"{stage} execution is incomplete or from another protocol")
        expected_raw = digest(directory / "raw.parquet")
        require(audit["raw"]["raw_sha256"] == execution["raw_sha256"] == summary["raw_sha256"] == expected_raw, f"{stage} source raw hash inconsistency")
        require(summary["protocol_id"] == root["protocol"]["protocol_id"] and summary["stage"] == stage, f"{stage} summary identity")
        require(audit["provenance"]["protocol_id"] == root["protocol"]["protocol_id"], f"{stage} audit protocol identity")
        if stage == "stage2":
            require(bool(audit.get("audited_artifact_sha256")), "Stage2 audit must pin audited analysis artifacts")
        for name, expected in audit.get("audited_artifact_sha256", {}).items():
            require(digest(directory / name) == expected, f"Audited artifact changed after {stage} certification: {name}")
        result[stage] = {"raw_sha256": expected_raw, "summary_sha256": digest(directory / "summary.json"),
                         "independent_audit_sha256": digest(directory / "independent_audit.json"),
                         "raw_rows": execution["raw_rows"], "strata": len(summary["strata"]),
                         "analysis_checks": audit["analysis"]}
    return result


def audit_disjoint_raw(run, protocol):
    keys = ["phase", "calibration_round", "task", "n", "rel", "replicate"]
    columns = [*keys, "base_seed", "dataset_sha256", "klass", "split", "group_unit"]
    frames = {stage: pd.read_parquet(run / stage / "raw.parquet", columns=columns) for stage in ("stage1", "stage2")}
    first, second = frames["stage1"], frames["stage2"]
    seeds = set(first.base_seed) & set(second.base_seed)
    hashes = set(first.dataset_sha256.dropna()) & set(second.dataset_sha256.dropna())
    require(not seeds and not hashes, "Stage1/stage2 share seeds or generated datasets")
    planned = {(n, float(rel)) for n in protocol["stages"]["stage2"]["n_grid"] for rel in protocol["stages"]["stage2"]["noise_grid"]}
    strata = {stage: set(zip(frame.n, frame.rel, strict=True)) for stage, frame in frames.items()}
    require(not strata["stage1"] & strata["stage2"], "The stages repeat an existing n/noise stratum")
    assert_same(strata["stage1"] | strata["stage2"], planned, "The two stages do not cover exactly the 12 registered strata")
    combined = pd.concat(frames.values(), ignore_index=True)
    require(not combined.duplicated(keys).any(), "Repeated complete-grid experiment identities")
    require(combined.base_seed.nunique() == len(combined), "Repeated seeds inside complete grid")
    require(not combined.dataset_sha256.dropna().duplicated().any(), "Repeated dataset hashes inside complete grid")
    require(not combined.split.eq("reserve").any(), "Reserve tasks have outcomes in the complete grid")
    # Each stratum must retain the same task/class/group identities. Repetition
    # over sample size or noise does not increase the structural library size.
    identity = ["task", "klass", "split", "group_unit"]
    expected_identities = None
    for _, part in combined.groupby(["n", "rel"], sort=True):
        observed = set(part[identity].itertuples(index=False, name=None))
        if expected_identities is None:
            expected_identities = observed
        assert_same(observed, expected_identities, "Task/group identities changed across strata")
    test = combined[combined.phase.eq("test")]
    return {"raw_rows": len(combined), "test_rows": len(test), "strata": len(planned),
            "stage_strata": {stage: [list(value) for value in sorted(values)] for stage, values in strata.items()},
            "stage_seed_overlap": len(seeds), "stage_dataset_hash_overlap": len(hashes),
            "unique_seeds": int(combined.base_seed.nunique()),
            "unique_nonmissing_dataset_hashes": int(combined.dataset_sha256.nunique()),
            "registered_test_tasks": int(test.task.nunique()),
            "global_test_H1_structural_groups": int(test.loc[test.klass.eq("H1"), "group_unit"].nunique()),
            "confirmation_H1_structural_groups": int(test.loc[test.klass.eq("H1") & test.split.eq("confirmation"), "group_unit"].nunique()),
            "interpretation": "One fixed library evaluated in 12 separate sample-size/noise strata; stages are execution partitions, not independent replications."}


def audit_summary(run, complete, raw_audit):
    combined = read_json(complete / "summary.json")
    source = {stage: read_json(run / stage / "summary.json") for stage in ("stage1", "stage2")}
    expected = sorted([entry for saved in source.values() for entry in saved["strata"]], key=lambda s: (s["n"], s["rel"]))
    require(combined["stage"] == "complete_grid", "Combined summary stage identity")
    require(combined["raw_rows"] == raw_audit["raw_rows"], "Combined summary raw row count")
    require(combined["protocol_id"] == read_json(run / "manifest.json")["protocol"]["protocol_id"], "Combined summary protocol identity")
    expected_stages = {stage: {"raw_sha256": value["raw_sha256"], "raw_rows": value["raw_rows"]} for stage, value in source.items()}
    assert_same(combined["source_stages"], expected_stages, "Combined summary changes stage source references")
    assert_same(combined["strata"], expected, "Combined summary alters per-stratum estimates or intervals")
    strata = {(entry["n"], entry["rel"]) for entry in combined["strata"]}
    require(len(strata) == len(expected) == raw_audit["strata"], "Duplicate/missing summary strata")
    for key in ("primary_outcome", "primary_danger_threshold", "primary_budget", "calibration_rounds", "group_inference", "point_estimator", "uncertainty_names"):
        require(source["stage1"][key] == source["stage2"][key] == combined[key], f"Combined summary changes statistical metadata: {key}")
    return {"exact_source_strata": len(expected), "per_stratum_estimates_preserved": True,
            "cross_stratum_reestimation_performed": False}


def audit_csvs(run, complete):
    result = {}
    for name in ("primary_metrics.csv", "sensitivity_metrics.csv", "round_metrics.csv", "reference_scales.csv"):
        actual = pd.read_csv(complete / name, float_precision="round_trip")
        require("source_stage" in actual, f"{name}: missing source-stage provenance")
        sources = []
        for stage in ("stage1", "stage2"):
            frame = pd.read_csv(run / stage / name, float_precision="round_trip")
            frame["source_stage"] = stage
            sources.append(frame)
        expected = pd.concat(sources, ignore_index=True)
        require(set(actual.columns) == set(expected.columns), f"{name}: added/dropped columns")
        actual = actual[expected.columns]
        pd.testing.assert_frame_equal(actual.reset_index(drop=True), expected.reset_index(drop=True),
                                      check_dtype=False, check_exact=False, rtol=2e-15, atol=0,
                                      obj=f"Read-only concatenation: {name}")
        result[name] = {"rows": len(actual), "sha256": digest(complete / name), "source_rows_preserved": True}
    return result


def audit_policies(run, complete):
    expected = [entry for stage in ("stage1", "stage2") for entry in read_json(run / stage / "policies.json")]
    actual = read_json(complete / "policies.json")
    assert_same(actual, expected, "Combined policies change source calibrations or ordering")
    keys = [(p["n"], p["rel"], p["budget"], p["policy"], p["calibration_round"]) for p in actual]
    require(len(keys) == len(set(keys)), "Repeated combined policy identities")
    return {"policy_configurations": len(actual), "exact_source_equality": True}


def _digests(value):
    """Read source/output hash fields without relying on reporting key labels."""
    if isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value):
        return {value}
    if isinstance(value, dict):
        return set().union(*(_digests(v) for v in value.values())) if value else set()
    if isinstance(value, list):
        return set().union(*(_digests(v) for v in value)) if value else set()
    return set()


def audit_manifest(run, complete, certificates, csvs):
    manifest = read_json(complete / "manifest.json")
    root = read_json(run / "manifest.json")
    require(manifest["complete"] is True and manifest["stage"] == "complete_grid", "Combined manifest incomplete or wrong stage")
    require(manifest["protocol_fingerprint"] == root["fingerprint"], "Combined manifest refers to a different frozen protocol")
    require(manifest["dataset_count"] == sum(records["raw_rows"] for records in certificates.values()), "Combined manifest dataset cardinality")
    require(manifest["strata"] == sum(records["strata"] for records in certificates.values()), "Combined manifest stratum cardinality")
    require(manifest["test_count"] == sum(read_json(run / stage / "independent_audit.json")["raw"]["test_rows"] for stage in certificates), "Combined manifest test cardinality")
    require(manifest["policy_configurations"] == len(read_json(complete / "policies.json")), "Combined manifest policy cardinality")
    require(manifest["reporting_source_sha256"] == digest(ROOT / "experiments/reporting/p02_complete.py"), "Reporting source changed since output generation")
    sources = read_json(complete / "source_records.json")
    require(sources == manifest["source_records"], "Source-reference JSON differs from combined manifest")
    require(set(sources) == {"stage1", "stage2"}, "Combined source reference contains wrong stages")
    source_files = 0
    for stage, files in sources.items():
        for name, record in files.items():
            expected_path = (run / stage / name).relative_to(ROOT).as_posix()
            require(record["path"] == expected_path and digest(ROOT / record["path"]) == record["sha256"], f"Changed combined source reference: {stage}/{name}")
            source_files += 1
    for name, expected in manifest["output_hashes"].items():
        require(digest(complete / name) == expected, f"Changed reporting output: {name}")
    values = _digests(manifest)
    required = {"summary.json": digest(complete / "summary.json"),
                "policies.json": digest(complete / "policies.json")}
    for stage, records in certificates.items():
        for field in ("raw_sha256", "summary_sha256", "independent_audit_sha256"):
            required[f"{stage}/{field}"] = records[field]
    required.update({name: record["sha256"] for name, record in csvs.items()})
    for name, expected in required.items():
        require(expected in values, f"Combined manifest omits or changes critical digest: {name}")
    return {"manifest_sha256": digest(complete / "manifest.json"),
            "root_frozen_manifest_sha256": digest(run / "manifest.json"),
            "critical_source_and_output_hashes_checked": len(required),
            "all_referenced_source_files": source_files, "all_hashed_reporting_outputs": len(manifest["output_hashes"])}


def audit_reporting_tables(run, complete):
    """Check the compact tables used in the manuscript against raw/summary."""
    summary = read_json(complete / "summary.json")
    risk, pairs = [], []
    metric_names = ("dangerous_H1_release_given_screen", "H0_withholding_given_screen", "H0_retention_overall")
    for stratum in summary["strata"]:
        budget = str(float(summary["primary_budget"]))
        pool = stratum["pools"]["confirmation_all"]["budgets"][budget]
        for policy, record in pool["policies"].items():
            row = {"n": stratum["n"], "rel": stratum["rel"], "policy": policy, "feasible_rounds": record["feasible_rounds"]}
            for name in metric_names:
                stats = record["metrics"][name]
                interval = stats.get("crossed")
                row[name] = stats["point_mean"]
                row[f"{name}_ci_low"] = interval["ci95"][0] if interval and interval["ci95"] else None
                row[f"{name}_ci_high"] = interval["ci95"][1] if interval and interval["ci95"] else None
                row[f"{name}_defined_draws"] = interval["defined_draws"] if interval else 0
            row["danger_contributing_groups"] = json.dumps(record["metrics"][metric_names[0]]["denominator_contributing_test_groups_by_round"])
            risk.append(row)
        for pair, stats in pool["paired"].items():
            value = stats[metric_names[0]]
            interval = value.get("crossed")
            pairs.append({"n": stratum["n"], "rel": stratum["rel"], "pair": pair,
                          "difference": value["point_mean"],
                          "ci_low": interval["ci95"][0] if interval and interval["ci95"] else None,
                          "ci_high": interval["ci95"][1] if interval and interval["ci95"] else None,
                          "defined_draws": interval["defined_draws"] if interval else 0})
    coverage = []
    columns = ["phase", "n", "rel", "klass", "split", "source_kind", "group_unit", "r2", "weighted_mean_relative_error", "weighted_p95_relative_error"]
    for stage in ("stage1", "stage2"):
        raw = pd.read_parquet(run / stage / "raw.parquet", columns=columns)
        for (n, rel), frame in raw[raw.phase.eq("test")].groupby(["n", "rel"], sort=True):
            h0 = frame[frame.klass.eq("H0") & frame.split.eq("confirmation")]
            h1 = frame[frame.klass.eq("H1") & frame.split.eq("confirmation")]
            physical = h1[h1.source_kind.eq("new_source_positive_restriction")]
            h0_screen = np.isfinite(h0.r2) & h0.r2.ge(.99)
            h1_danger = np.isfinite(h1.r2) & h1.r2.ge(.99) & np.isfinite(h1.weighted_mean_relative_error) & h1.weighted_mean_relative_error.gt(.1)
            physical_screen = np.isfinite(physical.r2) & physical.r2.ge(.99)
            physical_danger = physical_screen & np.isfinite(physical.weighted_mean_relative_error) & physical.weighted_mean_relative_error.gt(.1)
            physical_p95 = physical_screen & np.isfinite(physical.weighted_p95_relative_error) & physical.weighted_p95_relative_error.gt(.1)
            coverage.append({"n": n, "rel": rel, "confirmation_h0_total": len(h0), "confirmation_h0_screened": int(h0_screen.sum()),
                             "confirmation_h0_screen_rate": float(h0_screen.mean()),
                             "confirmation_dangerous_screened": int(h1_danger.sum()),
                             "confirmation_danger_structure_groups": int(h1.loc[h1_danger, "group_unit"].nunique()),
                             "physical_total": len(physical), "physical_screened": int(physical_screen.sum()),
                             "physical_dangerous_screened": int(physical_danger.sum()),
                             "physical_danger_structure_groups": int(physical.loc[physical_danger, "group_unit"].nunique()),
                             "physical_p95_dangerous_screened": int(physical_p95.sum()),
                             "physical_p95_danger_structure_groups": int(physical.loc[physical_p95, "group_unit"].nunique())})
    checked = {}
    for name, rows, keys in [("primary_confirmation.csv", risk, ["n", "rel", "policy"]),
                              ("paired_confirmation.csv", pairs, ["n", "rel", "pair"]),
                              ("coverage.csv", coverage, ["n", "rel"])]:
        actual = pd.read_csv(complete / name, float_precision="round_trip").sort_values(keys).reset_index(drop=True)
        expected = pd.DataFrame(rows).sort_values(keys).reset_index(drop=True)
        require(set(actual.columns) == set(expected.columns), f"Compact reporting table has wrong columns: {name}")
        pd.testing.assert_frame_equal(actual[expected.columns], expected, check_dtype=False,
                                      check_exact=False, rtol=2e-15, atol=0, obj=name)
        checked[name] = len(actual)
    return checked


def audit_physical_selection(directory):
    """Review metadata-only anti-selection rules; not an execution certificate."""
    manifest = read_json(directory / "review_manifest.json")
    for name, expected in manifest["source_and_input_sha256"].items():
        require(digest(Path(name)) == expected, f"Physical registration input changed: {name}")
        if "/runs/" in Path(name).as_posix():
            require(Path(name).as_posix().endswith(("/runs/p02_inventory/inventory.json", "/runs/p02_fixed_domain/task_registry.json")), f"Unexpected performance-bearing registration input: {name}")
    for name, expected in manifest["metadata_artifact_sha256"].items():
        require(digest(directory / name) == expected, f"Physical metadata artifact changed: {name}")
    protocol = read_json(directory / "protocol_rules.json")
    require(protocol == read_json(ROOT / "experiments/physical_coverage_v1/protocol.json"), "Physical rule snapshot differs from source")
    require(protocol["registration_status"] == "pre_execution_registration_incomplete" and manifest["not_an_execution_freeze"], "Physical metadata must not be labelled an execution freeze")
    require(manifest["forbidden_performance_inputs_read"] == [], "Forbidden performance inputs reported in physical registration")
    records = read_json(directory / "candidates.json")["equations"]
    inventory = read_json(ROOT / "runs/p02_inventory/inventory.json")["equations"]
    source = {record["eq_name"]: record for record in inventory}
    require(len(records) == 120 and {r["eq_name"] for r in records} == set(source), "Physical census must retain all120 unique source equations")
    valid_status = set(protocol["scientific_decisions"])
    full_domains = 0
    type_pending, documented_relaxations = [], 0
    for record in records:
        name = record["eq_name"]
        original = source[name]
        for field in ("expression", "sympy_srepr", "sampling_objects", "n_vars"):
            require(record[field] == original[field], f"{name}: source formula/coordinates changed")
        require(record["scientific_status"] in valid_status, f"{name}: unknown scientific status")
        require(record["analytic_math_status"] == record["scientific_status"], f"{name}: mathematical status changed by contract readiness")
        require(record["physical_source_contract_review_status"].startswith("pending_"), f"{name}: pending physical semantics represented as certified")
        require("independent_physical_source_contract_review_pending" in record["execution_gates"], f"{name}: missing physical-contract execution block")
        expected_types = []
        for match in re.finditer(r"- x\[(\d+)\]:[^\n]*?\(([^)]*)\)", original["upstream_docstring"] or ""):
            if "integer" not in match.group(2).lower():
                continue
            index = int(match.group(1))
            if index >= record["n_vars"]:
                expected_types.append((None, "source_type_review_pending"))
                continue
            sampler = original["sampling_objects"][index]["class_name"]
            status = ("source_integer_sampling_consistent" if sampler == "IntegerSampling" else
                      "source_documented_integer_to_real_relaxation" if "real" in match.group(2).lower() else
                      "source_type_review_pending")
            expected_types.append((index, status))
            documented_relaxations += status == "source_documented_integer_to_real_relaxation"
        saved_types = [(item.get("variable_index"), item["status"]) for item in record["source_variable_type_review"]]
        require(saved_types == expected_types, f"{name}: integer descriptions not reviewed consistently across source catalog")
        needs_type_review = any(status == "source_type_review_pending" for _, status in expected_types)
        require(("source_type_review_pending" in record["execution_gates"]) == needs_type_review, f"{name}: unresolved type contract not gated")
        if needs_type_review:
            type_pending.append(name)
        require(record["inclusion"] == (record["scientific_status"] == "eligible_analytic"), f"{name}: inconsistent admission")
        require(isinstance(record["positive_laurent_sum_tag"]["matches"], bool), f"{name}: missing all-catalog fine grammar tag")
        require(isinstance(record["nested_additive_radical_tag"], bool), f"{name}: missing all-catalog coarse structure tag")
        expected, scales, support_error = [], [], None
        for sampling in original["sampling_objects"]:
            a = sampling["attributes"]
            lo, hi, kind = a["min_value"], a["max_value"], sampling["class_name"]
            positive, negative = a.get("uses_positive", True), a.get("uses_negative", True)
            if kind == "IntegerSampling":
                raw = set(range(int(lo), int(hi)))
                values = (raw if positive else set()) | ({-v for v in raw} if negative else set())
                values = sorted(v for v in values if v > 0)
                if len(values) < 2:
                    support_error = "degenerate_positive_integer_support"
                    break
                bounds = [values[max(0, math.ceil(q*len(values))-1)] for q in (.05,.95)]
                scale = "integer"
            elif kind == "DefaultSampling":
                if not positive:
                    support_error = "no_original_positive_support"
                    break
                bounds = [math.exp(math.log(lo)+q*(math.log(hi)-math.log(lo))) for q in (.05,.95)]
                scale = "log"
            else:
                require(kind == "SimpleSampling", f"{name}: unknown original sampler")
                intervals = []
                if positive and hi > 0:
                    intervals.append((max(float(lo),0.), float(hi)))
                if negative and lo < 0:
                    intervals.append((max(float(-hi),0.), float(-lo)))
                if not intervals:
                    support_error = "no_original_positive_support"
                    break
                require(len(set(intervals)) == 1, f"{name}: unequal positive mixture needs explicit CDF review")
                a0, b0 = intervals[0]
                bounds = [a0+q*(b0-a0) for q in (.05,.95)]
                scale = "linear"
            expected.append(bounds)
            scales.append(scale)
        if support_error:
            require(not record["inclusion"] and record["outer_bounds"] is None, f"{name}: inadmissible original support silently repaired")
        else:
            require(record["coordinate_scales"] == scales and np.allclose(record["outer_bounds"], expected, rtol=3e-14, atol=0), f"{name}: domain is not the uniform original-support quantile box")
            require(set(record["train_bounds_by_width"]) == {str(w) for w in protocol["common_factorial_design"]["training_width_fractions"]}, f"{name}: missing training width")
            for width in protocol["common_factorial_design"]["training_width_fractions"]:
                train = []
                for (lo, hi), scale in zip(expected, scales, strict=True):
                    a0, b0 = (math.log(lo),math.log(hi)) if scale == "log" else (lo,hi)
                    ends = [(a0+b0)/2-width*(b0-a0)/2, (a0+b0)/2+width*(b0-a0)/2]
                    if scale == "log":
                        ends = [math.exp(value) for value in ends]
                    elif scale == "integer":
                        ends = [math.ceil(ends[0]), math.floor(ends[1])]
                    train.append(ends)
                require(np.allclose(record["train_bounds_by_width"][str(width)], train, rtol=3e-14, atol=0), f"{name}: source-specific training window rescue")
            full_domains += 1
        if record["inclusion"]:
            require(record["proof_obligations"] == [] and all(c["certified"] for c in record["source_constraint_checks"]), f"{name}: admitted with unresolved proof obligations")
            certificate = record["truth_interval_certificate"]
            require(mp.mpf(certificate["lower"]) > 0 and mp.isfinite(mp.mpf(certificate["upper"])), f"{name}: admitted without positive finite interval certificate")
            require(certificate["random_or_grid_point_evaluations"] == 0, f"{name}: sample-dependent domain admission")
            require(record["execution_status"] == "await_adapter_and_registration", f"{name}: incomplete physical task advertised runnable")
    groups = read_json(directory / "groups.json")["groups"]
    members = [name for group in groups for name in group["members"]]
    require(len(members) == 120 and set(members) == set(source), "Primary physical groups do not partition all120 sources")
    by_name = {r["eq_name"]:r for r in records}
    for group in groups:
        for name in group["members"]:
            record = by_name[name]
            require(record["structural_group"] == group["structural_group"], f"{name}: primary group mismatch")
            if group["exposure"] == "development":
                require(record["exposure"] in {"development", "development_conservative"}, f"{name}: prior group exposure not inherited")
    ledger = read_json(directory / "exposure_ledger.json")
    old_tasks = read_json(ROOT / "runs/p02_fixed_domain/task_registry.json")
    direct = {t["task"] for t in old_tasks if t["family"] == "upstream_source"}
    require(set(ledger["P02_physical_sources"]) == direct and len(direct) == 38, "Physical exposure ledger changes prior38 source identities")
    for name in direct:
        require(by_name[name]["source_exposure_level"] == "prior_study_outcomes" and by_name[name]["exposure"].startswith("development"), f"{name}: prior source advertised new")
    for item in ledger["P02_designed_mechanisms"]:
        if item["registered_role"] == "reserve":
            require(item["exposure_level"] == "mechanism_design_only", "Reserved mechanism misrepresented as evaluated")
    summary = read_json(directory / "summary.json")
    require(summary["scientific_status"] == dict(Counter(r["scientific_status"] for r in records)), "Physical status summary mismatch")
    prospective = [r for r in records if r["inclusion"] and r["klass"] == "H1" and r["exposure"] == "prospective_candidate"]
    require(summary["eligible_H1_prospective_source_equations"] == len(prospective) == 0, "Primary undisputed prospective count must remain zero")
    disputed = [r for r in records if r["inclusion"] and r["klass"] == "H1" and r["structural_assignment_status"] == "candidate_pending_structural_review"]
    require(summary["eligible_H1_disputed_fine_structure_candidates"] == len(disputed) == 1, "Disputed fine candidate count mismatch")
    require(all(summary[name] == 0 for name in ("study_samples_generated", "study_fits_generated", "diagnostic_outcomes_generated")), "Physical metadata includes study outcomes")
    require(summary["source_type_review_pending_equations"] == len(type_pending), "Type readiness summary mismatch")
    metadata_names = {"protocol_rules.json", "candidates.json", "groups.json", "summary.json", "exposure_ledger.json", "candidates.csv", "REGISTRATION.md"}
    require(set(manifest["metadata_artifact_sha256"]) == metadata_names, "Metadata hash inputs include audit certificates or archived snapshots")
    source_path = ROOT / "experiments/physical_coverage_v1/registration.py"
    source_tree = ast.parse(source_path.read_text(encoding="utf-8"))
    functions = {node.name: node for node in source_tree.body if isinstance(node, ast.FunctionDef)}
    writer = functions["write_registration"]
    exists_branches = [node for node in writer.body if isinstance(node, ast.If) and any(
        isinstance(call, ast.Call) and isinstance(call.func, ast.Attribute) and call.func.attr == "exists"
        for call in ast.walk(node.test))]
    require(bool(exists_branches), "Registration writer lacks an existing-snapshot guard")
    guarded = exists_branches[0]
    require(any(isinstance(node, ast.Return) for node in guarded.body), "Existing registration does not return without rebuilding")
    require(any(isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "verify_registration" for node in ast.walk(guarded)), "Existing registration does not verify before returning")
    forbidden_writes = {"_write_json", "write_text", "write_bytes", "mkdir", "copy2", "copyfile", "replace", "unlink"}
    for node in ast.walk(guarded):
        if isinstance(node, ast.Call):
            name = node.func.id if isinstance(node.func, ast.Name) else node.func.attr if isinstance(node.func, ast.Attribute) else ""
            require(name not in forbidden_writes, "Existing-registration guard contains a mutation")
    verifier = functions["verify_registration"]
    require(any(isinstance(node, ast.Raise) for node in ast.walk(verifier)), "Registration verifier does not reject drift")
    tests_path = source_path.with_name("test_registration.py")
    if not tests_path.exists():
        tests_path = ROOT / "tests/test_physical_registration.py"
    test_tree = ast.parse(tests_path.read_text(encoding="utf-8"))
    test_names = [node.name for node in test_tree.body if isinstance(node, ast.FunctionDef) and node.name.startswith("test_")]
    require(len(test_names) == 3 and any("read_only" in name and "drift" in name for name in test_names), "Missing saved-registration idempotence/drift tests")
    return {"passed": True, "scope": "Independent anti-selection, source-fidelity, domain-rule and exposure-metadata review only",
            "rules_sha256": manifest["rules_sha256"], "equations": 120, "complete_original_rule_domains": full_domains,
            "scientific_status_counts": summary["scientific_status"], "primary_groups": len(groups),
            "eligible_uncontested_prospective_H1_groups": 0,
            "disputed_fine_candidates": [r["eq_name"] for r in disputed],
            "source_type_pending_equations": type_pending, "documented_integer_to_real_relaxations": documented_relaxations,
            "physical_source_contract_certified_equations": 0,
            "overwrite_guard": {"review": "Independent static AST/source/test review; no study implementation imported or executed",
                                "same_snapshot": "verify then return without writes", "drift": "reject",
                                "tests_present": test_names, "compiler_sha256": digest(source_path)},
            "not_an_execution_freeze": True,
            "not_verified_here": ["Independent interval-arithmetic reevaluation of every source constraint and target bound", "Independent physical/source/type/unit semantics certification", "Float64/log-truth/high-precision adapters", "High-dimensional diagnostics and frozen grids", "Future null controls, seed namespaces and execution roles", "Absence of unlogged human information exposure cannot be established by saved metadata alone"],
            "selection_evidence": "All120 source rows retained; expressions/constants/original coordinates unchanged; common positive-support quantiles and centered widths independently reconstructed; admission rules have no fit or diagnostic outcomes; primary/fine grouping uncertainty explicitly retained."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", default="runs/p02_fixed_domain")
    parser.add_argument("--physical-review-only", action="store_true")
    args = parser.parse_args()
    if args.physical_review_only:
        directory = ROOT / "runs/p03_physical_protocol"
        result = audit_physical_selection(directory)
        result["audit_source_sha256"] = digest(Path(__file__))
        (directory / "independent_selection_review.json").write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")
        print(json.dumps(result, ensure_ascii=False))
        return
    run = ROOT / args.run
    complete = run / "complete_grid"
    protocol = read_json(run / "manifest.json")["protocol"]
    result = {"audit": "independent_p02_complete_grid", "audit_source_sha256": digest(Path(__file__)),
              "imports_experimental_implementation": False, "preservation": audit_preservation(run),
              "stage_certificates": audit_stage_certificates(run)}
    result["raw_disjointness"] = audit_disjoint_raw(run, protocol)
    result["summary"] = audit_summary(run, complete, result["raw_disjointness"])
    result["csv_concatenation"] = audit_csvs(run, complete)
    result["policies"] = audit_policies(run, complete)
    result["manuscript_tables"] = audit_reporting_tables(run, complete)
    result["manifest"] = audit_manifest(run, complete, result["stage_certificates"], result["csv_concatenation"])
    result["boundaries"] = ["No A-G single-gate statistic is recomputed from sampled observations.",
                             "Per-stage independent audits reconstruct policy calibration, fixed-grid errors, group event counts, primary/equal-group and paired uncertainty.",
                             "The complete-grid audit preserves separate stratum estimates; no meta-analysis or new independent experiment is implied."]
    result["passed"] = True
    output = complete / "independent_audit.json"
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({"passed": True, "output": str(output), "raw_rows": result["raw_disjointness"]["raw_rows"],
                      "strata": result["summary"]["exact_source_strata"], "policies": result["policies"]["policy_configurations"]}))


if __name__ == "__main__":
    main()
