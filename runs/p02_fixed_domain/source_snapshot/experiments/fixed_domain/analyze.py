"""Frozen P0-2 policy evaluation with crossed structure/design-group uncertainty."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from experiments.release_validation.analyze import PAIRS, RATIOS, per_task_counts
from experiments.release_validation.analyze_v2 import (
    aggregate_policy,
    bootstrap_design,
    divide,
    mean_uncertainty,
    paired_aggregate,
)
from experiments.release_validation.calibration_v2 import (
    apply_policy,
    calibrate,
    fit_reference,
)
from experiments.release_validation.inventory import ROOT, digest
from experiments.release_validation.runner import write_json

from .grouping import group_unit

PRIMARY_OUTCOME = "weighted_mean_relative_error"
RISK_METRICS = (
    "dangerous_H1_release_given_screen",
    "dangerous_H1_release_overall",
    "all_dangerous_share_of_known_outcome_releases",
    "outcome_unavailable_rate",
)
POOLS = (
    "confirmation_all", "confirmation_physical_sources",
    "confirmation_designed_controls", "development", "stress",
)


def _verify_frozen(run):
    # Import at execution time so helper tests do not need a complete run.
    from .freeze import verify_frozen
    return verify_frozen(run)


def pool_mask(frame, pool):
    """Named pools always retain the independent confirmation H0 controls."""
    h0 = frame.klass.eq("H0") & frame.split.eq("confirmation")
    confirm_h1 = frame.klass.eq("H1") & frame.split.eq("confirmation")
    if pool == "confirmation_all":
        return h0 | confirm_h1
    if pool == "confirmation_physical_sources":
        return h0 | (confirm_h1 & frame.source_kind.eq("new_source_positive_restriction"))
    if pool == "confirmation_designed_controls":
        return h0 | (confirm_h1 & frame.source_kind.eq("synthetic_positive_mechanism"))
    if pool == "development":
        return frame.split.eq("development")
    if pool == "stress":
        return h0 | (frame.klass.eq("H1") & frame.split.eq("stress"))
    raise ValueError(f"Unknown registered evaluation pool: {pool}")


def outcome_specs(protocol, columns):
    specs = [(PRIMARY_OUTCOME, float(t)) for t in protocol["danger_thresholds"]]
    primary = float(protocol["primary_danger_threshold"])
    specs += [(name, primary) for name in protocol["outcome"]["secondary"] if name in columns]
    if (PRIMARY_OUTCOME, primary) not in specs:
        raise ValueError("Primary outcome and danger threshold not registered")
    return list(dict.fromkeys(specs))


def outcome_task_counts(evaluated, outcome, threshold):
    """Reuse audited event accounting while replacing only the truth outcome."""
    frame = evaluated.copy()
    frame["legacy_endpoint_error"] = frame[outcome]
    frame.attrs["danger_threshold"] = threshold
    return per_task_counts(frame)


def group_counts(task_counts, task_metadata, tasks, groups):
    """Subset tasks BEFORE grouping, because a group may span source kinds."""
    selected = task_counts.reindex(tasks, fill_value=0)
    labels = task_metadata.loc[selected.index, "group_unit"]
    labels.index = selected.index
    result = selected.groupby(labels, sort=True).sum().reindex(groups, fill_value=0)
    result.index.name = "group_unit"
    return result


def _rename_group_fields(value):
    if isinstance(value, list):
        return [_rename_group_fields(v) for v in value]
    if not isinstance(value, dict):
        return value
    renamed = {}
    for key, item in value.items():
        output_key = {"conditional_task": "conditional_group",
               "denominator_contributing_test_tasks_by_round":
                   "denominator_contributing_test_groups_by_round"}.get(key, key)
        renamed[output_key] = _rename_group_fields(item)
    return renamed


def equal_group_risk(cube, columns, weights, round_draws, feasible):
    """Secondary average over defined group-specific dangerous release ratios.

    H0 and H1 groups with no screened-dangerous H1 records have zero denominator
    and are explicitly omitted. Bootstrap draws with no eligible group remain
    undefined. No undefined group or bootstrap ratio is replaced with zero.
    """
    numerator, denominator = RATIOS["dangerous_H1_release_given_screen"]
    num = cube[:, :, columns.index(numerator)]
    den = cube[:, :, columns.index(denominator)]
    ratios = divide(num, den)
    defined = np.isfinite(ratios)
    per_round = divide(np.nansum(ratios, axis=1), defined.sum(axis=1))
    weighted = divide(weights @ np.nan_to_num(ratios, nan=0).T,
                      weights @ defined.astype(float).T)
    result = _rename_group_fields(mean_uncertainty(per_round, weighted, round_draws,
                                                  bool(np.all(feasible))))
    result["defined_groups_by_round"] = defined.sum(axis=1).astype(int).tolist()
    result["omitted_groups_by_round"] = (~defined).sum(axis=1).astype(int).tolist()
    result["omitted_group_rule"] = "No defined screened-dangerous H1 denominator; includes all H0 groups"
    result["estimand"] = "Equal mean of defined group-specific screened-dangerous H1 release ratios, then equal round mean"
    return result


def _aggregate(cube, columns, weights, round_draws, feasible, budget):
    result = _rename_group_fields(aggregate_policy(cube, columns, weights, round_draws, feasible, budget))
    result["equal_group_dangerous_release"] = equal_group_risk(
        cube, columns, weights, round_draws, feasible)
    return result


def _sensitivity(cube, columns, weights, round_draws, feasible):
    # Keep all outcomes in group_counts, but avoid repeating invariant H0 policy
    # accounting in each sensitivity summary.
    result = {"all_rounds_feasible": bool(np.all(feasible)), "metrics": {}}
    for metric in RISK_METRICS:
        numerator, denominator = RATIOS[metric]
        num = cube[:, :, columns.index(numerator)]
        den = cube[:, :, columns.index(denominator)]
        original = divide(num.sum(axis=1), den.sum(axis=1))
        weighted = divide(weights @ num.T, weights @ den.T)
        stats = _rename_group_fields(mean_uncertainty(original, weighted, round_draws,
                                                       bool(np.all(feasible))))
        stats["denominator_contributing_test_groups_by_round"] = (den > 0).sum(axis=1).astype(int).tolist()
        result["metrics"][metric] = stats
    result["equal_group_dangerous_release"] = equal_group_risk(
        cube, columns, weights, round_draws, feasible)
    return result


def _interval(stats, key="crossed"):
    return stats[key]["ci95"] if stats.get(key) else None


def _flat_stats(stats):
    interval = _interval(stats)
    crossed = stats.get("crossed") or {}
    sensitivity = stats["calibration_sensitivity"]
    return {"point_mean": stats["point_mean"], "primary_defined": stats["primary_defined"],
            "calibration_round_sd": sensitivity["sample_sd"],
            "calibration_round_min": sensitivity["min"],
            "calibration_round_max": sensitivity["max"],
            "crossed_ci_low": interval[0] if interval else None,
            "crossed_ci_high": interval[1] if interval else None,
            "crossed_defined_draws": crossed.get("defined_draws"),
            "crossed_total_draws": crossed.get("total_draws")}


def _registry_metadata(run, manifest, raw):
    path = run / "task_registry.json"
    expected = manifest.get("task_registry_sha256", manifest.get("registry_sha256"))
    if expected is None or digest(path) != expected:
        raise ValueError("Frozen task registry hash mismatch")
    registry = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(registry, dict):
        registry = registry.get("tasks", registry.get("registry"))
    if not isinstance(registry, list):
        raise ValueError("Frozen registry must contain a task list")
    records = []
    for task in registry:
        records.append({"task": task["task"], "klass": task["klass"],
                        "split": task["split"], "source_kind": task["source_kind"],
                        "group_unit": group_unit(task)})
    metadata = pd.DataFrame(records).set_index("task")
    if not metadata.index.is_unique:
        raise ValueError("Repeated task identity in frozen registry")
    if not set(raw.task).issubset(set(metadata.index)):
        raise ValueError("Raw record contains an unregistered task")
    frame = raw.copy()
    for column in metadata:
        expected_values = frame.task.map(metadata[column])
        if column in frame and not frame[column].eq(expected_values).all():
            raise ValueError(f"Raw task metadata changed: {column}")
        frame[column] = expected_values
    return frame, metadata


def _pool_design(test, metadata, pool, rounds, repetitions, seed):
    selected = test.loc[pool_mask(test, pool)]
    tasks = sorted(selected.task.unique())
    groups = sorted(selected.group_unit.unique())
    group_classes = selected.groupby("group_unit").klass.agg(lambda x: tuple(sorted(set(x))))
    if any(len(x) != 1 for x in group_classes):
        raise ValueError("One bootstrap group contains both H0 and H1")
    classes = np.array([group_classes.loc[g][0] for g in groups])
    if not groups:
        raise ValueError(f"Registered pool {pool} is empty")
    weights, round_draws = bootstrap_design(classes, rounds, repetitions, seed)
    info = {"test_rows": len(selected), "test_tasks": len(tasks), "test_groups": len(groups),
            "H0_tasks": int(metadata.loc[tasks].klass.eq("H0").sum()),
            "H1_tasks": int(metadata.loc[tasks].klass.eq("H1").sum()),
            "H0_design_groups": int(np.sum(classes == "H0")),
            "H1_structure_groups": int(np.sum(classes == "H1")),
            "group_units": groups,
            "task_ids": tasks,
            "group_task_counts": metadata.loc[tasks].groupby("group_unit").size().astype(int).to_dict()}
    return {"tasks": tasks, "groups": groups, "weights": weights,
            "round_draws": round_draws, "info": info}


def analyze_run(run, stage="stage1"):
    run = Path(run).resolve()
    _verify_frozen(run)
    manifest = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    protocol = manifest["protocol"]
    if stage not in protocol["stages"]:
        raise ValueError("Stage is not registered")
    out = run / stage
    execution = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    if not execution["complete"] or digest(out / "raw.parquet") != execution["raw_sha256"]:
        raise ValueError("Incomplete stage or changed raw data")
    raw = pd.read_parquet(out / "raw.parquet")
    if len(raw) != execution["raw_rows"]:
        raise ValueError("Raw row count differs from completed stage manifest")
    raw, metadata = _registry_metadata(run, manifest, raw)
    keys = ["phase", "calibration_round", "task", "n", "rel", "replicate"]
    if raw.duplicated(keys).any():
        raise ValueError("Repeated dataset record key")
    non_test = raw.phase.ne("test")
    if (not raw.loc[non_test, "klass"].eq("H0").all()
            or not raw.loc[non_test, "split"].eq("calibration").all()
            or not raw.phase.isin(["reference", "calibration", "test"]).all()):
        raise ValueError("Reference/calibration leakage or unknown phase")
    if raw.split.eq("reserve").any():
        raise ValueError("Reserved tasks must not have generated records")
    raw["screen"] = np.isfinite(raw.r2) & raw.r2.ge(protocol["screen_r2"])
    rounds = int(protocol["calibration_rounds"])
    expected_rounds = set(range(rounds))
    specs = outcome_specs(protocol, raw.columns)
    primary_spec = (PRIMARY_OUTCOME, float(protocol["primary_danger_threshold"]))
    decisions, policy_records, references, group_tables = [], [], [], []
    sensitivity_rows, primary_rows, round_rows = [], [], []
    summary = {
        "protocol_id": protocol["protocol_id"], "stage": stage,
        "raw_sha256": execution["raw_sha256"], "raw_rows": len(raw),
        "analysis_sha256": digest(Path(__file__)),
        "primary_outcome": PRIMARY_OUTCOME, "primary_danger_threshold": primary_spec[1],
        "primary_budget": float(protocol["primary_budget"]),
        "calibration_rounds": rounds, "group_inference": protocol["group_inference"],
        "point_estimator": "Equal mean of round-specific pooled-record ratios within each named pool",
        "scope": "Registered structure/design-group library with fixed pre-sample external grids. Old tasks are development; H0 transfer concerns held-out design groups within the monomial family. Intervals condition on the fixed calibration-design pool and shared realized test bank.",
        "uncertainty_names": {
            "conditional_group": "Class-stratified structure/design-group bootstrap with fitted calibration rounds fixed",
            "round_mean": "Bootstrap of the empirical independent reference/calibration rounds with test groups fixed",
            "crossed": "One shared group-weight vector across all selected rounds and policies, crossed with one shared calibration-round vector",
            "round_percentile_range": "Empirical calibration sensitivity range, not a confidence interval"},
        "sensitivity_table": "sensitivity_metrics.csv",
        "strata": [],
    }
    for (n, rel), stratum in raw.groupby(["n", "rel"], sort=True):
        test = stratum.loc[stratum.phase.eq("test")].copy()
        if not test.calibration_round.eq(-1).all():
            raise ValueError("Shared test bank must have calibration_round=-1")
        pool_designs = {}
        for pool in protocol["evaluation_pools"]:
            pool_seed = int.from_bytes(hashlib.sha256(pool.encode()).digest()[:4], "little")
            pool_designs[pool] = _pool_design(
                test, metadata, pool, rounds,
                int(protocol["group_inference"]["bootstrap_repetitions"]),
                int(protocol["master_seed"]) + int(n) + round(float(rel) * 1e6) + pool_seed)
        fitted_references, calibration_frames = {}, {}
        for phase in ("reference", "calibration"):
            if set(stratum.loc[stratum.phase.eq(phase), "calibration_round"].unique()) != expected_rounds:
                raise ValueError(f"Missing or extra {phase} calibration round")
        for calibration_round in range(rounds):
            ref = stratum.loc[stratum.phase.eq("reference") & stratum.calibration_round.eq(calibration_round) & stratum.screen]
            cal = stratum.loc[stratum.phase.eq("calibration") & stratum.calibration_round.eq(calibration_round) & stratum.screen]
            fitted_references[calibration_round] = fit_reference(ref, protocol["methods"])
            calibration_frames[calibration_round] = cal
            for method, reference in fitted_references[calibration_round].items():
                references.append({"n": int(n), "rel": float(rel), "calibration_round": calibration_round,
                                   "method": method, **reference})
        result_stratum = {"n": int(n), "rel": float(rel), "test_rows": len(test),
                          "pools": {pool: {**design["info"], "budgets": {}}
                                    for pool, design in pool_designs.items()}}
        for budget in protocol["budgets"]:
            policy_cubes, feasible_maps = {}, {}
            for policy, members in protocol["policies"].items():
                lists = {(pool, spec): [] for pool in pool_designs for spec in specs}
                feasible = []
                # Registered event column schema, even if no policy is feasible.
                empty = test.copy()
                empty["alarm"], empty["decline"], empty["released"] = False, False, False
                columns = outcome_task_counts(empty, *primary_spec).columns.tolist()
                for calibration_round in range(rounds):
                    reference = fitted_references[calibration_round]
                    fitted = calibrate(calibration_frames[calibration_round], members, reference, float(budget))
                    feasible.append(bool(fitted["feasible"]))
                    policy_records.append({"n": int(n), "rel": float(rel), "budget": float(budget),
                                           "policy": policy, "members": members,
                                           "calibration_round": calibration_round,
                                           "reference": {m: reference[m] for m in members}, "fitted": fitted})
                    if fitted["feasible"]:
                        score, alarm, decline = apply_policy(test, members, reference, fitted)
                        evaluated = test.copy()
                        evaluated["policy_score"] = score
                        evaluated["alarm"], evaluated["decline"] = alarm, decline
                        evaluated["released"] = evaluated.screen & ~alarm & ~decline
                        keep = ["task", "group_unit", "klass", "split", "source_kind", "n", "rel", "replicate",
                                "dataset_sha256", "screen", "policy_score", "alarm", "decline", "released"]
                        keep += list(dict.fromkeys(spec[0] for spec in specs))
                        decision = evaluated[keep].copy()
                        decision["calibration_round"] = calibration_round
                        decision["budget"], decision["policy"] = float(budget), policy
                        decisions.append(decision)
                        task_tables = {spec: outcome_task_counts(evaluated, *spec) for spec in specs}
                    for pool, design in pool_designs.items():
                        for spec in specs:
                            if fitted["feasible"]:
                                counts = group_counts(task_tables[spec], metadata,
                                                      design["tasks"], design["groups"])
                                values = counts.to_numpy(dtype=float)
                            else:
                                values = np.full((len(design["groups"]), len(columns)), np.nan)
                                counts = pd.DataFrame(values, index=pd.Index(design["groups"], name="group_unit"), columns=columns)
                            lists[(pool, spec)].append(values)
                            table = counts.reset_index()
                            table["pool"], table["outcome"] = pool, spec[0]
                            table["danger_threshold"] = spec[1]
                            table["n"], table["rel"] = int(n), float(rel)
                            table["calibration_round"] = calibration_round
                            table["budget"], table["policy"] = float(budget), policy
                            table["feasible"] = bool(fitted["feasible"])
                            group_tables.append(table)
                            if spec == primary_spec:
                                total = counts.sum(min_count=1)
                                item = {"pool": pool, "n": int(n), "rel": float(rel), "budget": float(budget),
                                        "policy": policy, "calibration_round": calibration_round,
                                        "feasible": bool(fitted["feasible"]),
                                        "calibration_withheld": fitted.get("withheld_rate"),
                                        "threshold": fitted.get("threshold")}
                                for metric, (num, den) in RATIOS.items():
                                    item[metric] = float(total[num] / total[den]) if np.isfinite(total[den]) and total[den] > 0 else None
                                round_rows.append(item)
                feasible_maps[policy] = np.asarray(feasible, dtype=bool)
                policy_cubes[policy] = {key: np.stack(values) for key, values in lists.items()}
            for pool, design in pool_designs.items():
                weights, round_draws = design["weights"], design["round_draws"]
                budget_result = {"policies": {}, "paired": {}}
                for policy in protocol["policies"]:
                    cube = policy_cubes[policy][(pool, primary_spec)]
                    result = _aggregate(cube, columns, weights, round_draws, feasible_maps[policy], float(budget))
                    budget_result["policies"][policy] = result
                    for metric, stats in result["metrics"].items():
                        primary_rows.append({"pool": pool, "n": int(n), "rel": float(rel),
                                             "budget": float(budget), "policy": policy,
                                             "metric": metric, **_flat_stats(stats)})
                    for spec in specs:
                        if spec == primary_spec:
                            continue
                        secondary = _sensitivity(policy_cubes[policy][(pool, spec)], columns, weights,
                                                 round_draws, feasible_maps[policy])
                        for metric, stats in secondary["metrics"].items():
                            sensitivity_rows.append({"pool": pool, "n": int(n), "rel": float(rel),
                                                     "budget": float(budget), "policy": policy,
                                                     "outcome": spec[0], "danger_threshold": spec[1],
                                                     "metric": metric, **_flat_stats(stats),
                                                     "contributing_groups_by_round": json.dumps(stats["denominator_contributing_test_groups_by_round"])})
                        sensitivity_rows.append({"pool": pool, "n": int(n), "rel": float(rel),
                                                 "budget": float(budget), "policy": policy,
                                                 "outcome": spec[0], "danger_threshold": spec[1],
                                                 "metric": "equal_group_dangerous_release",
                                                 **_flat_stats(secondary["equal_group_dangerous_release"]),
                                                 "contributing_groups_by_round": json.dumps(secondary["equal_group_dangerous_release"]["defined_groups_by_round"])})
                for first, second in PAIRS:
                    if first not in policy_cubes or second not in policy_cubes:
                        continue
                    both = feasible_maps[first] & feasible_maps[second]
                    budget_result["paired"][f"{first}_minus_{second}"] = _rename_group_fields(paired_aggregate(
                        policy_cubes[first][(pool, primary_spec)], policy_cubes[second][(pool, primary_spec)],
                        columns, weights, round_draws, both))
                result_stratum["pools"][pool]["budgets"][str(float(budget))] = budget_result
        summary["strata"].append(result_stratum)
    write_json(out / "summary.json", summary)
    write_json(out / "policies.json", policy_records)
    pd.DataFrame(references).to_csv(out / "reference_scales.csv", index=False)
    pd.DataFrame(primary_rows).to_csv(out / "primary_metrics.csv", index=False)
    pd.DataFrame(sensitivity_rows).to_csv(out / "sensitivity_metrics.csv", index=False)
    pd.DataFrame(round_rows).to_csv(out / "round_metrics.csv", index=False)
    pd.concat(group_tables, ignore_index=True).to_parquet(out / "group_counts.parquet", index=False)
    if decisions:
        pd.concat(decisions, ignore_index=True).to_parquet(out / "decisions.parquet", index=False)
    else:
        pd.DataFrame(columns=["task", "group_unit", "policy", "budget", "calibration_round"]).to_parquet(out / "decisions.parquet", index=False)
    report(out, summary)
    _verify_frozen(run)
    return summary


def _pct(value):
    return "undefined" if value is None else f"{100 * value:.2f}%"


def _interval_text(stats):
    values = _interval(stats)
    return "undefined" if values is None else f"[{_pct(values[0])}, {_pct(values[1])}]"


def report(out, summary):
    rounds = summary["calibration_rounds"]
    lines = ["# P0-2: fixed external regions and held-out structural groups", "", summary["scope"], "",
             f"Stage: {summary['stage']}. Independent reference/calibration rounds: {rounds}.",
             "Primary danger: fixed weighted mean pointwise relative error >10%; fit screen: training log-space R2>=0.99.",
             "The table uses the registered 5% joint withholding budget. Actual independent-test H0 withholding is measured; transfer is not guaranteed.",
             "Physical-source and designed-control subpools overlap in H0 controls and may overlap in structural groups. They are not independent replications.", "",
             "| n | pool | H0/H1 tasks | H0/H1 groups | policy | feasible rounds | H0 withheld | screened-dangerous H1 released | crossed group/round interval | contributing H1 groups | defined crossed draws | equal-group release |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    budget_key = str(float(summary["primary_budget"]))
    for stratum in summary["strata"]:
        for pool, details in stratum["pools"].items():
            for policy, record in details["budgets"][budget_key]["policies"].items():
                h0 = record["metrics"]["H0_withholding_given_screen"]
                risk = record["metrics"]["dangerous_H1_release_given_screen"]
                groups = sorted(set(risk["denominator_contributing_test_groups_by_round"]))
                defined = risk.get("crossed") or {}
                equal = record["equal_group_dangerous_release"]
                lines.append(f"| {stratum['n']} | {pool} | {details['H0_tasks']}/{details['H1_tasks']} | "
                             f"{details['H0_design_groups']}/{details['H1_structure_groups']} | {policy} | "
                             f"{record['feasible_rounds']}/{rounds} | {_pct(h0['point_mean'])} | {_pct(risk['point_mean'])} | "
                             f"{_interval_text(risk)} | {','.join(map(str, groups))} | "
                             f"{defined.get('defined_draws', 0)}/{defined.get('total_draws', 0)} | {_pct(equal['point_mean'])} |")
    lines += ["", "## Estimands and uncertainty", "",
              "- Primary points are equal means of round-specific pooled-record ratios. Reusing one test bank across rounds does not multiply its independent sample size.",
              "- H1 resampling units are structural groups; H0 units are frozen design groups. Parameter windows and near-twin formulas remain inside their groups.",
              "- Conditional-group, round-only, and crossed intervals use paired draws across every policy within a pool. Crossed group weights stay fixed across selected rounds.",
              "- Undefined bootstrap ratios remain undefined. summary.json records defined/total draws, contributing group counts, feasibility, and calibration sensitivity ranges.",
              "- The equal-group risk is secondary: it averages only groups with a defined screened-dangerous H1 denominator. Omitted group counts are explicit; H0 groups are excluded from this risk average.",
              f"- {rounds}-round empirical percentile ranges describe calibration sensitivity, not confidence intervals or coverage guarantees.",
              "- All intervals condition on this fixed calibration-design pool and realized test bank. No record-level Monte Carlo bootstrap or universal physical-family guarantee is claimed.",
              "- All targets include method declines. Infeasible rounds are not silently discarded from primary aggregates.",
              "- Truth/unavailable outcomes remain in operational counts; unknown errors are excluded only from explicitly known-outcome risk denominators.", "",
              "## Additional tables", "",
              "- primary_metrics.csv contains all three budgets and all registered operational metrics.",
              "- sensitivity_metrics.csv contains 2%/5%/20% mean-error cutoffs and registered weighted-p95/max/legacy outcomes at the primary cutoff.",
              "- group_counts.parquet retains every endpoint, threshold, pool, policy and calibration round for independent reconstruction.",
              "- round_metrics.csv, policies.json and reference_scales.csv expose threshold transfer and repeated-calibration variation.",
              "- decisions.parquet stores one decision per shared test record/policy/budget/round, not duplicated by evaluation pool.",
              "- Reserved mechanism groups have no generated samples, fits or diagnostic results. V1 and V2 remain separate preserved runs."]
    (out / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", default="runs/p02_fixed_domain")
    parser.add_argument("--stage", default="stage1")
    args = parser.parse_args()
    run = (ROOT / args.run).resolve()
    if run.parent != (ROOT / "runs").resolve() or not run.name.startswith("p02_fixed_domain"):
        raise ValueError("Analysis writes only to a dedicated P0-2 fixed-domain run")
    summary = analyze_run(run, args.stage)
    print(json.dumps({"run": str(run), "stage": args.stage, "raw_rows": summary["raw_rows"],
                      "strata": len(summary["strata"]), "frozen_and_preserved": True}))


if __name__ == "__main__":
    main()
