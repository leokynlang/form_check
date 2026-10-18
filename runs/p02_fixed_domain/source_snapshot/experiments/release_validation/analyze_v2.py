"""Repeated-calibration v2 estimates and paired crossed empirical uncertainty."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from experiments.release_validation.analyze import PAIRS, RATIOS, per_task_counts
from experiments.release_validation.calibration_v2 import (
    apply_policy,
    calibrate,
    fit_reference,
)
from experiments.release_validation.inventory import ROOT, digest
from experiments.release_validation.runner import write_json
from experiments.release_validation.runner_v2 import V1, tree_hashes, verify_new_records


def bootstrap_design(task_classes, rounds, repetitions, seed):
    """Share one stratified task vector across rounds/policies in each draw."""
    rng = np.random.default_rng(seed)
    weights = np.zeros((repetitions, len(task_classes)), dtype=int)
    for klass in ("H0", "H1"):
        positions = np.flatnonzero(np.asarray(task_classes) == klass)
        sampled = rng.choice(positions, size=(repetitions, len(positions)), replace=True)
        for task in positions:
            weights[:, task] = np.sum(sampled == task, axis=1)
    return weights, rng.integers(0, rounds, size=(repetitions, rounds))


def divide(num, den):
    result = np.full(np.broadcast_shapes(np.shape(num), np.shape(den)), np.nan)
    return np.divide(num, den, out=result, where=np.asarray(den) > 0)


def empirical_interval(values):
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    return {"ci95": np.quantile(values[finite], [.025, .975], method="linear").tolist() if finite.any() else None,
            "defined_draws": int(finite.sum()), "total_draws": len(values)}


def ratio_arrays(cube, columns, numerator, denominator, weights):
    num, den = cube[:, :, columns.index(numerator)], cube[:, :, columns.index(denominator)]
    original = divide(num.sum(axis=1), den.sum(axis=1))
    weighted = divide(weights @ num.T, weights @ den.T)
    return original, weighted


def mean_uncertainty(original, weighted, round_draws, all_feasible):
    finite = np.isfinite(original)
    primary_defined = bool(all_feasible and finite.all())
    finite_values = original[finite]
    sensitivity = {
        "defined_rounds": int(finite.sum()),
        "sample_sd": float(np.std(finite_values, ddof=1)) if len(finite_values) > 1 else None,
        "min": float(np.min(finite_values)) if len(finite_values) else None,
        "max": float(np.max(finite_values)) if len(finite_values) else None,
        "empirical_percentile_range95": np.quantile(finite_values, [.025, .975], method="linear").tolist() if len(finite_values) else None,
        "interpretation": "Calibration-round sensitivity range conditional on fixed realized test bank; NOT a confidence interval",
    }
    result = {"point_mean": float(np.mean(original)) if primary_defined else None,
              "primary_defined": primary_defined,
              "feasible_defined_round_only_sensitivity_mean": float(np.mean(finite_values)) if len(finite_values) else None,
              "calibration_sensitivity": sensitivity}
    if primary_defined:
        conditional = np.mean(weighted, axis=1)
        round_only = np.mean(original[round_draws], axis=1)
        crossed = np.mean(np.take_along_axis(weighted, round_draws, axis=1), axis=1)
        result.update({"conditional_task": empirical_interval(conditional),
                       "round_mean": empirical_interval(round_only),
                       "crossed": empirical_interval(crossed)})
    else:
        result.update({"conditional_task": None, "round_mean": None, "crossed": None})
    return result


def aggregate_policy(cube, columns, weights, round_draws, feasible, budget):
    all_feasible = bool(np.all(feasible))
    result = {"feasible_rounds": int(np.sum(feasible)), "total_rounds": len(feasible),
              "all_rounds_feasible": all_feasible, "metrics": {}}
    for metric, (numerator, denominator) in RATIOS.items():
        original, weighted = ratio_arrays(cube, columns, numerator, denominator, weights)
        stats = mean_uncertainty(original, weighted, round_draws, all_feasible)
        den = cube[:, :, columns.index(denominator)]
        stats["denominator_contributing_test_tasks_by_round"] = np.sum(den > 0, axis=1).astype(int).tolist()
        if metric == "H0_withholding_given_screen":
            stats["test_rounds_above_nominal"] = int(np.sum(np.isfinite(original) & (original > budget)))
            stats["test_fraction_rounds_above_nominal"] = float(np.mean(original > budget)) if all_feasible and np.isfinite(original).all() else None
        result["metrics"][metric] = stats
    return result


def paired_aggregate(first, second, columns, weights, round_draws, feasible):
    result = {}
    for metric in ("dangerous_H1_release_given_screen", "H0_retention_overall", "H0_withholding_given_screen"):
        numerator, denominator = RATIOS[metric]
        den_idx = columns.index(denominator)
        if not np.array_equal(first[feasible, :, den_idx], second[feasible, :, den_idx], equal_nan=True):
            raise ValueError("Paired per-task denominators differ")
        a, aw = ratio_arrays(first, columns, numerator, denominator, weights)
        b, bw = ratio_arrays(second, columns, numerator, denominator, weights)
        result[metric] = mean_uncertainty(a - b, aw - bw, round_draws, bool(np.all(feasible)))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", default="runs/p01_v2")
    args = parser.parse_args()
    run = (ROOT / args.run).resolve()
    if run.parent != (ROOT / "runs").resolve() or not run.name.startswith("p01_v2"):
        raise ValueError("Analysis writes only into a dedicated v2 directory")
    manifest = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    if not manifest["complete"] or digest(run / "raw_policy.parquet") != manifest["raw_sha256"]:
        raise ValueError("Incomplete run or changed raw records")
    if tree_hashes(V1) != manifest["v1_recursive_snapshot"]:
        raise ValueError("v1 preservation failed")
    protocol = manifest["protocol"]
    raw = pd.read_parquet(run / "raw_policy.parquet")
    verify_new_records(raw, protocol)
    raw["screen"] = np.isfinite(raw.r2) & raw.r2.ge(protocol["screen_r2"])
    rounds = protocol["calibration_rounds"]
    decisions, policies, task_tables, round_metrics, references = [], [], [], [], []
    summary = {"protocol_id": protocol["protocol_id"], "raw_sha256": manifest["raw_sha256"],
               "analysis_sha256": digest(Path(__file__)), "uncertainty": protocol["uncertainty"],
               "point_estimator": protocol["point_estimator"], "strata": [],
               "scope": "Known 32 equations; fully fresh data and20 independent reference/calibration pairs; shared independent test bank; legacy endpoint outcome."}
    for (n, rel), frame in raw.groupby(["n", "rel"], sort=True):
        test = frame[frame.phase == "test"].copy()
        tasks = sorted(test.task.unique())
        task_classes = test.groupby("task").klass.first().reindex(tasks).to_numpy()
        weights, round_draws = bootstrap_design(task_classes, rounds,
                                                protocol["uncertainty"]["bootstrap_repetitions"],
                                                protocol["master_seed"] + int(n) + 9001)
        fitted_refs, calibration_frames = {}, {}
        for calibration_round in range(rounds):
            reference_rows = frame[(frame.phase == "reference") & frame.calibration_round.eq(calibration_round) & frame.screen]
            fitted_refs[calibration_round] = fit_reference(reference_rows, protocol["methods"])
            calibration_frames[calibration_round] = frame[(frame.phase == "calibration") & frame.calibration_round.eq(calibration_round) & frame.screen]
            for method, reference in fitted_refs[calibration_round].items():
                references.append({"n": int(n), "rel": float(rel), "calibration_round": calibration_round,
                                   "method": method, **reference})
        stratum = {"n": int(n), "rel": float(rel), "test_rows": len(test),
                   "test_H0_tasks": int(np.sum(task_classes == "H0")), "test_H1_tasks": int(np.sum(task_classes == "H1")),
                   "budgets": {}}
        for budget in protocol["budgets"]:
            cubes, feasible_maps, columns = {}, {}, None
            for policy, members in protocol["policies"].items():
                round_counts, feasible = [], []
                for calibration_round in range(rounds):
                    reference = fitted_refs[calibration_round]
                    cal = calibration_frames[calibration_round]
                    fitted = calibrate(cal, members, reference, budget)
                    policies.append({"n": int(n), "rel": float(rel), "calibration_round": calibration_round,
                                     "budget": budget, "policy": policy, "members": members,
                                     "reference": {m: reference[m] for m in members}, "fitted": fitted})
                    feasible.append(fitted["feasible"])
                    if not fitted["feasible"]:
                        round_counts.append(None)
                        round_metrics.append({"n": n, "rel": rel, "calibration_round": calibration_round,
                                              "budget": budget, "policy": policy, "feasible": False})
                        continue
                    score, alarm, decline = apply_policy(test, members, reference, fitted)
                    evaluated = test.copy()
                    evaluated["policy_score"], evaluated["alarm"], evaluated["decline"] = score, alarm, decline
                    evaluated["released"] = evaluated.screen & ~alarm & ~decline
                    evaluated.attrs["danger_threshold"] = protocol["danger_threshold"]
                    counts = per_task_counts(evaluated).reindex(tasks, fill_value=0)
                    columns = counts.columns.tolist()
                    round_counts.append(counts.to_numpy(dtype=float))
                    total = counts.sum()
                    row = {"n": n, "rel": rel, "calibration_round": calibration_round,
                           "budget": budget, "policy": policy, "feasible": True,
                           "calibration_withheld": fitted["withheld_rate"], "threshold": fitted["threshold"],
                           **{key: int(value) for key, value in total.items()}}
                    for metric, (num, den) in RATIOS.items():
                        row[metric] = float(total[num] / total[den]) if total[den] else None
                    round_metrics.append(row)
                    keep = ["task", "klass", "n", "rel", "replicate", "dataset_sha256", "screen",
                            "policy_score", "alarm", "decline", "released", "legacy_endpoint_error"]
                    decision = evaluated[keep].copy()
                    decision["calibration_round"] = calibration_round
                    decision["budget"], decision["policy"] = budget, policy
                    decisions.append(decision)
                    counts = counts.reset_index()
                    counts["n"], counts["rel"] = n, rel
                    counts["calibration_round"] = calibration_round
                    counts["budget"], counts["policy"] = budget, policy
                    task_tables.append(counts)
                if columns is None:
                    # Fixed event schema even if every policy was infeasible.
                    empty = test.copy()
                    for key in ("alarm", "decline", "released"):
                        empty[key] = False
                    empty.attrs["danger_threshold"] = protocol["danger_threshold"]
                    columns = per_task_counts(empty).columns.tolist()
                cube = np.stack([value if value is not None else np.full((len(tasks), len(columns)), np.nan)
                                 for value in round_counts])
                cubes[policy], feasible_maps[policy] = cube, np.array(feasible)
            budget_result = {"policies": {policy: aggregate_policy(cubes[policy], columns, weights, round_draws,
                                                                    feasible_maps[policy], budget)
                                          for policy in protocol["policies"]}, "paired": {}}
            for first, second in PAIRS:
                budget_result["paired"][f"{first}_minus_{second}"] = paired_aggregate(
                    cubes[first], cubes[second], columns, weights, round_draws,
                    feasible_maps[first] & feasible_maps[second])
            stratum["budgets"][str(budget)] = budget_result
        summary["strata"].append(stratum)
    write_json(run / "policies.json", policies)
    write_json(run / "summary.json", summary)
    pd.DataFrame(references).to_csv(run / "reference_scales.csv", index=False)
    pd.DataFrame(round_metrics).to_csv(run / "round_metrics.csv", index=False)
    pd.concat(task_tables, ignore_index=True).to_csv(run / "per_task_counts.csv", index=False)
    pd.concat(decisions, ignore_index=True).to_parquet(run / "decisions.parquet", index=False)
    score_diagnostics(run, raw, protocol)
    report(run, summary, pd.DataFrame(round_metrics))
    plots(run, summary, pd.DataFrame(round_metrics))
    if tree_hashes(V1) != manifest["v1_recursive_snapshot"]:
        raise AssertionError("v1 changed during analysis")
    print(json.dumps({"run": str(run), "rounds": rounds, "raw_rows": len(raw),
                      "decision_rows": sum(len(d) for d in decisions), "v1_preserved": True}))


def score_diagnostics(run, raw, protocol):
    diagnostics = []
    for (phase, n), frame in raw.groupby(["phase", "n"]):
        for method in protocol["methods"]:
            score = frame[f"score_{method}"]
            diagnostics.append({"phase": phase, "n": int(n), "method": method, "N": len(frame),
                                "raw_nonfinite": int((~np.isfinite(score)).sum()),
                                "typed_declines": int(frame[f"fail_{method}"].notna().sum()),
                                "p_zero": int(score.eq(0).sum()) if method in "ABCD" else 0,
                                "A_permutation_floor": int(score.eq(-.001).sum()) if method == "A" else 0,
                                "float_p_clamped": int((score.le(0) & score.ge(-np.finfo(float).tiny)).sum()) if method in "ABCD" else 0})
    write_json(run / "score_diagnostics.json", diagnostics)


def pct(value):
    return "undefined" if value is None else f"{100*value:.2f}%"


def report(run, summary, round_frame):
    lines = ["# P0-1 v2: fresh data and repeated independent calibration", "", summary["scope"], "",
             "20 independent H0 reference/calibration pairs; mean of round-specific ratios.",
             "All targets include declines. Independent test error budgets are measured, not guaranteed.",
             "The affine score has no reference-CDF ceiling; intrinsic p-value floors/ties remain.", "",
             "| n | target | policy | feasible rounds | mean test H0 withheld | round sensitivity range | mean screened-dangerous H1 release | crossed interval |",
             "|---|---|---|---|---|---|---|---|"]
    for stratum in summary["strata"]:
        for budget, record in stratum["budgets"].items():
            for policy, result in record["policies"].items():
                h0 = result["metrics"]["H0_withholding_given_screen"]
                danger = result["metrics"]["dangerous_H1_release_given_screen"]
                spread = h0["calibration_sensitivity"]["empirical_percentile_range95"]
                ci = danger["crossed"]["ci95"] if danger["crossed"] else None
                interval = lambda values: "undefined" if values is None else f"[{pct(values[0])}, {pct(values[1])}]"
                lines.append(f"| {stratum['n']} | {pct(float(budget))} | {policy} | {result['feasible_rounds']}/{result['total_rounds']} | "
                             f"{pct(h0['point_mean'])} | {interval(spread)} | {pct(danger['point_mean'])} | {interval(ci)} |")
    lines += ["", "## Uncertainty interpretation", "",
              "- Calibration SD/min/max/percentile range measures sensitivity to independently redrawn reference/calibration samples on the same13 H0 equations; the range is not a CI.",
              "- Conditional task intervals keep all20 fitted policies fixed and resample test equations stratified13 H0/19 H1.",
              "- Round-only intervals resample the20 empirical calibration rounds with test tasks fixed.",
              "- Crossed intervals jointly resample rounds and class-stratified test-task weights; the same vectors are paired across policies.",
              "- All intervals condition on the original H0 calibration task identities and shared realized test bank. No new-equation guarantee or complete record-level Monte Carlo interval is claimed.",
              "- Primary aggregates require every round feasible and its ratio defined. Feasible-only estimates are explicitly secondary.",
              "- Legacy endpoint and known tasks remain; fixed domains and unseen structural-group confirmation are pending.",
              "- V1 remains byte-for-byte unchanged; v2 uses no v1 sample seeds or dataset hashes.", "",
              "## Calibration implementation diagnostics", "",
              f"AEF at n=50, 1% target has nonzero test alarms in {int(((round_frame.policy == 'AEF_panel') & (round_frame.n == 50) & (round_frame.budget == .01) & round_frame.h1_alarm.gt(0)).sum())}/20 rounds.",
              "This removes the v1 CDF-ceiling pathology, without claiming that all remaining p-value ties disappear."]
    (run / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def plots(run, summary, rounds):
    fig, axes = plt.subplots(1, len(summary["strata"]), figsize=(6*len(summary["strata"]), 4.8), squeeze=False)
    for ax, stratum in zip(axes[0], summary["strata"], strict=True):
        names = next(iter(stratum["budgets"].values()))["policies"]
        for name in names:
            xs, ys = [], []
            for record in stratum["budgets"].values():
                result = record["policies"][name]["metrics"]
                x, y = result["H0_retention_overall"]["point_mean"], result["dangerous_H1_release_given_screen"]["point_mean"]
                if x is not None and y is not None:
                    xs.append(100*x)
                    ys.append(100*y)
            if xs:
                ax.plot(xs, ys, "o-", label=name)
        ax.set(title=f"n={stratum['n']}; means over20 calibrations", xlabel="Correct H0 retained overall (%)",
               ylabel="Screened dangerous H1 released (%)")
        ax.legend(fontsize=8)
        ax.grid(alpha=.2)
    fig.suptitle("V2 / fresh samples, known equations, legacy endpoint", fontsize=11)
    fig.tight_layout()
    fig.savefig(run / "tradeoffs.png", dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(1, len(summary["strata"]), figsize=(6*len(summary["strata"]), 4.8), squeeze=False)
    for ax, stratum in zip(axes[0], summary["strata"], strict=True):
        frame = rounds[(rounds.n == stratum["n"]) & (rounds.budget == .05)]
        names = sorted(frame.policy.unique())
        arrays = [100*frame.loc[frame.policy == name, "H0_withholding_given_screen"].dropna().to_numpy() for name in names]
        ax.boxplot(arrays, tick_labels=names)
        ax.axhline(5, color="black", linestyle="--", linewidth=1)
        ax.set(title=f"n={stratum['n']};20 independent calibrations", ylabel="Test H0 withholding (%)")
        ax.tick_params(axis="x", rotation=20)
        ax.grid(axis="y", alpha=.2)
    fig.suptitle("Calibration sensitivity at5% target / conditional on shared test bank", fontsize=11)
    fig.tight_layout()
    fig.savefig(run / "calibration_sensitivity.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
