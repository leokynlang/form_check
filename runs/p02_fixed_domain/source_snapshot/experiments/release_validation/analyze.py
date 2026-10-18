"""H0-only policy fitting followed by one frozen independent-sample evaluation."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from experiments.release_validation.calibration import (
    apply_policy,
    calibrate,
    fit_reference,
)
from experiments.release_validation.inventory import ROOT, digest
from experiments.release_validation.runner import assert_dataset_isolation, write_json

RATIOS = {
    "H0_withholding_given_screen": ("h0_withheld_screened", "h0_screened"),
    "H0_alarm_given_screen": ("h0_alarm_screened", "h0_screened"),
    "H0_decline_given_screen": ("h0_decline_screened", "h0_screened"),
    "H0_retention_overall": ("h0_released", "h0_total"),
    "H1_alarm_yield": ("h1_alarm", "h1_total"),
    "H1_decline_rate": ("h1_decline", "h1_total"),
    "dangerous_H1_release_given_screen": ("h1_danger_released", "h1_danger_screened"),
    "dangerous_H1_release_overall": ("h1_danger_released", "h1_danger_total"),
    "safe_H1_retention_given_screen": ("h1_safe_released", "h1_safe_screened"),
    "release_overall": ("released", "total"),
    "dangerous_H1_share_of_known_outcome_releases": ("h1_danger_released", "released_known_outcome"),
    "all_dangerous_share_of_known_outcome_releases": ("all_danger_released", "released_known_outcome"),
    "outcome_unavailable_rate": ("outcome_unavailable", "total"),
}
PAIRS = [("AD_panel", "D"), ("AD_panel", "F"), ("AEF_panel", "F"),
         ("G_rbf_cv", "F"), ("D", "F")]


def per_task_counts(frame):
    h0, h1 = frame.klass.eq("H0"), frame.klass.eq("H1")
    screen, alarm, decline, release = [frame[key] for key in ("screen", "alarm", "decline", "released")]
    known = np.isfinite(frame.legacy_endpoint_error)
    danger = known & frame.legacy_endpoint_error.gt(frame.attrs["danger_threshold"])
    safe = known & ~danger
    events = {
        "total": pd.Series(True, index=frame.index), "released": release,
        "h0_total": h0, "h0_screened": h0 & screen, "h0_released": h0 & release,
        "h0_withheld_screened": h0 & screen & (alarm | decline),
        "h0_alarm_screened": h0 & screen & alarm,
        "h0_decline_screened": h0 & screen & decline,
        "h1_total": h1, "h1_alarm": h1 & alarm, "h1_decline": h1 & decline,
        "h1_danger_total": h1 & danger, "h1_danger_screened": h1 & danger & screen,
        "h1_danger_released": h1 & danger & release,
        "h1_safe_screened": h1 & safe & screen, "h1_safe_released": h1 & safe & release,
        "all_danger_released": danger & release,
        "h0_danger_released": h0 & danger & release,
        "released_known_outcome": release & known,
        "outcome_unavailable": ~known,
        "released_unknown_outcome": release & ~known,
        "data_failures": frame.data_failure.notna(),
        "fit_failures": frame.fit_failure.notna(),
    }
    table = pd.DataFrame({key: value.astype(int) for key, value in events.items()})
    return table.groupby(frame.task, sort=True).sum()


def ratio(a, b):
    return float(a / b) if b else None


def bootstrap_summary(counts, draws):
    totals = counts.sum()
    result = {"counts": {k: int(v) for k, v in totals.items()}, "metrics": {}}
    for name, (numerator, denominator) in RATIOS.items():
        num = counts[numerator].to_numpy()
        den = counts[denominator].to_numpy()
        bn, bd = num[draws].sum(axis=1), den[draws].sum(axis=1)
        ok = bd > 0
        ci = np.quantile(bn[ok] / bd[ok], [.025, .975]).tolist() if ok.any() else None
        result["metrics"][name] = {
            "value": ratio(totals[numerator], totals[denominator]),
            "numerator": int(totals[numerator]), "denominator": int(totals[denominator]),
            "task_bootstrap_ci95": ci, "bootstrap_defined_draws": int(ok.sum()),
            "denominator_contributing_tasks": int(np.sum(den > 0)),
        }
    return result


def paired_effect(a, b, draws):
    effects = {}
    if not a.index.equals(b.index):
        raise ValueError("Paired task populations differ")
    for metric in ("dangerous_H1_release_given_screen", "H0_retention_overall"):
        numerator, denominator = RATIOS[metric]
        an, ad = a[numerator].to_numpy(), a[denominator].to_numpy()
        bn, bd = b[numerator].to_numpy(), b[denominator].to_numpy()
        if not np.array_equal(ad, bd):
            raise ValueError("Paired denominators differ")
        den = ad[draws].sum(axis=1)
        ok = den > 0
        diff = (an[draws].sum(axis=1)[ok] - bn[draws].sum(axis=1)[ok]) / den[ok]
        effects[metric] = {
            "difference_first_minus_second": ratio(an.sum() - bn.sum(), ad.sum()),
            "paired_task_bootstrap_ci95": np.quantile(diff, [.025, .975]).tolist() if len(diff) else None,
            "denominator_contributing_tasks": int(np.sum(ad > 0)),
            "bootstrap_defined_draws": int(ok.sum()),
        }
    return effects


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", default="runs/p01_pilot_v1")
    args = parser.parse_args()
    run = (ROOT / args.run).resolve()
    if run.parent != (ROOT / "runs").resolve() or not run.name.startswith("p01_"):
        raise ValueError("Analysis writes only to its dedicated runs/p01_* directory")
    manifest = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    if not manifest["complete"]:
        raise ValueError("Run is incomplete")
    path = run / "raw_policy.parquet"
    if digest(path) != manifest["raw_sha256"]:
        raise ValueError("Raw record hash mismatch")
    frame = pd.read_parquet(path)
    assert_dataset_isolation(frame)
    protocol = manifest["protocol"]
    frame["screen"] = np.isfinite(frame.r2) & frame.r2.ge(protocol["screen_r2"])
    summary = {"protocol_id": protocol["protocol_id"], "scope": protocol["scope"],
               "raw_sha256": manifest["raw_sha256"], "analysis_sha256": digest(Path(__file__)),
               "inference": protocol["inference"], "strata": [],
               "risk_mixture": "Descriptive counts in a fixed 13-H0/19-H1 task mixture; unknown outcomes excluded from risk denominator and reported separately."}
    fitted_policies, decisions, task_tables = [], [], []
    for (n, rel), stratum in frame.groupby(["n", "rel"], sort=True):
        reference = fit_reference(stratum[(stratum.phase == "reference") & stratum.screen], protocol["methods"])
        cal = stratum[(stratum.phase == "calibration") & stratum.screen]
        test = stratum[stratum.phase == "test"].copy()
        tasks = sorted(test.task.unique())
        draws = np.random.default_rng(protocol["master_seed"] + int(n)).integers(
            0, len(tasks), size=(protocol["bootstrap_repetitions"], len(tasks)))
        stratum_result = {"n": int(n), "rel": float(rel), "reference_screened":
                          int(((stratum.phase == "reference") & stratum.screen).sum()),
                          "calibration_screened": len(cal), "test_rows": len(test),
                          "budgets": {}}
        for budget in protocol["budgets"]:
            budget_result, counts_by_policy = {}, {}
            for policy, members in protocol["policies"].items():
                fitted = calibrate(cal, members, reference, budget)
                # A finite empirical CDF has an attainable ceiling. A threshold
                # at that ceiling disables all valid alarms; disclose this
                # calibration limitation without changing the frozen v1 rule.
                ceilings = [len(reference[m]) / (len(reference[m]) + 1) for m in members]
                ceiling = max(ceilings)
                fitted["score_ceiling"] = ceiling
                fitted["all_valid_alarms_disabled_by_ceiling"] = bool(
                    fitted["feasible"] and fitted["threshold"] >= ceiling)
                fitted_policies.append({"n": int(n), "rel": float(rel), "budget": budget,
                                        "policy": policy, "members": members, "fitted": fitted,
                                        "reference": {method: reference[method] for method in members}})
                if not fitted["feasible"]:
                    budget_result[policy] = {"feasible": False, "calibration": fitted}
                    continue
                score, alarm, decline = apply_policy(test, members, reference, fitted)
                evaluated = test.copy()
                evaluated["policy_score"], evaluated["alarm"], evaluated["decline"] = score, alarm, decline
                evaluated["released"] = evaluated.screen & ~alarm & ~decline
                evaluated.attrs["danger_threshold"] = protocol["danger_threshold"]
                counts = per_task_counts(evaluated).reindex(tasks, fill_value=0)
                counts_by_policy[policy] = counts
                budget_result[policy] = {"feasible": True, "calibration": fitted,
                                         **bootstrap_summary(counts, draws)}
                keep = ["task", "klass", "n", "rel", "replicate", "dataset_sha256", "screen",
                        "policy_score", "alarm", "decline", "released", "legacy_endpoint_error"]
                decision = evaluated[keep].copy()
                decision["budget"], decision["policy"] = budget, policy
                decisions.append(decision)
                task_counts = counts.reset_index()
                task_counts["n"], task_counts["rel"] = n, rel
                task_counts["budget"], task_counts["policy"] = budget, policy
                task_tables.append(task_counts)
            paired = {}
            for first, second in PAIRS:
                if first in counts_by_policy and second in counts_by_policy:
                    paired[f"{first}_minus_{second}"] = paired_effect(
                        counts_by_policy[first], counts_by_policy[second], draws)
            stratum_result["budgets"][str(budget)] = {"policies": budget_result, "paired": paired}
        summary["strata"].append(stratum_result)
    write_json(run / "policies.json", fitted_policies)
    write_json(run / "summary.json", summary)
    pd.concat(decisions, ignore_index=True).to_parquet(run / "decisions.parquet", index=False)
    pd.concat(task_tables, ignore_index=True).to_csv(run / "per_task_counts.csv", index=False)
    failures = []
    for method in protocol["methods"]:
        failed = frame[frame[f"fail_{method}"].notna()]
        for (phase, task, reason), group in failed.groupby(["phase", "task", f"fail_{method}"]):
            failures.append({"method": method, "phase": phase, "task": task,
                             "reason": reason, "count": len(group)})
    write_json(run / "failures.json", failures)
    write_report(run, summary)
    plot_tradeoffs(run, summary)
    print(json.dumps({"output": str(run / "summary.json"), "strata": len(summary["strata"]),
                      "status": "development pilot; not unseen-equation validation"}))


def fmt(value):
    return "undefined" if value is None else f"{100 * value:.2f}%"


def write_report(run, summary):
    lines = ["# P0-1 independent-sample development pilot", "", summary["scope"], "",
             "Reference and calibration use H0 only. Every method sees identical data.",
             "Budget = rejection OR decline among screened H0; calibration and test are separate.",
             "Test size is observed, not guaranteed equal to its target. Legacy endpoint only.",
             "Task-bootstrap intervals freeze calibration; they omit calibration uncertainty.", "",
             "| n | noise | target | policy | calibration withheld | test H0 withheld | H0 retained overall | dangerous H1 released / screened dangerous | unavailable outcomes |",
             "|---|---|---|---|---|---|---|---|---|"]
    for stratum in summary["strata"]:
        for budget, result in stratum["budgets"].items():
            for policy, record in result["policies"].items():
                prefix = f"| {stratum['n']} | {stratum['rel']} | {fmt(float(budget))} | {policy} |"
                if not record["feasible"]:
                    lines.append(prefix + " infeasible | - | - | - | - |")
                    continue
                metrics = record["metrics"]
                dangerous = metrics["dangerous_H1_release_given_screen"]
                lines.append(prefix + f" {fmt(record['calibration']['withheld_rate'])} | "
                             + fmt(metrics["H0_withholding_given_screen"]["value"]) + " | "
                             + fmt(metrics["H0_retention_overall"]["value"]) + " | "
                             + f"{dangerous['numerator']}/{dangerous['denominator']} ({fmt(dangerous['value'])}) | "
                             + f"{record['counts']['outcome_unavailable']} |")
    lines += ["", "## Interpretation limits", "",
              "- The equations were previously studied. Fresh seeds give sample independence, not untouched-equation confirmation.",
              "- A budget-compliant calibration fit can exceed its budget on independent test data; report that gap without retuning.",
              "- Endpoint-unavailable records remain in operational metrics, with risk unlabelled.",
              "- Fixed external domains, structural-group holdouts and mechanism confirmation remain pending.",
              "- Finite reference CDFs can saturate: at n=50 and 1% target the AEF threshold reaches the score ceiling, disabling all valid alarms. This diagnoses the v1 score mapping, not inherent panel power. Do not tune v1 from its test outcomes; revise the protocol and use new samples.",
              "- No winner is declared from this development pilot. Per-task counts and paired intervals are in summary.json."]
    (run / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot_tradeoffs(run, summary):
    fig, axes = plt.subplots(1, len(summary["strata"]), figsize=(6 * len(summary["strata"]), 4.8), squeeze=False)
    for ax, stratum in zip(axes[0], summary["strata"], strict=True):
        policies = next(iter(stratum["budgets"].values()))["policies"]
        for policy in policies:
            points = []
            for _budget, result in stratum["budgets"].items():
                record = result["policies"][policy]
                if record["feasible"]:
                    x = record["metrics"]["H0_retention_overall"]["value"]
                    y = record["metrics"]["dangerous_H1_release_given_screen"]["value"]
                    if x is not None and y is not None:
                        points.append((x, y))
            if points:
                ax.plot([100*x for x, _ in points], [100*y for _, y in points], "o-", label=policy)
        ax.set(title=f"n={stratum['n']}, noise={stratum['rel']}", xlabel="Correct-form H0 retained overall (%)",
               ylabel="Screened dangerous H1 released (%)")
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
    fig.suptitle("Development pilot / fresh samples, known equations / legacy endpoint", fontsize=11)
    fig.tight_layout()
    fig.savefig(run / "tradeoffs.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
