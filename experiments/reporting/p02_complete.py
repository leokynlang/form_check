"""Combine completed P0-2 strata without pooling estimands or changing policies."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from experiments.fixed_domain.freeze import verify_frozen
from experiments.release_validation.inventory import ROOT, digest
from experiments.release_validation.runner import write_json
from experiments.release_validation.runner_v2 import tree_hashes

POLICIES = {"D": "D", "F": "F", "AD_panel": "A-D", "AEF_panel": "A/E/F", "G_rbf_cv": "RBF CV"}
COLORS = {"D": "#777777", "F": "#007f6d", "AD_panel": "#8b71a9", "AEF_panel": "#c05431", "G_rbf_cv": "#457ba0"}
CSVS = ["primary_metrics.csv", "sensitivity_metrics.csv", "round_metrics.csv", "reference_scales.csv"]


def verify_preservation(run):
    saved = json.loads((run / "stage2_preservation_before.json").read_text(encoding="utf-8"))
    if tree_hashes(run / "stage1") != saved["stage1_files"]:
        raise ValueError("Stage1 changed during stage2")
    if tree_hashes(ROOT / "paper/research_R2") != saved["paper_R2_files"]:
        raise ValueError("Original R2 changed during R3 development")
    return {"stage1_files": len(saved["stage1_files"]), "paper_R2_files": len(saved["paper_R2_files"])}


def combine_summaries(inputs):
    summaries = [value for _, value in inputs]
    first = summaries[0]
    invariant = ["protocol_id", "primary_outcome", "primary_danger_threshold", "primary_budget",
                 "calibration_rounds", "group_inference", "point_estimator", "scope", "uncertainty_names"]
    for other in summaries[1:]:
        if any(first[key] != other[key] for key in invariant):
            raise ValueError("Stages do not have identical estimands and inference rules")
    strata = sorted([st for summary in summaries for st in summary["strata"]], key=lambda s: (s["n"], s["rel"]))
    keys = [(st["n"], st["rel"]) for st in strata]
    expected = {(n, rel) for n in (50, 200, 1000) for rel in (.005, .01, .03, .05)}
    if len(keys) != 12 or set(keys) != expected:
        raise ValueError("Complete grid requires exactly12 disjoint registered strata")
    return {**{key: first[key] for key in invariant}, "stage": "complete_grid",
            "raw_rows": sum(s["raw_rows"] for s in summaries), "strata": strata,
            "source_stages": {name: {"raw_sha256": value["raw_sha256"], "raw_rows": value["raw_rows"]}
                              for name, value in inputs},
            "combination_rule": "Exact concatenation of12 distinct registered strata; no cross-stratum pooled effect, no independent-replication claim"}


def _pct(value):
    return "undefined" if value is None or not np.isfinite(value) else f"{100*value:.2f}%"


def _mean(pool, metric, policy="F", budget="0.05"):
    return pool["budgets"][budget]["policies"][policy]["metrics"][metric]["point_mean"]


def primary_tables(summary):
    risk, pairs = [], []
    for st in summary["strata"]:
        pool = st["pools"]["confirmation_all"]
        for policy, record in pool["budgets"]["0.05"]["policies"].items():
            metrics = record["metrics"]
            item = {"n": st["n"], "rel": st["rel"], "policy": policy,
                    "feasible_rounds": record["feasible_rounds"]}
            for name in ("dangerous_H1_release_given_screen", "H0_withholding_given_screen", "H0_retention_overall"):
                stat = metrics[name]
                interval = stat.get("crossed")
                item[name] = stat["point_mean"]
                item[f"{name}_ci_low"] = interval["ci95"][0] if interval and interval["ci95"] else None
                item[f"{name}_ci_high"] = interval["ci95"][1] if interval and interval["ci95"] else None
                item[f"{name}_defined_draws"] = interval["defined_draws"] if interval else 0
            item["danger_contributing_groups"] = json.dumps(metrics["dangerous_H1_release_given_screen"]["denominator_contributing_test_groups_by_round"])
            risk.append(item)
        for pair, values in pool["budgets"]["0.05"]["paired"].items():
            stat = values["dangerous_H1_release_given_screen"]
            interval = stat.get("crossed")
            pairs.append({"n": st["n"], "rel": st["rel"], "pair": pair,
                          "difference": stat["point_mean"],
                          "ci_low": interval["ci95"][0] if interval and interval["ci95"] else None,
                          "ci_high": interval["ci95"][1] if interval and interval["ci95"] else None,
                          "defined_draws": interval["defined_draws"] if interval else 0})
    return pd.DataFrame(risk), pd.DataFrame(pairs)


def coverage(test):
    test = test.copy()
    test["screen"] = np.isfinite(test.r2) & test.r2.ge(.99)
    test["danger"] = np.isfinite(test.weighted_mean_relative_error) & test.weighted_mean_relative_error.gt(.1)
    test["p95_danger"] = np.isfinite(test.weighted_p95_relative_error) & test.weighted_p95_relative_error.gt(.1)
    test["danger_screened"] = test.screen & test.danger
    test["p95_danger_screened"] = test.screen & test.p95_danger
    records = []
    for (n, rel), frame in test.groupby(["n", "rel"]):
        h0 = frame[frame.split.eq("confirmation") & frame.klass.eq("H0")]
        h1 = frame[frame.split.eq("confirmation") & frame.klass.eq("H1")]
        physical = h1[h1.source_kind.eq("new_source_positive_restriction")]
        records.append({"n": int(n), "rel": float(rel), "confirmation_h0_total": len(h0),
                        "confirmation_h0_screened": int(h0.screen.sum()),
                        "confirmation_h0_screen_rate": float(h0.screen.mean()),
                        "confirmation_dangerous_screened": int(h1.danger_screened.sum()),
                        "confirmation_danger_structure_groups": h1.loc[h1.danger_screened, "group_unit"].nunique(),
                        "physical_total": len(physical), "physical_screened": int(physical.screen.sum()),
                        "physical_dangerous_screened": int(physical.danger_screened.sum()),
                        "physical_danger_structure_groups": physical.loc[physical.danger_screened, "group_unit"].nunique(),
                        "physical_p95_dangerous_screened": int(physical.p95_danger_screened.sum()),
                        "physical_p95_danger_structure_groups": physical.loc[physical.p95_danger_screened, "group_unit"].nunique()})
    return pd.DataFrame(records)


def figures(out, risk, covered):
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.dpi": 180})
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 7.2), sharex=True, sharey=True)
    for ax, rel in zip(axes.flat, (.005, .01, .03, .05), strict=True):
        subset = risk[risk.rel.eq(rel)]
        for policy, label in POLICIES.items():
            points = subset[subset.policy.eq(policy)].sort_values("n")
            y = points.dangerous_H1_release_given_screen.to_numpy(dtype=float)
            ax.plot(points.n, 100*y, marker="o", color=COLORS[policy], label=label,
                    linewidth=2.2 if policy in {"F", "AEF_panel"} else 1.2,
                    linestyle="--" if policy == "AEF_panel" else "-")
            if policy in {"F", "AEF_panel"}:
                low = points.dangerous_H1_release_given_screen_ci_low.to_numpy(dtype=float)
                high = points.dangerous_H1_release_given_screen_ci_high.to_numpy(dtype=float)
                # The point need not fall inside a percentile interval; draw
                # endpoint segments directly rather than imposing that premise.
                for n, lo, hi in zip(points.n, low, high, strict=True):
                    if np.isfinite(lo) and np.isfinite(hi):
                        ax.vlines(n, lo*100, hi*100, color=COLORS[policy], alpha=.6, linewidth=2)
        ax.set_title(f"Log-noise SD = {rel:g}")
        ax.set_xscale("log")
        ax.set_xticks([50, 200, 1000], ["50", "200", "1000"])
        ax.set_ylim(-2, 102)
        ax.grid(alpha=.15)
    fig.supxlabel("Training sample size", y=.055)
    fig.supylabel("Screened-dangerous H1 released (%)")
    fig.suptitle("Fixed external regions: confirmation groups, 5% calibration budget", y=.98, fontsize=13)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5, bbox_to_anchor=(.5, -.005), frameon=False)
    fig.tight_layout(rect=(.025, .10, 1, .95))
    fig.savefig(out / "p02_full_grid.png", bbox_inches="tight")
    fig.savefig(out / "p02_full_grid.pdf", bbox_inches="tight")
    plt.close(fig)
    f = risk[risk.policy.eq("F")][["n", "rel", "H0_withholding_given_screen"]]
    covered = covered.merge(f, on=["n", "rel"], validate="one_to_one")
    specs = [("confirmation_h0_screen_rate", "H0 passed fit screen (%)", True),
             ("H0_withholding_given_screen", "F: H0 withheld after screen (%)", True),
             ("physical_dangerous_screened", "Physical: primary dangerous + screened (count)", False),
             ("physical_p95_dangerous_screened", "Physical: p95-dangerous + screened (count)", False)]
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.8))
    for ax, (column, title, percent) in zip(axes.flat, specs, strict=True):
        matrix = covered.pivot(index="n", columns="rel", values=column).reindex(index=[50, 200, 1000], columns=[.005, .01, .03, .05]).to_numpy(dtype=float)
        if percent:
            matrix = matrix*100
        maximum = 100 if percent else max(float(np.nanmax(matrix)), 1)
        ax.imshow(matrix, cmap="YlGnBu", vmin=0, vmax=maximum, aspect="auto")
        for (row, col), value in np.ndenumerate(matrix):
            label = "NA" if not np.isfinite(value) else f"{value:.1f}" if percent else f"{value:.0f}"
            ax.text(col, row, label, ha="center", va="center", color="white" if value > maximum*.55 else "#1b2838", fontsize=11)
        ax.set_title(title, fontsize=10.5)
        ax.set_xticks(range(4), [".005", ".01", ".03", ".05"])
        ax.set_yticks(range(3), ["50", "200", "1000"])
        ax.set_xlabel("Log-noise SD")
        ax.set_ylabel("n")
    fig.suptitle("Fit-screen transfer and physical endpoint coverage", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, .95))
    fig.savefig(out / "p02_coverage.png", bbox_inches="tight")
    fig.savefig(out / "p02_coverage.pdf", bbox_inches="tight")
    plt.close(fig)


def report(out, summary, risk, covered):
    lines = ["# P0-2：冻结协议的完整样本量与噪声网格", "",
             "两阶段共246,240份独立数据、108,000条共享测试记录；3种样本量×4档噪声全部完成。每格10轮独立参考／校准，所有政策保留1%／5%／10%预算。本汇总只拼接互不重叠的分层结果，不将不同噪声或阶段当成同条件独立重复。", "",
             "主风险：训练对数空间R²≥0.99且固定外域加权平均相对误差>10%的H1，被完整政策放行的比例。下表预算为相同5%校准预算，实际测试H0成本另行报告。", "",
             "| n | 对数噪声SD | D | F | AD | AEF | G | 主危险分母／结构组 | 物理主危险分母／组 |", "|---|---|---|---|---|---|---|---|---|"]
    for st in summary["strata"]:
        row = covered[covered.n.eq(st["n"]) & covered.rel.eq(st["rel"])].iloc[0]
        points = risk[risk.n.eq(st["n"]) & risk.rel.eq(st["rel"])].set_index("policy")
        values = [_pct(points.loc[p, "dangerous_H1_release_given_screen"]) for p in POLICIES]
        lines.append(f"| {st['n']} | {st['rel']:g} | " + " | ".join(values) + f" | {int(row.confirmation_dangerous_screened)}/{int(row.confirmation_danger_structure_groups)} | {int(row.physical_dangerous_screened)}/{int(row.physical_danger_structure_groups)} |")
    lines += ["", "## 解读与边界", "",
              "- 逐式、窗口及结构组的身份在采样前冻结；这次扩展没有按第一阶段结果改域、挑任务或改变主终点。",
              "- 相同校准预算不等于相同测试误拒绝率；同时看H0筛选通过率、条件拒放率和整体保留率。",
              "- 无危险筛选分母的物理格子，其主风险是undefined。不可报成0%，也不可用次要p95终点替代主终点。",
              "- 物理来源的3600个已观察主危险拟合全部被R²筛选拦住，整个注册工作流的次要overall危险放行是0/3600；这不能归功于后续检验，也不是物理总体安全保证。缺的是高拟合且仍危险的物理候选。",
              "- F在11格主危险放行点估计最低，另1格五种政策均为0；点估计排序不是逐格显著优势、相同测试成本的支配或通用最优。最低噪声n=1000的零值及退化区间不表示总体风险为零。",
              "- 区间按H1结构组／H0设计组与校准轮交叉重采样。参数窗口不是独立机制，测试库跨轮共享。",
              "- 区间条件于固定校准设计池和已实现测试记录；不包含新抽记录或新抽校准结构族的全部不确定性。单组或零分母可产生退化或无定义区间。",
              "- 物理库扩展协议是下一项前瞻研究设计；没有把尚未执行的新物理任务写成验证结果。", "",
              "## 文件", "",
              "- summary.json：完整12格，逐格保留原始点估计、区间、有效抽样次数。",
              "- primary_confirmation.csv、paired_confirmation.csv：5%预算主要确认结果及配对差异。",
              "- coverage.csv：H0筛选与物理主／次要危险分母；primary_metrics.csv与sensitivity_metrics.csv保留全部注册预算和终点。",
              "- source_records.json：原始数据、逐组计数、决策的两阶段路径与哈希，避免复制巨量决策而改变原文件。",
              "- independent_audit.json：联合网格独立审计；stage1和stage2各自审计仍保留在原目录。"]
    (out / "RUN_NOTES.md").write_text("\n".join(lines)+"\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", default="runs/p02_fixed_domain")
    args = parser.parse_args()
    run = (ROOT / args.run).resolve()
    frozen = verify_frozen(run)
    preserved = verify_preservation(run)
    out = run / "complete_grid"
    if (out / "manifest.json").exists():
        manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
        if all(digest(out / name) == expected for name, expected in manifest["output_hashes"].items()):
            print("Existing complete-grid reporting artifacts verified; no overwrite.")
            return
        raise ValueError("Previously completed reporting artifacts changed")
    inputs, frames, sources, policies = [], [], {}, []
    cols = ["base_seed", "dataset_sha256", "phase", "split", "klass", "source_kind", "group_unit", "task", "n", "rel", "r2", "weighted_mean_relative_error", "weighted_p95_relative_error"]
    for stage in ("stage1", "stage2"):
        path = run / stage
        execution = json.loads((path / "manifest.json").read_text(encoding="utf-8"))
        audit = json.loads((path / "independent_audit.json").read_text(encoding="utf-8"))
        if not execution["complete"] or not audit["passed"] or digest(path / "raw.parquet") != execution["raw_sha256"]:
            raise ValueError(f"Stage {stage} is incomplete, unaudited or changed")
        inputs.append((stage, json.loads((path / "summary.json").read_text(encoding="utf-8"))))
        frames.append(pd.read_parquet(path / "raw.parquet", columns=cols))
        policies.extend(json.loads((path / "policies.json").read_text(encoding="utf-8")))
        sources[stage] = {name: {"path": (path / name).relative_to(ROOT).as_posix(), "sha256": digest(path / name)}
                          for name in ["manifest.json", "summary.json", "raw.parquet", "independent_audit.json", "decisions.parquet", "group_counts.parquet", "policies.json", *CSVS]}
    raw = pd.concat(frames, ignore_index=True)
    if raw.base_seed.duplicated().any() or raw.dataset_sha256.dropna().duplicated().any():
        raise ValueError("Stage1/stage2 data overlap")
    if len(raw) != 246240 or len(raw[raw.phase.eq("test")]) != 108000:
        raise ValueError("Incorrect full-grid dataset/test cardinality")
    summary = combine_summaries(inputs)
    covered = coverage(raw[raw.phase.eq("test")])
    risk, pairs = primary_tables(summary)
    out.mkdir(exist_ok=True)
    write_json(out / "summary.json", summary)
    write_json(out / "policies.json", policies)
    write_json(out / "source_records.json", sources)
    for name in CSVS:
        pd.concat([pd.read_csv(run / stage / name, float_precision="round_trip").assign(source_stage=stage) for stage in ("stage1", "stage2")], ignore_index=True).to_csv(out / name, index=False)
    risk.to_csv(out / "primary_confirmation.csv", index=False)
    pairs.to_csv(out / "paired_confirmation.csv", index=False)
    covered.to_csv(out / "coverage.csv", index=False)
    figures(out, risk, covered)
    report(out, summary, risk, covered)
    verify_frozen(run)
    verify_preservation(run)
    manifest = {"complete": True, "protocol_fingerprint": frozen["fingerprint"], "stage": "complete_grid",
                "dataset_count": len(raw), "test_count": int(raw.phase.eq("test").sum()), "strata": 12,
                "policy_configurations": len(policies), "source_records": sources, "preservation": preserved,
                "reporting_source_sha256": digest(Path(__file__)),
                "output_hashes": {p.name: digest(p) for p in sorted(out.iterdir()) if p.is_file() and p.name != "manifest.json"}}
    write_json(out / "manifest.json", manifest)
    print(json.dumps({"complete": True, "datasets": len(raw), "test_records": manifest["test_count"], "strata": 12, "out": str(out)}))


if __name__ == "__main__":
    main()
