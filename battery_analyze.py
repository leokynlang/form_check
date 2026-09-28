# -*- coding: utf-8 -*-
"""battery_analyze.py — 电池正式跑汇总:校准/判别/升级三层,逐场景表。

读 runs/<tag>/raw_battery.parquet,产出:
- summary_battery.json:size(合成 H0+噪声层上 power 族)、判别率
  (各真值族的 containment/假包含/假升级,合成+SRSD 分开)、
  viable 组合分布;
- report_battery.md:逐场景表 + 三层结论 + 与探针/正式跑的一致性。

语义对齐 battery.py docstring:假排除(漏真族/空集)、假包含(错族留
viable)、正确升级(truth=none 时空集)分开计数,绝不合并成一个数字。
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = os.path.dirname(os.path.abspath(__file__))

FAMS = ["power", "exponential", "logarithmic", "rational"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="battery_full")
    args = ap.parse_args()
    run_dir = os.path.join(HERE, "runs", args.tag)
    raw = pd.read_parquet(os.path.join(run_dir, "raw_battery.parquet"))
    ok = raw[raw.data_fail.isna()].copy()
    ok["n_viable"] = ok["n_viable"].astype(int)
    ok["viable_set"] = ok.viable.apply(
        lambda s: set(s.split(",")) if s else set())
    ok["escalation"] = ok.n_viable == 0

    def contained(row):
        t = row["truth_family"]
        return (t in row["viable_set"]) if t != "none" else None

    ok["containment"] = ok.apply(contained, axis=1)

    def false_included(row):
        t = row["truth_family"]
        if t == "none":
            return None
        return sorted(row["viable_set"] - {t})

    ok["false_included"] = ok.apply(false_included, axis=1)

    summary = {"tag": args.tag, "n_rows": int(len(raw)),
               "n_data_fail": int(raw.data_fail.notna().sum())}

    # ---------- 层 1:合成校准(结构 H0 + 噪声层,power 族尺寸) ----------
    h0 = ok[(ok.side == "synth") & (ok.truth_family == "power")]
    size_by_scen = {name: float(g.rej_power.mean())
                    for name, g in h0.groupby("name")}
    summary["synth_power_size"] = size_by_scen
    summary["synth_power_size_pooled"] = float(h0.rej_power.mean())

    # ---------- 层 2:合成判别(EXP/LOG/RAT 真值族) ----------
    synth_h1 = ok[(ok.side == "synth") & (ok.truth_family != "power")]
    disc = {}
    for name, g in synth_h1.groupby("name"):
        truth = g.truth_family.iloc[0]
        fi = {}
        for _, gg in g.groupby("rel"):
            combos = {}
            for lst in gg.false_included.dropna():
                key = "+".join(lst) if lst else "(none)"
                combos[key] = combos.get(key, 0) + 1
            fi[str(gg["rel"].iloc[0])] = {
                "containment": float(gg.containment.dropna().astype(bool).mean()),
                "escalation": float(gg.escalation.mean()),
                "false_included_top": dict(sorted(
                    combos.items(), key=lambda kv: -kv[1])[:3])}
        disc[name] = {"truth": truth, "by_rel": fi}
    summary["synth_discrimination"] = disc

    # truth=none 场景:升级是正确结局 + 假包含分布
    none_synth = ok[(ok.side == "synth") & (ok.truth_family == "none")]
    esc_ok, fi_none = {}, {}
    for name, g in none_synth.groupby("name"):
        esc_ok[name] = float(g.escalation.mean())
        combos = {}
        for s in g.viable_set:
            if s:
                key = "+".join(sorted(s))
                combos[key] = combos.get(key, 0) + 1
        fi_none[name] = {"false_inclusion_rate":
                         float((g.n_viable > 0).mean()),
                         "top_combos": dict(sorted(combos.items(),
                                                   key=lambda kv: -kv[1])[:3])}
    summary["synth_none_escalation"] = esc_ok
    summary["synth_none_false_inclusion"] = fi_none

    # ---------- 层 3:SRSD(每方程 containment 假设=none,只报 viable/升级) ----------
    srsd = ok[ok.side == "srsd"]
    srsd_summary = {}
    if len(srsd):
        for name, g in srsd.groupby("name"):
            srsd_summary[name] = {
                "escalation": float(g.escalation.mean()),
                "viable_top": {k: round(float(v) / len(g), 3)
                               for k, v in g.viable.value_counts()
                               .head(3).items()},
                "median_n_viable": float(g.n_viable.median())}
    summary["srsd"] = srsd_summary

    with open(os.path.join(run_dir, "summary_battery.json"), "w",
              encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=1)

    # ---------- 报告 ----------
    L, A = [], None
    def A_(s=""):
        L.append(s)
    A = A_
    A_(f"# 形式合格性电池 — {args.tag}\n")
    A_(f"行数 {summary['n_rows']},数据层失败 {summary['n_data_fail']}(显式登记)。")
    A_("wild bootstrap 零分布(Rademacher);α=0.05 每族独立;viable=未拒绝族。\n")
    A_("## 层1 合成校准(power 族尺寸,结构 H0+噪声层)\n")
    A_("| 场景 | 尺寸 |")
    A_("|---|---|")
    for name, v in size_by_scen.items():
        A_(f"| {name} | {v:.3f} |")
    A_(f"\n池化:**{summary['synth_power_size_pooled']:.3f}**(目标 0.05)\n")
    A_("## 层2 合成判别(EXP/LOG/RAT)\n")
    for name, d in disc.items():
        A_(f"### {name}(truth={d['truth']})")
        for rel, v in d["by_rel"].items():
            A_(f"- rel={rel}: containment={v['containment']:.3f} "
               f"假升级={v['escalation']:.3f} "
               f"假包含={v['false_included_top']}")
    A_("\n## 层2b truth=none 场景(升级正确结局)\n")
    A_("| 场景 | 正确升级率 | 假包含率 | top 假包含组合 |")
    A_("|---|---|---|---|")
    for name in esc_ok:
        v = fi_none[name]
        top = next(iter(v["top_combos"].items()), ("-", 0))
        A_(f"| {name} | {esc_ok[name]:.3f} | {v['false_inclusion_rate']:.3f} "
           f"| {top[0]}({top[1]}) |")
    if srsd_summary:
        A_("\n## 层3 SRSD(逐方程,H1 语义=必升级)\n")
        A_("| 方程 | 升级率 | 中位|V| | top viable |")
        A_("|---|---|---|---|")
        for name, v in sorted(srsd_summary.items()):
            top = next(iter(v["viable_top"].items()), ("-", 0))
            A_(f"| {name} | {v['escalation']:.3f} | "
               f"{v['median_n_viable']:.1f} | {top[0]}({top[1]}) |")
    out = os.path.join(run_dir, "report_battery.md")
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(L))
    print("written", out)
    print(json.dumps({k: summary[k] for k in
                      ("synth_power_size_pooled",)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
