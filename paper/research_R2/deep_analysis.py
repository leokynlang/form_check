"""deep_analysis.py — 两组补充分析，为桌面报告提供数据。

A. 论文 C.1 的 far 阈值敏感性
   「481 个危险且沉默」在不同 far 阈值下的规模与任务构成。
   检验 bonus.15 贡献是否为阈值切分效应。

B. PySR 实验的门卫贡献分解
   每个门卫单独拒绝多少、联合抓住多少、哪个门卫是唯一防线。

输出：<wide_svd>/deep_analysis.json + 控制台表格
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = Path(__file__).resolve().parent
assert OUT_DIR == ROOT / "paper" / "research_R2", "Unexpected script directory"
WIDE = ROOT / "runs/pysr_gate/wide_svd"
OUT = OUT_DIR / "evidence" / "deep_analysis.json"

# ---------------------------------------------------------------- A. 论文侧
def paper_side():
    raw = pd.read_parquet(ROOT / "runs/srsd_full/raw_srsd.parquet")
    q = pd.read_parquet(ROOT / "runs/srsd_full/quality_recompute.parquet")
    cd = pd.read_csv(ROOT / "runs/cd_svd/s2_cd_svd.csv", low_memory=False)
    df = raw.merge(q[["task", "n", "rel", "seed", "r2"]],
                   on=["task", "n", "rel", "seed"], how="left")
    df = df.merge(cd[["task", "n", "rel", "seed", "rej_C_svd", "rej_D_svd",
                      "fail_C_svd", "fail_D_svd"]],
                  on=["task", "n", "rel", "seed"], how="left")
    gates = ["rej_fc", "rej_fitted_f", "rej_C_svd", "rej_D_svd"]
    fail = df[["fail_fc", "fail_fitted_f", "fail_C_svd", "fail_D_svd"]].notna()
    df["_sil"] = (~df[gates].fillna(False).any(axis=1)) & (~fail.any(axis=1))

    base = ((df.klass == "H1") & (df.rel > 0) & df.r2.notna()
            & df.extrap_far_relerr.notna())
    s99 = df[base & (df.r2 >= 0.99)]

    result = {"screened": int(len(s99)), "thresholds": []}
    for th in (0.05, 0.075, 0.10, 0.125, 0.15, 0.20):
        d = s99[s99.extrap_far_relerr > th]
        js = d[d._sil]
        by_task = js.groupby("task").size().sort_values(ascending=False)
        entry = {
            "threshold": th,
            "dangerous": int(len(d)),
            "dangerous_frac_of_screened": round(len(d) / len(s99), 4),
            "silent": int(len(js)),
            "release_rate_of_dangerous": round(len(js) / len(d), 4) if len(d) else 0.0,
            "by_task": {k: int(v) for k, v in by_task.items()},
        }
        # bonus.15 占比
        b15 = int(by_task.get("feynman-bonus.15", 0))
        entry["bonus15_share"] = round(b15 / len(js), 4) if len(js) else 0.0
        result["thresholds"].append(entry)

    # 全部危险 H1（不限 R2）在基准阈值下的分解
    dall = df[base & (df.extrap_far_relerr > 0.10)]
    by_task_all = dall.groupby("task").agg(
        dangerous=("task", "size"), silent=("_sil", "sum"))
    by_task_all["release_rate"] = (by_task_all.silent / by_task_all.dangerous).round(4)
    result["dangerous_all_h1_at_10pct"] = by_task_all.reset_index().to_dict("records")
    return result


# ---------------------------------------------------------------- B. PySR 侧
GATES = {"A": "p_fc", "B": "p_fitted_f",
         "C": "p_fitted_hc3_svd", "D": "p_exog_hc3_svd"}


def pysr_side():
    h1 = pd.DataFrame([json.loads(line) for line in
                       open(WIDE / "h1_runs_wide_svd_far.jsonl",
                            encoding="utf-8")])
    h0 = pd.DataFrame([json.loads(line) for line in
                       open(WIDE / "h0_pool_wide_svd_far.jsonl",
                            encoding="utf-8")])
    for df in (h1, h0):
        for c in list(GATES.values()) + ["r2_log", "far_paper"]:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    h1 = h1[np.isclose(h1.rel, 0.01)].drop_duplicates(["task", "seed"])
    h0 = h0[np.isclose(h0.rel, 0.01)].drop_duplicates(["task", "seed"])

    tau = {g: float(np.quantile(h0[c].dropna(), 0.05)) for g, c in GATES.items()}
    rej = {g: (h1[c] < tau[g]).fillna(True) for g, c in GATES.items()}
    silent = ~np.logical_or.reduce([r.values for r in rej.values()])

    n = len(h1)
    contrib = {}
    for g in GATES:
        contrib[g] = {
            "reject_n": int(rej[g].sum()),
            "reject_rate": round(float(rej[g].mean()), 4),
            "sole_defense_n": int(((rej[g].values) & ~np.logical_or.reduce(
                [rej[k].values for k in GATES if k != g])).sum()),
        }
    contrib["joint_silent_n"] = int(silent.sum())
    contrib["joint_silent_rate"] = round(float(silent.mean()), 4)
    contrib["tau"] = {g: round(v, 4) for g, v in tau.items()}

    # 沉默样本的分层构成
    good = (h1["r2_log"] >= 0.99).values
    danger = (h1["far_paper"] > 0.10).values
    layers = {}
    for name, mask in (("全部", np.ones(n, bool)),
                       ("R2>=0.99", good),
                       ("far>10%", danger),
                       ("R2>=0.99且far>10%", good & danger)):
        if mask.sum() == 0:
            layers[name] = {"n": 0}
            continue
        layers[name] = {
            "n": int(mask.sum()),
            "silent_n": int(silent[mask].sum()),
            "silent_rate": round(float(silent[mask].mean()), 4),
        }
    # 危险候选明细
    dm = good & danger
    detail = h1[dm].assign(
        silent=silent[dm],
        reject_A=(rej["A"].values[dm]), reject_B=(rej["B"].values[dm]),
        reject_C=(rej["C"].values[dm]), reject_D=(rej["D"].values[dm]),
    )[["task", "seed", "r2_log", "far_paper", "p_fc", "p_fitted_f",
       "p_fitted_hc3_svd", "p_exog_hc3_svd", "silent",
       "reject_A", "reject_B", "reject_C", "reject_D", "pysr_expr"]]

    # far>10% 的所有样本里有多少被各门卫单独抓住
    danger_all = h1[h1.far_paper > 0.10]
    danger_contrib = {}
    for g in GATES:
        v = (danger_all[GATES[g]] < tau[g]).fillna(True)
        danger_contrib[g] = int(v.sum())
    danger_contrib["n"] = int(len(danger_all))

    return {"gate_contribution": contrib, "layers": layers,
            "danger_candidates": detail.to_dict("records"),
            "danger_layer_gate_catch": danger_contrib}


if __name__ == "__main__":
    res = {"paper_side_threshold_sensitivity": paper_side(),
           "pysr_side": pysr_side()}
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2, default=str)

    print("=" * 70)
    print("A. 论文 481 的 far 阈值敏感性（R²≥0.99 层）")
    print(f'  通过筛选 {res["paper_side_threshold_sensitivity"]["screened"]}')
    print(f'  {"阈值":>6} {"危险":>8} {"危险占比":>8} {"沉默放行":>8} '
          f'{"放行率":>8} {"bonus.15占比":>12}')
    for e in res["paper_side_threshold_sensitivity"]["thresholds"]:
        print(f'  {e["threshold"]:>6.3f} {e["dangerous"]:>8} '
              f'{e["dangerous_frac_of_screened"]:>8.3f} {e["silent"]:>8} '
              f'{e["release_rate_of_dangerous"]:>8.4f} '
              f'{e["bonus15_share"]:>12.3f}')
    print()
    print("  阈值 0.10 时各任务的沉默放行数：")
    for e in res["paper_side_threshold_sensitivity"]["thresholds"]:
        if e["threshold"] == 0.10:
            for k, v in e["by_task"].items():
                print(f"    {k:22s} {v}")

    print()
    print("=" * 70)
    print("B. PySR 实验的门卫贡献（rel=0.01）")
    c = res["pysr_side"]["gate_contribution"]
    for g in GATES:
        print(f'  {g}: 单独拒绝 {c[g]["reject_n"]:>3} '
              f'({c[g]["reject_rate"]:.1%})  唯一防线 {c[g]["sole_defense_n"]:>3}')
    print(f'  联合沉默 {c["joint_silent_n"]}/{c["joint_silent_n"] + int(sum(c[g]["reject_n"] for g in GATES) > 0 and 0 or 0)}')
    print(f'  τ = {c["tau"]}')
    print()
    print("  分层：")
    for k, v in res["pysr_side"]["layers"].items():
        if v.get("n", 0):
            print(f'    {k:22s} n={v["n"]:>3}  沉默 {v["silent_n"]:>3} '
                  f'({v["silent_rate"]:.1%})')
    print()
    dl = res["pysr_side"]["danger_layer_gate_catch"]
    print(f'  far>10% 层（全部 {dl["n"]} 个）被各门卫抓住数：')
    for g in GATES:
        print(f'    {g}: {dl[g]}')
    print(f'\n已写入 {OUT}')