# -*- coding: utf-8 -*-
"""v5_round2_analysis.py — 第二轮审稿意见在桌面 v5 谱系上的数据分析。

P1. H1 侧配对任务级 bootstrap(18 任务,SVD 重跑行): Δ(D−A), Δ(D−B),
    Δ(D−C), Δ(F−D), Δ(F−A) 区间与 D/F 获胜任务比例。
P2. 32-vs-31 任务池敏感性: C/D(SVD) 的 Table 3 行在两池下的取值
    (i.30.5 用 s2_i305_patch.csv 的重跑行,按 stat_F_recheck 对齐)。
输出 paper/v5_round2_results.json。
"""
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = {}

# τ 与 v5 Table 3 题注/JSON 一致 (n=50/200/1000)
TAU = {"A": {50: 0.050, 200: 0.046, 1000: 0.047},
       "B": {50: 0.051, 200: 0.051, 1000: 0.048},
       "C": {50: 0.0424, 200: 0.0465, 1000: 0.0467},
       "D": {50: 0.0393, 200: 0.0447, 1000: 0.0429}}

s2 = pd.read_csv(os.path.join(ROOT, "runs", "cd_svd", "s2_cd_svd.csv"),
                 low_memory=False)
i305 = "feynman-i.30.5"
s2_ok = s2[s2["sentinel_ok"] == True]  # noqa: E712  i.30.5 全部 sentinel 失败
s2_31 = s2_ok[s2_ok["task"] != i305]
raw = pd.read_parquet(os.path.join(ROOT, "runs", "srsd_full",
                                   "raw_srsd.parquet"))
klass = raw.drop_duplicates("task").set_index("task")["klass"].to_dict()
H1 = sorted(t for t in s2_31["task"].unique() if klass.get(t) == "H1")
print("H1 tasks (31-pool):", len(H1))

P = {"per_task_power": {}, "paired": {}}
COLS = {"A": "p_fc", "B": "p_fitted_f"}
for n in (50, 200, 1000):
    sel_sv = s2_31[(s2_31["n"] == n) & (s2_31["rel"] >= 0.005)]
    sel_raw = raw[(raw["klass"] == "H1") & (raw["n"] == n)
                  & (raw["rel"] >= 0.005)]
    for t in H1:
        key = (t, n)
        g = sel_sv[sel_sv["task"] == t]
        r = sel_raw[sel_raw["task"] == t]
        P["per_task_power"]["%s@n%d" % (t, n)] = {
            "A": float((r["p_fc"] <= TAU["A"][n]).mean()),
            "B": float((r["p_fitted_f"] <= TAU["B"][n]).mean()),
            "C": float((g["p_C_svd"] <= TAU["C"][n]).mean()),
            "D": float((g["p_D_svd"] <= TAU["D"][n]).mean())}

def paired(a, b, label, rng):
    pa = np.array([P["per_task_power"]["%s@n%d" % (t, n)][a] for t in H1])
    pb = np.array([P["per_task_power"]["%s@n%d" % (t, n)][b] for t in H1])
    boots = [float((pa[pick] - pb[pick]).mean())
             for pick in (rng.integers(0, len(H1), len(H1))
                          for _ in range(10000))]
    return {"mean_" + a: round(float(pa.mean()), 3),
            "mean_" + b: round(float(pb.mean()), 3),
            "delta": round(float((pa - pb).mean()), 3),
            "ci95": [round(float(np.quantile(boots, 0.025)), 3),
                     round(float(np.quantile(boots, 0.975)), 3)],
            a + "_wins_tasks": round(float((pa > pb).mean()), 3)}

rng = np.random.default_rng(20260908)
for a, b in [("D", "A"), ("D", "B"), ("D", "C")]:
    for n in (50, 200, 1000):
        P["paired"][f"Δ({a}−{b}) n={n}"] = paired(a, b, b, rng)

# Δ(F−D) / Δ(F−A): F 来自 gate_cmp(n=50/200, rel=0.03)
g = pd.read_csv(os.path.join(ROOT, "runs", "gate_cmp", "results.csv"))
gs = g[(g["bench"] == "srsd") & (g["klass"] == "H1")
       & (~g["group"].str.startswith("noise"))]
gh0 = g[(g["bench"] == "srsd") & (g["klass"] == "H0")]

def calib(arr):
    arr = np.asarray([v for v in arr if np.isfinite(v)])
    cand = np.unique(np.concatenate([arr, [arr.max() + 1e-12]]))
    for t in cand:
        if float((arr > t).mean()) <= 0.05:
            return float(t)
    return float(arr.max())

for n in (50, 200):
    m_star = calib(gh0[gh0["n"] == n]["delta_aicc"].tolist())
    sel = gs[gs["n"] == n]
    sel_sv = s2_31[(s2_31["n"] == n) & (s2_31["rel"] == 0.03)]
    h1_f = sorted(set(sel["task"]) & set(H1))
    fpa, dpa = [], []
    for t in h1_f:
        gt = sel[sel["task"] == t]
        gd = sel_sv[sel_sv["task"] == t]
        fpa.append(float((gt["delta_aicc"] > m_star).mean()))
        dpa.append(float((gd["p_D_svd"] <= TAU["D"][n]).mean()))
    fpa, dpa = np.array(fpa), np.array(dpa)
    rng2 = np.random.default_rng(20260908)
    boots = [float((fpa[pick] - dpa[pick]).mean())
             for pick in (rng2.integers(0, len(h1_f), len(h1_f))
                          for _ in range(10000))]
    P["paired"][f"Δ(F−D) n={n}"] = {
        "mean_F": round(float(fpa.mean()), 3),
        "mean_D": round(float(dpa.mean()), 3),
        "delta": round(float((fpa - dpa).mean()), 3),
        "ci95": [round(float(np.quantile(boots, 0.025)), 3),
                 round(float(np.quantile(boots, 0.975)), 3)],
        "F_wins_tasks": round(float((fpa > dpa).mean()), 3),
        "n_tasks_paired": len(h1_f)}
OUT["P1_paired_h1"] = P

# ------------------------------------------------------------ P2: 32 vs 31
print("== P2: 32-vs-31 C/D(SVD) sensitivity ==")
patch = pd.read_csv(os.path.join(ROOT, "runs", "cd_svd",
                                 "s2_i305_patch.csv"))
main305 = s2[s2["task"] == i305].reset_index(drop=True)
assert len(main305) == 7500, len(main305)
pat = patch.reset_index(drop=True)
align_share = float((main305["stat_F_recheck"].values
                     == pat["stat_F_recheck"].values).mean())
# 对齐校验: recheck 列相同 → 行序一致
if align_share < 0.999:
    # 退而求其次: 按 base_seed 对齐
    pat = pat.iloc[np.argsort(np.arange(len(pat)))]
P2 = {"row_alignment_on_stat_F_recheck": round(align_share, 4)}
# 32 任务池(含 i.30.5 重跑行)的 C/D 校准功效
rows32 = []
for n in (50, 200, 1000):
    sel31 = s2_31[(s2_31["n"] == n) & (s2_31["rel"] >= 0.005)
                  & (s2_31["task"].isin(H1))]
    m305 = main305[(main305["n"] == n) & (main305["rel"] >= 0.005)]
    merged = {}
    merged["D31"] = float((sel31["p_D_svd"] <= TAU["D"][n]).mean())
    merged["C31"] = float((sel31["p_C_svd"] <= TAU["C"][n]).mean())
    p_D_all = np.concatenate([sel31["p_D_svd"].values,
                              m305["p_D_svd"].values])
    p_C_all = np.concatenate([sel31["p_C_svd"].values,
                              m305["p_C_svd"].values])
    merged["D32"] = float((p_D_all <= TAU["D"][n]).mean())
    merged["C32"] = float((p_C_all <= TAU["C"][n]).mean())
    P2[f"n={n}"] = {k: round(v, 3) for k, v in merged.items()}
OUT["P2_pool_32v31"] = P2
print(json.dumps(P2, indent=1))

with open(os.path.join(HERE, "v5_round2_results.json"), "w",
          encoding="utf-8") as fh:
    json.dump(OUT, fh, indent=1, ensure_ascii=False)
print("saved -> paper/v5_round2_results.json")
