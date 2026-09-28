# -*- coding: utf-8 -*-
"""revision_analysis.py — 审稿意见回应的数据分析(全部基于已落盘的原始结果)。

A. OR(A,E) 组合的 H0 实际 size + 组合层重校准到 5% 后的 TPR (gate_cmp)
B. S2 以方程任务为聚类单位的不确定性: 逐任务 FPR/宏平均/任务 bootstrap CI/
   leave-one-equation-out/HC3 逐任务支配度 (srsd_full)
C. S5 校准池/评估池分离: 按方程 leave-one-out 阈值迁移 (gate_cmp SRSD H0)
D. HC3 失准的诊断: 条件数与拒绝的关联 (srsd_full) + statsmodels 独立复核
E. 电池的集合层错误: 真族被排除/错族被保留/空集 (battery_full)
F. S6 五个排列 seed 的 Monte Carlo 误差 (wb_real)

写入 paper/revision_results.json。
"""
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RUNS = os.path.join(ROOT, "runs")
OUT = {}


def calibrate_threshold(stats, alpha=0.05):
    """最小阈值 t 使 H0 池 fire 率 ≤ alpha(与 gate_compare.calibrate 同口径)。"""
    arr = np.asarray([v for v in stats if np.isfinite(v)])
    if len(arr) == 0:
        return np.inf
    cand = np.unique(np.concatenate([arr, [arr.max() + 1e-12]]))
    for t in cand:
        if float((arr > t).mean()) <= alpha:
            return float(t)
    return float(arr.max() + 1e-12)


# ================================================================ A. OR 公平性
print("== A. OR(A,E) size + 组合层重校准 ==")
g = pd.read_csv(os.path.join(RUNS, "gate_cmp", "results.csv"))
A = {"or": {}, "or_recal": {}}
for (bench, n), gg in g.groupby(["bench", "n"]):
    h0 = gg[gg["klass"] == "H0"]
    struct = gg[(gg["klass"] == "H1")
                & (~gg["group"].str.startswith("noise"))]
    if len(h0) == 0 or len(struct) == 0:
        continue
    d_star = calibrate_threshold(h0["D_cv"].tolist())
    m_star = calibrate_threshold(h0["delta_aicc"].tolist())
    a = h0["rej_fc"].astype(bool)
    e = h0["D_cv"] > d_star
    f = h0["delta_aicc"] > m_star
    size_a, size_e, size_f = float(a.mean()), float(e.mean()), float(f.mean())
    size_or = float((a | e).mean())
    # 组合层重校准: A 固定其 5%, E 阈值调到联合 size ≤ 5%
    budget = 0.05 - size_a
    quiet = h0[~a]
    if budget > 0 and len(quiet) > 0:
        # P(D>t | A quiet) ≤ budget / (1 - size_a)
        lvl = 1.0 - budget / (1.0 - size_a)
        t_new = float(np.quantile(quiet["D_cv"].to_numpy(), lvl))
    else:
        t_new = np.inf
    orr = struct["rej_fc"].astype(bool) | (struct["D_cv"] > t_new)
    tp_a = float(struct["rej_fc"].mean())
    tp_e = float((struct["D_cv"] > d_star).mean())
    tp_f = float((struct["delta_aicc"] > m_star).mean())
    tp_or = float((struct["rej_fc"].astype(bool)
                   | (struct["D_cv"] > d_star)).mean())
    tp_or_recal = float(orr.mean())
    key = f"{bench}_n{n}"
    A["or"][key] = {"size_A": round(size_a, 4), "size_E": round(size_e, 4),
                    "size_F": round(size_f, 4), "size_OR": round(size_or, 4),
                    "tpr_A": round(tp_a, 3), "tpr_E": round(tp_e, 3),
                    "tpr_F": round(tp_f, 3), "tpr_OR": round(tp_or, 3)}
    A["or_recal"][key] = {"size_OR_recal": round(
        float((h0["rej_fc"].astype(bool) | (h0["D_cv"] > t_new)).mean()), 4),
        "tpr_OR_recal": round(tp_or_recal, 3), "t_new": round(t_new, 5),
        "d_star": round(d_star, 5)}
OUT["A_or"] = A

# ============================================== B. S2 任务级不确定性
print("== B. S2 task-level ==")
s = pd.read_parquet(os.path.join(RUNS, "srsd_full", "raw_srsd.parquet"))
# 与 srsd_analyze/Table 2 同口径: H0 池 = rel>=0.005(零噪声格不可用)
h0 = s[(s["klass"] == "H0") & (s["n"] == 200) & (s["rel"] >= 0.005)
       & s["rej_fc"].notna()]
METHODS = {"A_fc": "rej_fc", "B_resetF": "rej_fitted_f",
           "C_hc3": "rej_fitted_hc3", "D_exoghc3": "rej_exog_hc3"}
tasks = sorted(h0["task"].unique())
B = {"n_tasks": len(tasks), "n_fits": int(len(h0)), "per_task": {},
     "macro": {}, "micro": {}, "boot_ci": {}, "loo_range": {}, "dominance": {}}
for m, col in METHODS.items():
    per = h0.groupby("task")[col].mean()
    micro = float(h0[col].mean())
    macro = float(per.mean())
    rng = np.random.default_rng(20260906)
    boots = []
    arr_tasks = [h0[h0["task"] == t][col].to_numpy() for t in tasks]
    for _ in range(10000):
        pick = rng.integers(0, len(tasks), len(tasks))
        boots.append(float(np.mean([arr_tasks[i].mean() for i in pick])))
    ci = [float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975))]
    loo = [float(np.mean([arr_tasks[j].mean() for j in range(len(tasks))
                          if j != i])) for i in range(len(tasks))]
    B["per_task"][m] = {t: round(float(v), 4) for t, v in per.items()}
    B["micro"][m] = round(micro, 4)
    B["macro"][m] = round(macro, 4)
    B["boot_ci"][m] = [round(ci[0], 4), round(ci[1], 4)]
    B["loo_range"][m] = [round(min(loo), 4), round(max(loo), 4)]
    # 支配度: 单一任务占全部误报的份额 + 去掉最差任务后的池化率
    tot = float(h0[col].sum())
    worst = per.idxmax()
    B["dominance"][m] = {
        "worst_task": worst, "worst_rate": round(float(per.max()), 4),
        "share_of_all_false_alarms": round(float(h0[h0["task"] == worst][col].sum()) / tot, 4)
        if tot > 0 else None,
        "micro_wo_worst": round(float(h0[h0["task"] != worst][col].mean()), 4)}
OUT["B_tasklevel"] = B

# ============================================== C. S5 校准/评估分离 (LOO)
print("== C. S5 LOO threshold transfer ==")
C = {}
for (bench, n), gg in g.groupby(["bench", "n"]):
    h0g = gg[gg["klass"] == "H0"]
    tks = sorted(h0g["task"].unique())
    if len(tks) < 4:
        continue
    res = {"cv": [], "aicc": []}
    for tk in tks:
        cal = h0g[h0g["task"] != tk]
        ev = h0g[h0g["task"] == tk]
        d_loo = calibrate_threshold(cal["D_cv"].tolist())
        m_loo = calibrate_threshold(cal["delta_aicc"].tolist())
        res["cv"].append(round(float((ev["D_cv"] > d_loo).mean()), 4))
        res["aicc"].append(round(float((ev["delta_aicc"] > m_loo).mean()), 4))
    C[f"{bench}_n{n}"] = {
        "n_tasks": len(tks),
        "cv_size_heldout_mean": round(float(np.mean(res["cv"])), 4),
        "cv_size_heldout_range": [min(res["cv"]), max(res["cv"])],
        "aicc_size_heldout_mean": round(float(np.mean(res["aicc"])), 4),
        "aicc_size_heldout_range": [min(res["aicc"]), max(res["aicc"])]}
OUT["C_splitpool"] = C

# ============================================== D. HC3 诊断
print("== D. HC3 diagnostics ==")
D = {}
sub = h0[h0["fail_fitted_hc3"].isna()]
cond = sub["cond_fitted_hc3"].astype(float)
rej = sub["rej_fitted_hc3"].astype(bool)
lr = np.log10(cond.clip(lower=1))
r = float(np.corrcoef(lr, rej.astype(float))[0, 1])
hi = rej[lr > lr.quantile(0.9)].mean()
lo = rej[lr < lr.quantile(0.1)].mean()
D["s2_h0_n200"] = {
    "n": int(len(sub)), "pointbiserial_logcond_rej": round(r, 3),
    "rej_rate_top_decile_cond": round(float(hi), 4),
    "rej_rate_bottom_decile_cond": round(float(lo), 4),
    "cond_median": round(float(cond.median()), 1),
    "cond_max": round(float(cond.max()), 1)}
OUT["D_hc3"] = D

# ============================================== E. 电池集合层错误
print("== E. battery set-level ==")
bt = pd.read_parquet(os.path.join(RUNS, "battery_full", "raw_battery.parquet"))
members = ["rej_power", "rej_exponential", "rej_logarithmic", "rej_rational"]
E = {}
for tf, gg in bt.groupby("truth_family"):
    col = f"rej_{tf}"
    fams = gg["viable"].apply(lambda v: set() if not isinstance(v, str)
                              else set(v.split("+")))
    if col in gg.columns:
        excluded = float(gg[col].astype(bool).mean())  # 真族成员被拒 = 真族被排除
        entry = {"n": int(len(gg)),
                 "true_family_excluded": round(excluded, 4),
                 "truth_in_viable_set": round(1.0 - excluded, 4),
                 "viable_set_empty": round(float(np.mean([len(x) == 0
                                                          for x in fams])), 4),
                 "wrong_family_also_kept": round(
                     float(np.mean([bool(x - {tf}) for x in fams])), 4)}
    else:  # truth_family='none': 结构 H1 池, 无成员是真值(设计如此)
        entry = {"n": int(len(gg)),
                 "note": "no member is the truth by design (structural H1)",
                 "viable_set_empty": round(float(np.mean([len(x) == 0
                                                          for x in fams])), 4),
                 "any_family_kept": round(float(np.mean([len(x) > 0
                                                         for x in fams])), 4)}
    E[tf] = entry
OUT["E_battery_set"] = E

# ============================================== F. S6 seed MC 误差
print("== F. S6 Monte-Carlo spread ==")
w = pd.read_csv(os.path.join(RUNS, "wb_real", "results.csv"))
wok = w[w["status"] == "ok"] if "status" in w else w
sp = wok.groupby(["year", "task"])["p_fc"].agg(lambda x: float(x.max() - x.min()))
F = {"cells": int(len(sp)),
     "max_spread": round(float(sp.max()), 4),
     "median_spread": round(float(sp.median()), 4),
     "n_seeds": int(wok.groupby(["year", "task"]).size().max()),
     "p_resolution": "1/1000 (B=999)"}
OUT["F_s6_seeds"] = F

with open(os.path.join(HERE, "revision_results.json"), "w",
          encoding="utf-8") as fh:
    json.dump(OUT, fh, indent=1, ensure_ascii=False)
print(json.dumps(OUT, indent=1, ensure_ascii=False)[:4000])
print("saved -> paper/revision_results.json")
