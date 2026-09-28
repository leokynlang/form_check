"""srsd_analyze.py — 合并 SRSD 分片 → summary_srsd.json + task_registry.csv
+ report_srsd.md + figures。

四重点问题的量化落点:
Q1 单变量二次偏离低自由度功效 → 合成台架 QUAD_LIN(报告引用)+ SRSD 内
   form_check vs RESET 在 1-2 变量任务上的功效差;
Q2 平方平均是否遗漏特征特定曲率/交互 → 多变量 H1 任务逐任务功效差
   (form_check vs exog/fitted),点名交互类任务;
Q3 HC3 是否改善异方差 FPR → H0 池 fitted_f vs fitted_hc3 vs exog_hc3
   误报率按 n 分解(含 fitted Wald 共线脆弱性现象);
Q4 报警预测外推误差 → p 值对"远端外推相对误差>10%"的 rank-AUC。
"""
import argparse
import glob
import json
import os
import sys

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = os.path.dirname(os.path.abspath(__file__))
METHODS = {"fc": "form_check", "fitted_f": "RESET fitted F",
           "fitted_hc3": "RESET fitted HC3", "exog_hc3": "RESET exog HC3"}
NS = [50, 200, 1000]


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    den = 1 + z ** 2 / n
    c = (p + z ** 2 / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / den
    return (max(0.0, c - h), min(1.0, c + h))


def rank_auc(p, y):
    """Mann-Whitney rank AUC;y=1 为正类。"""
    ok = np.isfinite(p)
    p, y = p[ok], y[ok]
    n1, n0 = int(y.sum()), int((1 - y).sum())
    if not n1 or not n0:
        return np.nan
    order = np.argsort(p)
    ranks = np.empty(len(p), dtype=float)
    ranks[order] = np.arange(1, len(p) + 1)
    # 并列处理:平均秩
    sp = pd.Series(p).rank(method="average").values
    return float((sp[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="srsd_pilot",
                    help="合并输出目录名")
    ap.add_argument("--shards", default="srsd_p0,srsd_p1,srsd_p2,srsd_p3,srsd_p4")
    ap.add_argument("--seeds", type=int, default=50)
    args = ap.parse_args()

    frames = []
    for shard in args.shards.split(","):
        for f in glob.glob(os.path.join(HERE, "runs", shard, "raw_srsd.parquet")):
            frames.append(pd.read_parquet(f))
    if not frames:
        print("no shard data yet"); return
    raw = pd.concat(frames, ignore_index=True)
    out_dir = os.path.join(HERE, "runs", args.tag)
    os.makedirs(out_dir, exist_ok=True)
    raw.to_parquet(os.path.join(out_dir, "raw_srsd.parquet"), index=False)

    # 数据完整性:每格种子数
    cnt = raw.groupby(["task", "n", "rel"]).size()
    summary = {"n_rows": len(raw),
               "cells": int(cnt.groupby(level=[0, 1, 2]).ngroups),
               "min_seeds_per_cell": int(cnt.min()) if len(cnt) else 0}

    # ---------- H0 FPR(仅 rel>=0.005;按 n) ----------
    h0 = raw[(raw.klass == "H0") & (raw.rel >= 0.005)]
    h0_fpr, h0_fpr_ci = {}, {}
    for n in NS:
        h0_fpr[str(n)], h0_fpr_ci[str(n)] = {}, {}
        for m in METHODS:
            v = h0[h0.n == n][f"rej_{m}"].dropna().values.astype(bool)
            lo, hi = wilson(int(v.sum()), len(v))
            h0_fpr[str(n)][m] = float(v.mean()) if len(v) else np.nan
            h0_fpr_ci[str(n)][m] = [float(lo), float(hi)]
    summary["h0_fpr"] = h0_fpr
    summary["h0_fpr_ci"] = h0_fpr_ci

    # ---------- 经验 FPR=5% 校准阈值(每 n,H0 池) ----------
    tau = {}
    for n in NS:
        sub = h0[h0.n == n]
        tau[n] = {m: float(np.nanquantile(sub[f"p_{m}"].values, 0.05))
                  for m in METHODS}
    summary["tau_at_fpr5"] = {str(k): v for k, v in tau.items()}

    # ---------- 校准功效(每 n,全 H1 池;及逐任务) ----------
    h1 = raw[(raw.klass == "H1") & (raw.rel >= 0.005)]
    power_cal = {}
    for n in NS:
        power_cal[str(n)] = {}
        for m in METHODS:
            sub = h1[h1.n == n]
            p = sub[f"p_{m}"].values
            ok = np.isfinite(p)
            rej = (p[ok] <= tau[n][m])
            power_cal[str(n)][m] = float(rej.mean()) if ok.any() else np.nan
    summary["power_at_empirical_fpr5"] = power_cal

    # ---------- 逐任务功效 @α=0.05(n=200,rel=1%) ----------
    per_task = []
    for t in sorted(h1.task.unique()):
        row = {"task": t}
        for m in METHODS:
            sel = h1[(h1.task == t) & (h1.n == 200) & (h1.rel == 0.01)]
            v = sel[f"rej_{m}"].dropna().values.astype(bool)
            row[f"power_{m}"] = float(v.mean()) if len(v) else np.nan
        per_task.append(row)
    per_task = pd.DataFrame(per_task)
    per_task.to_csv(os.path.join(out_dir, "h1_power_per_task.csv"), index=False)

    # ---------- 外推 AUC(H1+H0 全部,rel>=0.005) ----------
    ex = raw[(raw.rel >= 0.005) & raw.extrap_far_relerr.notna()].copy()
    y = (ex.extrap_far_relerr > 0.10).values.astype(int)
    extrap_auc = {m: rank_auc(ex[f"p_{m}"].values, y) for m in METHODS}
    # 按 rel 分层(池化会混入"高噪声→高误差但低功效"的 Simpson 伪影)
    extrap_auc_by_rel = {}
    for rel in sorted(ex.rel.unique()):
        sub = ex[ex.rel == rel]
        yy = (sub.extrap_far_relerr > 0.10).values.astype(int)
        extrap_auc_by_rel[float(rel)] = {m: rank_auc(sub[f"p_{m}"].values, yy)
                                         for m in METHODS}
    extrap_auc_alarm = {}
    for m, m_name in METHODS.items():
        v = ex[f"rej_{m}"].values.astype(bool)
        tp = int((v & (y == 1)).sum()); fp = int((v & (y == 0)).sum())
        fn = int((~v & (y == 1)).sum()); _tn = int((~v & (y == 0)).sum())
        extrap_auc_alarm[m_name] = {
            "precision": tp / (tp + fp) if tp + fp else np.nan,
            "recall": tp / (tp + fn) if tp + fn else np.nan}
    summary["extrap_auc"] = extrap_auc
    summary["extrap_auc_by_rel"] = extrap_auc_by_rel
    summary["extrap_alarm"] = extrap_auc_alarm
    summary["extrap_pos_rate"] = float(y.mean())

    # ---------- Q3:HC3 的 FPR 改善 fitted_f vs fitted_hc3(按 n) ----------
    summary["hc3_vs_classical_fpr"] = {
        str(n): {"classical": h0_fpr[str(n)]["fitted_f"],
                 "hc3": h0_fpr[str(n)]["fitted_hc3"],
                 "exog_hc3": h0_fpr[str(n)]["exog_hc3"]}
        for n in NS}

    # ---------- 胜负矩阵(逐格,2pp) ----------
    wmat = pd.DataFrame(0, index=[METHODS[m] for m in METHODS],
                        columns=[METHODS[m] for m in METHODS])
    cells = raw[raw.rel >= 0.005].groupby(["task", "klass", "n", "rel"])
    for _, g in cells:
        rates = {m: g[f"rej_{m}"].dropna().mean() for m in METHODS}
        for a, a_name in METHODS.items():
            for b, b_name in METHODS.items():
                if a != b and rates[a] > rates[b] + 0.02:
                    wmat.loc[a_name, b_name] += 1
    wmat.to_csv(os.path.join(out_dir, "win_matrix_srsd.csv"))

    # ---------- 关键问题(数据驱动,报告里展开) ----------
    q2_gap = []
    multi = per_task[per_task.task.str.contains(r"i\.26|i\.27|i\.30|i\.37|i\.39|i\.41|i\.43|i\.48|ii\.|iii\.", regex=True)]
    for _, r in multi.iterrows():
        if np.isfinite(r["power_fc"]) and np.isfinite(r["power_exog_hc3"]):
            q2_gap.append((r["task"], r["power_fc"], r["power_exog_hc3"]))
    q2_gap.sort(key=lambda x: x[2] - x[1], reverse=True)
    summary["key_questions"] = {
        "Q1_univariate_low_df": (
            "SRSD 池无单变量非单项式任务(1 变量物理式几乎全是单项式);"
            "该问题由合成台架 QUAD_LIN 回答(见合成报告 §2)。"),
        "Q2_averaging_blindness": (
            "多变量 H1 任务中 exog-HC3 优于 form_check 最多的前 3:"
            + "; ".join(f"{t}(fc={a:.2f} vs exog={b:.2f})"
                        for t, a, b in q2_gap[:3])),
        "Q3_hc3_hetero_fpr": {
            str(n): summary["hc3_vs_classical_fpr"][str(n)] for n in NS},
        "Q4_alarm_extrapolation": {
            "auc": extrap_auc, "alarm_pr": extrap_auc_alarm},
    }

    with open(os.path.join(out_dir, "summary_srsd.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=1)

    # ---------- 任务注册表(入选/排除清单及理由,交付物) ----------
    from srsd_data import registry_with_meta
    reg_rows = []
    for r in registry_with_meta():
        reg_rows.append({"eq_name": r["eq_name"], "klass": r["klass"],
                         "reason": r["reason"], "n_vars": r["n_vars"],
                         "sampling": str(r["sampling_ranges"]),
                         "sympy": r["sympy_str"]})
    reg_df = pd.DataFrame(reg_rows)
    reg_df.to_csv(os.path.join(out_dir, "task_registry.csv"), index=False)
    summary["n_tasks"] = int((reg_df.klass != "EXCLUDE").sum())
    summary["n_h0"] = int((reg_df.klass == "H0").sum())
    summary["n_h1"] = int((reg_df.klass == "H1").sum())
    summary["n_excluded"] = int((reg_df.klass == "EXCLUDE").sum())
    summary["exclude_reasons"] = (reg_df[reg_df.klass == "EXCLUDE"]
                                  .reason.value_counts().to_dict())
    with open(os.path.join(out_dir, "summary_srsd.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=1)

    # ---------- 图 ----------
    fig_dir = os.path.join(HERE, "figures")
    # 校准曲线(按 n 分面)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    for ax, n in zip(axes, NS, strict=True):
        for m, m_name in METHODS.items():
            p = np.sort(h0[h0.n == n][f"p_{m}"].dropna().values)
            if len(p):
                ax.plot(p, np.arange(1, len(p) + 1) / len(p), label=m_name)
        ax.plot([0, 1], [0, 1], "k--", lw=0.8)
        ax.set_title(f"H0 (SRSD true monomials), n={n}")
        ax.set_xlabel("p"); ax.set_xlim(0, 1)
    axes[0].set_ylabel("empirical CDF"); axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, "srsd_calibration.png"), dpi=130)
    plt.close(fig)

    # 逐任务功效(n=200, rel=1%)横向条形
    if len(per_task):
        fig, ax = plt.subplots(figsize=(9, 0.32 * len(per_task) + 1.6))
        ypos = np.arange(len(per_task))
        w = 0.2
        for j, m in enumerate(METHODS):
            ax.barh(ypos + (j - 1.5) * w, per_task[f"power_{m}"].values,
                    height=w, label=METHODS[m])
        ax.set_yticks(ypos, per_task.task, fontsize=7)
        ax.set_xlabel("power @α=0.05 (n=200, rel=1%)")
        ax.legend(fontsize=7, loc="lower right")
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, "srsd_per_task_power.png"), dpi=130)
        plt.close(fig)

    print(json.dumps({k: summary[k] for k in
                      ["n_rows", "cells", "min_seeds_per_cell",
                       "h0_fpr", "power_at_empirical_fpr5", "extrap_auc"]},
                     ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
