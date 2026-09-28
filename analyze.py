"""analyze.py — 汇总 raw.parquet → summary.json + figures + report.md。

要点:
- 原生 alpha=0.05 拒绝率 = 各方法在同一显著性水平下的工作点;
- power@empFPR≤5%: 对每个方法在固定 n 下,用全部 H0(单项式)格的 p 值
  经验分位数校准阈值 τ(使 FPR=5%),再在 H1 格上算功效——
  这回答"相同经验 FPR 下谁更有功"(约束 10);
- 噪声层场景(结构=严格幂律,噪声假设违反)单独呈现:它们的拒绝
  对"形式"而言是误报,不能与结构 H1 功效混在一张总平均里(约束:
  不得仅依据总平均宣布胜者);
- Wilson 95% CI;胜-负矩阵按 (scenario,n,rel) 格逐对比较。
"""
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
H0_SET = {"H0_BASELINE", "H0_MULTIVAR", "H0_CORR"}
NOISE_H1 = {"H1_HEAVYTAIL", "H1_VAR_GROWING", "H1_AUTOCORR",
            "H1_OUTLIER", "H1_OUTLIER_GRP"}
STRUCT_H1_ORDER = ["QUAD_LIN", "CUBIC_LIN", "INTERACT", "THRESHOLD",
                   "EXP_MODEL", "LOG_MODEL", "RATIONAL", "OMITTED_M",
                   "DILUTION", "OPPOSING", "COLLINEAR_CURV"]


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    den = 1 + z ** 2 / n
    center = (p + z ** 2 / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / den
    return (max(0.0, center - half), min(1.0, center + half))


def rate_ci(mask):
    n = len(mask)
    k = int(np.sum(mask))
    lo, hi = wilson(k, n)
    return k / n, lo, hi, n


def main(tag, shards=None):
    # shards: 逗号分隔的多个 run tag,各自有 raw.parquet;合并后统一分析,
    # 输出目录用首个 tag(--shards full_r0,full_r1,... --tag full)。
    run_dir = os.path.join(HERE, "runs", tag)
    tags = [tag] if not shards else [t.strip() for t in shards.split(",") if t.strip()]
    raws, extraps = [], []
    for t in tags:
        p = os.path.join(HERE, "runs", t, "raw.parquet")
        if not os.path.exists(p):
            print(f"[warn] 缺少 {p},跳过该分片")
            continue
        raws.append(pd.read_parquet(p))
        ep = os.path.join(HERE, "runs", t, "extrap.parquet")
        if os.path.exists(ep):
            extraps.append(pd.read_parquet(ep))
    raw = pd.concat(raws, ignore_index=True)
    extrap = pd.concat(extraps, ignore_index=True) if extraps else None
    os.makedirs(run_dir, exist_ok=True)   # 合并模式下输出目录(=tag)可能不存在
    fig_dir = os.path.join(HERE, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    rej_cols = {m: f"rej_{m}" for m in METHODS}
    p_cols = {m: f"p_{m}" for m in METHODS}
    time_cols = {m: f"time_{m}" for m in METHODS}
    for c in list(rej_cols.values()) + list(p_cols.values()):
        assert c in raw.columns, f"missing column {c}"

    cells = raw.groupby(["scenario", "klass", "n", "rel"], as_index=False)
    cell_stats = []
    for key, g in cells:
        scen, klass, n, rel = key
        row = {"scenario": scen, "klass": klass, "n": int(n), "rel": float(rel)}
        for m in METHODS:
            r, lo, hi, _cnt = rate_ci(g[rej_cols[m]].values.astype(bool))
            row[f"rate_{m}"], row[f"lo_{m}"], row[f"hi_{m}"] = r, lo, hi
            row[f"meanp_{m}"] = float(np.nanmean(g[p_cols[m]].values))
            row[f"time_{m}"] = float(np.nanmean(g[time_cols[m]].values))
            row[f"failrate_{m}"] = float(np.mean(
                g[f"fail_{m}"].notna().values)) if f"fail_{m}" in g else 0.0
        cell_stats.append(row)
    cs = pd.DataFrame(cell_stats)

    # ---------- 经验 FPR 校准阈值(per n, pooling 全部 H0 格) ----------
    calib = {}
    h0 = raw[raw.scenario.isin(H0_SET)]
    for n in sorted(raw.n.unique()):
        sub = h0[h0.n == n]
        tau = {}
        for m in METHODS:
            p = sub[p_cols[m]].dropna().values
            tau[m] = float(np.quantile(p, 0.05)) if len(p) >= 50 else np.nan
        calib[int(n)] = tau

    def power_at_tau(scen_mask, n, m):
        sub = raw[(raw.n == n) & (raw.scenario.isin(scen_mask))]
        p = sub[p_cols[m]].values
        ok = np.isfinite(p)
        if not ok.any() or not np.isfinite(calib[int(n)][m]):
            return np.nan, 0
        rej = p[ok] <= calib[int(n)][m]
        r, _lo, _hi, cnt = rate_ci(rej)
        return r, cnt

    # ---------- 图 1: H0 p 值校准 (per n) ----------
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    for ax, n in zip(axes, [50, 200, 1000], strict=True):
        for m, m_name in METHODS.items():
            p = np.sort(h0[h0.n == n][p_cols[m]].dropna().values)
            if len(p) == 0:
                continue
            ax.plot(p, np.arange(1, len(p) + 1) / len(p), label=m_name)
        ax.plot([0, 1], [0, 1], "k--", lw=0.8)
        ax.set_title(f"H0 pooled, n={n}")
        ax.set_xlabel("p"); ax.set_xlim(0, 1)
    axes[0].set_ylabel("empirical CDF"); axes[0].legend(fontsize=8)
    fig.suptitle("p-value calibration under strict monomial H0 (uniform = calibrated)")
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, f"{tag}_calibration.png"), dpi=130)
    plt.close(fig)

    # ---------- 图 2: 结构 H1 功效热图 (n=200, α=0.05) ----------
    fams = STRUCT_H1_ORDER
    mat = np.zeros((len(fams), len(METHODS)))
    for i, s in enumerate(fams):
        for j, m in enumerate(METHODS):
            sel = cs[(cs.scenario == s) & (cs.n == 200) & (cs.rel == 0.01)]
            mat[i, j] = float(sel[f"rate_{m}"].iloc[0]) if len(sel) else np.nan
    fig, ax = plt.subplots(figsize=(7.5, 6))
    im = ax.imshow(mat, vmin=0, vmax=1, cmap="RdYlGn")
    ax.set_xticks(range(len(METHODS)), METHODS.values(), rotation=20, ha="right")
    ax.set_yticks(range(len(fams)), fams, fontsize=8)
    for i in range(len(fams)):
        for j in range(len(METHODS)):
            ax.text(j, i, f"{mat[i,j]:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, shrink=0.8)
    ax.set_title("Structural-H1 power at α=0.05 (n=200, rel=1%)")
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, f"{tag}_power_heat.png"), dpi=130)
    plt.close(fig)

    # ---------- 图 3: n 扩缩(关键场景) ----------
    key_scens = ["QUAD_LIN", "CUBIC_LIN", "INTERACT", "DILUTION", "OPPOSING"]
    fig, axes = plt.subplots(1, len(key_scens), figsize=(4 * len(key_scens), 3.6),
                             sharey=True)
    for ax, s in zip(axes, key_scens, strict=True):
        for m, m_name in METHODS.items():
            ys = [cs[(cs.scenario == s) & (cs.n == n) & (cs.rel == 0.01)][f"rate_{m}"]
                  for n in [50, 200, 1000]]
            ax.plot([50, 200, 1000], [float(y.iloc[0]) if len(y) else np.nan
                                      for y in ys], "o-", label=m_name)
        ax.set_xscale("log"); ax.set_title(s, fontsize=9)
        ax.set_xlabel("n"); ax.axhline(0.8, color="gray", ls=":", lw=0.7)
    axes[0].set_ylabel("power @α=0.05"); axes[0].legend(fontsize=7)
    fig.suptitle("Power vs n (rel=1%)")
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, f"{tag}_n_scaling.png"), dpi=130)
    plt.close(fig)

    # ---------- 图 4: 运行时间 ----------
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    means = {METHODS[m]: [] for m in METHODS}
    for m in METHODS:
        for n in [50, 200, 1000]:
            sel = cs[(cs.n == n)]
            means[METHODS[m]].append(float(sel[f"time_{m}"].mean()))
    xs = np.arange(3)
    for i, (name, vals) in enumerate(means.items()):
        ax.bar(xs + (i - 1.5) * 0.2, vals, width=0.2, label=name)
    ax.set_yscale("log"); ax.set_xticks(xs, ["n=50", "n=200", "n=1000"])
    ax.set_ylabel("seconds / dataset"); ax.legend(fontsize=8)
    ax.set_title("Runtime per method (mean over cells)")
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, f"{tag}_timing.png"), dpi=130)
    plt.close(fig)

    # ---------- 图 5: 报警 vs 外推误差 ----------
    auc_extrap = {}
    if extrap is not None and len(extrap):
        e1 = extrap[extrap.extrap_far_relerr.notna()]
        fig, axes = plt.subplots(1, 2, figsize=(11, 4))
        ax = axes[0]
        for m, col in [("fc", "rej_fc"), ("fitted_f", "rej_fitted_f")]:
            sel = e1[e1.rel >= 0.005]
            alarm = sel[col].values.astype(bool)
            err = sel.extrap_far_relerr.values
            ax.scatter(err[~alarm] + np.random.uniform(-.004, .004, (~alarm).sum()),
                       np.zeros((~alarm).sum()) + (0 if m == "fc" else 0.15),
                       s=6, alpha=0.3, label=f"{METHODS[m]} no-alarm")
            ax.scatter(err[alarm] + np.random.uniform(-.004, .004, alarm.sum()),
                       np.ones(alarm.sum()) + (0 if m == "fc" else 0.15),
                       s=6, alpha=0.3, label=f"{METHODS[m]} alarm")
        ax.set_xlabel("far extrapolation rel. error"); ax.set_yticks([0, 1])
        ax.set_title("alarm vs extrapolation error"); ax.legend(fontsize=7)
        ax = axes[1]
        for m, col in [("fc", "p_fc"), ("fitted_hc3", "p_fitted_hc3")]:
            if col not in e1.columns:
                continue
            sel = e1[(e1.rel >= 0.005) & e1[col].notna()]
            y_true = (sel.extrap_far_relerr > 0.10).values.astype(int)
            p = sel[col].values
            order = np.argsort(p)
            ranks = np.empty_like(order, dtype=float); ranks[order] = np.arange(1, len(p) + 1)
            n1, n0 = y_true.sum(), (1 - y_true).sum()
            auc = (ranks[y_true == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0) \
                if n1 and n0 else np.nan
            auc_extrap[METHODS[m]] = float(auc)
            ax.plot(p, np.cumsum(y_true[np.argsort(p)]) / max(n1, 1),
                    label=f"{METHODS[m]} AUC={auc:.2f}")
        ax.set_xscale("log"); ax.set_xlabel("p-value")
        ax.set_title("error>10% detection ROC (rank-based)")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(fig_dir, f"{tag}_extrap.png"), dpi=130)
        plt.close(fig)

    # ---------- 胜负矩阵(结构 H1 格逐对,按 n=200) ----------
    win = {m: dict.fromkeys(METHODS, 0) for m in METHODS}
    for s in STRUCT_H1_ORDER:
        for n in [50, 200, 1000]:
            for rel in [0.005, 0.01, 0.03, 0.05]:
                sel = cs[(cs.scenario == s) & (cs.n == n) & (cs.rel == rel)]
                if not len(sel):
                    continue
                rates = {m: float(sel[f"rate_{m}"].iloc[0]) for m in METHODS}
                for a in METHODS:
                    for b in METHODS:
                        if a != b:
                            if rates[a] > rates[b] + 0.02:
                                win[a][b] += 1
                            elif rates[a] < rates[b] - 0.02:
                                win[a][b] -= 0  # 负局只在 b vs a 计
    # 简单对称计:胜=rate 高 2pp 以上
    wmat = pd.DataFrame(0, index=[METHODS[m] for m in METHODS],
                        columns=[METHODS[m] for m in METHODS])
    for s in STRUCT_H1_ORDER:
        for n in [50, 200, 1000]:
            for rel in [0.005, 0.01, 0.03, 0.05]:
                sel = cs[(cs.scenario == s) & (cs.n == n) & (cs.rel == rel)]
                if not len(sel):
                    continue
                rates = {m: float(sel[f"rate_{m}"].iloc[0]) for m in METHODS}
                for a, a_name in METHODS.items():
                    for b, b_name in METHODS.items():
                        if a != b and rates[a] > rates[b] + 0.02:
                            wmat.loc[a_name, b_name] += 1

    # ---------- summary.json ----------
    summary = {
        "tag": tag, "n_rows": len(raw),
        "calibrated_tau_at_fpr5": calib,
        "h0_false_alarm_alpha05": {
            int(n): {m: {
                "rate": float(cs[(cs.scenario.isin(H0_SET)) & (cs.n == n)][f"rate_{m}"].mean()),
            } for m in METHODS} for n in sorted(raw.n.unique())},
        "structural_h1_power_alpha05": {
            s: {int(n): {m: float(cs[(cs.scenario == s) & (cs.n == n) & (cs.rel == 0.01)][f"rate_{m}"].iloc[0])
                         for m in METHODS} for n in [50, 200, 1000]}
            for s in STRUCT_H1_ORDER},
        "noise_layer_rejection_alpha05": {
            s: {int(n): {m: float(cs[(cs.scenario == s) & (cs.n == n) & (cs.rel == 0.01)][f"rate_{m}"].iloc[0])
                         for m in METHODS} for n in [50, 200, 1000]}
            for s in sorted(NOISE_H1)},
        "power_at_empirical_fpr5_rel1pct": {
            int(n): {m: power_at_tau(STRUCT_H1_ORDER, n, m)[0]
                     for m in METHODS} for n in [50, 200, 1000]},
        "extrap_auc": auc_extrap,
        "win_matrix_cells": wmat.to_dict(),
        "failure_rates": {
            m: float(raw[f"fail_{m}"].notna().mean()) for m in METHODS
            if f"fail_{m}" in raw},
    }
    with open(os.path.join(run_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=1)
    cs.to_csv(os.path.join(run_dir, "cell_stats.csv"), index=False)
    wmat.to_csv(os.path.join(run_dir, "win_matrix.csv"))
    print(json.dumps(summary["h0_false_alarm_alpha05"], indent=1))
    print(json.dumps(summary["power_at_empirical_fpr5_rel1pct"], indent=1))
    return summary, cs, wmat


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("tag", nargs="?", default="pilot")
    ap.add_argument("--shards", default=None,
                    help="逗号分隔分片 tag 列表;给出时合并各分片 raw/extrap 分析")
    a = ap.parse_args()
    main(a.tag, shards=a.shards)
