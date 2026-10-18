"""gate_compare.py — 关键对照: form_check 门卫 vs "二次族模型比较"门卫。

[背景] 审稿人/自己的核心质疑: "为什么不直接拿一个嵌套的二次族做
CV/AICc 模型比较来决定幂律形式好不好, 而要一个专门的检验?"
本实验正面回答: 三种门卫放在同一台架上比误报/检出/成本。

[预注册设计 — 运行前冻结, 不看结果改判据]
问题设定: workflow 已选中幂律模型(一次 OLS)。门卫的任务 = 判断
"幂律形式对当前数据是否失配"。三个门卫:
  G1 fc       : dimreflex form_check 本体 + 排列经验 p(999), p<0.05 报警。
                成本 = 2 次全量 lstsq + 1 次 pinv 行技巧(排列向量化)。
  G2 cv_quad  : 竞争族 = log 特征的全二次展开(线性+平方+两两交叉, 含截距,
                lstsq 最小范数兜秩亏)。统计量 D = RMSE_cv(power) − RMSE_cv(quad),
                5 折确定性分层(rng=数据种子); D>δ 报警。δ=0 朴素版 +
                δ* 在 H0 池上校准至合并 FPR≤5%(逐 bench×n)。
                成本 = 2 模型 × (5 折 + 1 全量) = 12 次拟合。
  G3 aicc_quad: Δ = AICc(power) − AICc(quad), k=rank(设计矩阵);
                Δ>m 报警。m=0 朴素版 + m* 同上校准。成本 = 2 次拟合。

台架:
  合成 19 场景(scenarios.py, 3 H0 + 11 结构 H1 + 5 噪声层 H1)
  SRSD 全部 32 任务(13 H0 + 19 H1, srsd_data.py 官方采样)
  网格: n ∈ {50, 200}, rel = 0.03, 200 种子。
  种子: 合成 300000+s*97+n_idx*13; SRSD 400000+s*97+n_idx*13(全新区间,
  与 run_experiment/srsd_full 不重叠)。

判读纪律(与正式跑一致):
  - H0: 门卫应安静(误报率)。噪声层 H1: 结构是幂律, 门卫也应安静
    (这层在测"形式门卫是否被噪声假设违规带偏", 不是检出目标)。
  - 结构 H1: 门卫应报警(检出率), 校准阈值(δ*/m*)下的 TPR 为主指标。
  - 成本按拟合次数计, 并列报告。
  - 逐场景/逐方程报告, 不宏平均掩盖盲区。

写入 runs/gate_cmp/ (results.csv + summary.json + report.md)。幂等:
已有 results.csv 时直接进入分析(--rerun 强制重跑)。
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

import common

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "runs", "gate_cmp")
SEEDS = 200
NS = [50, 200]
REL = 0.03
K_FOLDS = 5


# ---------------------------------------------------------------- 拟合小件
def _design(Xlog, quad=False):
    n, p = Xlog.shape
    parts = [np.ones(n), Xlog]
    if quad:
        parts += [Xlog[:, j] ** 2 for j in range(p)]
        parts += [Xlog[:, i] * Xlog[:, j]
                  for i in range(p) for j in range(i + 1, p)]
    return np.column_stack(parts)


def _fit_predict(Xtr, ytr, Xte):
    coef, _, _, _ = np.linalg.lstsq(Xtr, ytr, rcond=None)
    return Xte @ coef


def _aicc(rss, n, k):
    if rss <= 0 or n - k - 1 <= 0:
        return np.nan
    return n * np.log(rss / n) + 2 * k + 2 * k * (k + 1) / (n - k - 1)


# ---------------------------------------------------------------- 两个门卫
def gate_model_compare(Xlog, ly, seed):
    """G2(cv_quad) + G3(aicc_quad)。返回原始统计量与朴素报警, 失败留痕。"""
    n = len(ly)
    out = {"D_cv": np.nan, "d0_fire": False, "delta_aicc": np.nan,
           "m0_fire": False, "fail_mc": None}
    try:
        Xp = _design(Xlog, quad=False)
        Xq = _design(Xlog, quad=True)
        rss_p = float(((ly - _fit_predict(Xp, ly, Xp)) ** 2).sum())
        rss_q = float(((ly - _fit_predict(Xq, ly, Xq)) ** 2).sum())
        kp, kq = int(np.linalg.matrix_rank(Xp)), int(np.linalg.matrix_rank(Xq))
        a_p, a_q = _aicc(rss_p, n, kp), _aicc(rss_q, n, kq)
        if np.isfinite(a_p) and np.isfinite(a_q):
            out["delta_aicc"] = float(a_p - a_q)
            out["m0_fire"] = bool(a_p - a_q > 0)
        # 5 折 CV(确定性, 种子分流)
        rng = np.random.default_rng(int(seed) * 7919 + 13)
        idx = rng.permutation(n)
        folds = np.array_split(idx, K_FOLDS)
        sq_p = sq_q = 0.0
        cnt = 0
        for f in folds:
            mask = np.ones(n, bool)
            mask[f] = False
            tr, te = mask, ~mask
            pred_p = _fit_predict(Xp[tr], ly[tr], Xp[te])
            pred_q = _fit_predict(Xq[tr], ly[tr], Xq[te])
            sq_p += float(((ly[te] - pred_p) ** 2).sum())
            sq_q += float(((ly[te] - pred_q) ** 2).sum())
            cnt += len(f)
        rmse_p, rmse_q = np.sqrt(sq_p / cnt), np.sqrt(sq_q / cnt)
        if np.isfinite(rmse_p) and np.isfinite(rmse_q):
            out["D_cv"] = float(rmse_p - rmse_q)
            out["d0_fire"] = bool(rmse_p - rmse_q > 0)
    except np.linalg.LinAlgError as exc:
        out["fail_mc"] = f"numerical_linalg:{exc}"[:80]
    return out


# ---------------------------------------------------------------- 跑批
def run_synth(rows):
    import scenarios
    for sname, fn in sorted(scenarios.SCENARIOS.items()):
        klass, group, _, _ = scenarios.SCENARIO_META[sname]
        for n_i, n in enumerate(NS):
            for s in range(SEEDS):
                seed = 300000 + s * 97 + n_i * 13
                df, feats = fn(seed, n, REL)
                feats = feats.split(",")
                ly = np.log(df["target"].values)
                Xlog = np.log(df[feats].values)
                row = {"bench": "synth", "task": sname, "klass": klass,
                       "group": group, "n": n, "seed": seed}
                fc = common.form_check_with_p(df, "target", feats, seed)
                row.update({"p_fc": fc["p_fc"], "rej_fc": fc["rej_fc"],
                            "fail_fc": fc["fail_fc"]})
                row.update(gate_model_compare(Xlog, ly, seed))
                rows.append(row)
        print(f"synth {sname} done", flush=True)


def run_srsd(rows):
    import srsd_data
    reg = [t for t in srsd_data.registry_with_meta()
           if t["klass"] in ("H0", "H1")]
    for t in reg:
        fn = srsd_data.make_dataset_fn(t["cls"])
        for n_i, n in enumerate(NS):
            for s in range(SEEDS):
                seed = 400000 + s * 97 + n_i * 13
                row = {"bench": "srsd", "task": t["eq_name"], "klass": t["klass"],
                       "group": "srsd-" + t["klass"], "n": n, "seed": seed}
                try:
                    df, feats = fn(seed, n, REL)
                except Exception as exc:  # noqa: BLE001 — 生成失败留痕
                    row.update({"fail_gen": repr(exc)[:80]})
                    rows.append(row)
                    continue
                feats = feats.split(",")
                ly = np.log(df["target"].values)
                Xlog = np.log(df[feats].values)
                fc = common.form_check_with_p(df, "target", feats, seed)
                row.update({"p_fc": fc["p_fc"], "rej_fc": fc["rej_fc"],
                            "fail_fc": fc["fail_fc"]})
                row.update(gate_model_compare(Xlog, ly, seed))
                rows.append(row)
        print(f"srsd {t['eq_name']} done", flush=True)


# ---------------------------------------------------------------- 分析
def calibrate(h0_stats, alpha=0.05):
    """最小阈值 t 使 H0 池 fire 率 ≤ alpha(统计量越大越报警)。"""
    arr = np.asarray([v for v in h0_stats if np.isfinite(v)])
    if len(arr) == 0:
        return np.inf
    cand = np.unique(np.concatenate([arr, [arr.max() + 1e-12]]))
    for t in cand:
        if float((arr > t).mean()) <= alpha:
            return float(t)
    return float(arr.max() + 1e-12)


def analyze(res):
    lines = ["# 门卫对照: form_check vs 二次族模型比较(CV/AICc)", "",
             f"- 网格: n∈{NS}, rel={REL}, {SEEDS} 种子; 竞争族=log 全二次展开;"
             f" CV={K_FOLDS} 折; 校准=H0 池逐 bench×n 调阈值至 FPR≤5%", ""]
    summary = {"grid": {"n": NS, "rel": REL, "seeds": SEEDS}, "cells": {}}
    ok = res[~res.get("fail_gen", "").notna()] if "fail_gen" in res else res
    for (bench, n), g in ok.groupby(["bench", "n"]):
        h0 = g[g["klass"] == "H0"]
        d_star = calibrate(h0["D_cv"].tolist())
        m_star = calibrate(h0["delta_aicc"].tolist())
        g = g.copy()
        g["cv_fire"] = g["D_cv"] > d_star
        g["aicc_fire"] = g["delta_aicc"] > m_star
        # 原始(朴素)H0 fire 率 — 模型比较不做校准会怎样
        naive = {
            "H0_fire_fc": float(h0["rej_fc"].mean()),
            "H0_fire_cv_d0": float(h0["d0_fire"].mean()),
            "H0_fire_aicc_m0": float(h0["m0_fire"].mean()),
            "H0_fire_cv_cal": float((h0["D_cv"] > d_star).mean()),
            "H0_fire_aicc_cal": float((h0["delta_aicc"] > m_star).mean()),
            "d_star": d_star, "m_star": m_star,
            "n_h0": int(len(h0)),
        }
        # 噪声层 H1(结构=幂律, 门卫应安静): 仅合成
        noise = g[(g["klass"] == "H1") & (g["group"].str.startswith("noise"))]
        struct = g[(g["klass"] == "H1") & (~g["group"].str.startswith("noise"))]
        per_task = []
        for task, gt in struct.groupby("task"):
            per_task.append({
                "task": task, "group": gt["group"].iloc[0], "n": int(n),
                "fc": float(gt["rej_fc"].mean()),
                "cv": float(gt["cv_fire"].mean()),
                "aicc": float(gt["aicc_fire"].mean()),
                "or_fc_cv": float((gt["rej_fc"] | gt["cv_fire"]).mean()),
            })
        dfp = pd.DataFrame(per_task)
        cell = {
            "naive": naive,
            "noise_h1_fire": ({"fc": float(noise["rej_fc"].mean()),
                               "cv": float(noise["cv_fire"].mean()),
                               "aicc": float(noise["aicc_fire"].mean()),
                               "n_cells": int(len(noise))}
                              if len(noise) else None),
            "struct_tpr": {"fc": float(dfp["fc"].mean()),
                           "cv": float(dfp["cv"].mean()),
                           "aicc": float(dfp["aicc"].mean()),
                           "or_fc_cv": float(dfp["or_fc_cv"].mean()),
                           "n_tasks": int(len(dfp))},
            "win_matrix": _wins(dfp),
            "per_task": per_task,
        }
        summary["cells"][f"{bench}_n{n}"] = cell
        lines.append(f"## {bench} n={n}")
        lines.append(f"- H0(n={naive['n_h0']}): fc {naive['H0_fire_fc']:.3f} | "
                     f"cv 朴素δ0 {naive['H0_fire_cv_d0']:.3f} → 校准 "
                     f"{naive['H0_fire_cv_cal']:.3f} | aicc 朴素m0 "
                     f"{naive['H0_fire_aicc_m0']:.3f} → 校准 "
                     f"{naive['H0_fire_aicc_cal']:.3f}")
        if cell["noise_h1_fire"]:
            q = cell["noise_h1_fire"]
            lines.append(f"- 噪声层H1(应安静): fc {q['fc']:.3f} | "
                         f"cv {q['cv']:.3f} | aicc {q['aicc']:.3f}")
        t = cell["struct_tpr"]
        lines.append(f"- 结构H1 校准TPR(均值/任务): fc {t['fc']:.3f} | "
                     f"cv {t['cv']:.3f} | aicc {t['aicc']:.3f} | "
                     f"OR(fc,cv) {t['or_fc_cv']:.3f}  (任务数 {t['n_tasks']})")
        lines.append("- 逐任务: " + ", ".join(
            f"{r.task}[{r.group}]: fc{r.fc:.2f}/cv{r.cv:.2f}/aicc{r.aicc:.2f}"
            for r in dfp.itertuples()))
        lines.append("")
    cost = {"fc": "2 全量拟合 + 1 pinv(999 排列向量化)",
            "cv_quad": "12 次拟合(2 模型×6)", "aicc_quad": "2 次拟合"}
    summary["cost_fits"] = cost
    lines.append("## 成本(拟合次数): " + "; ".join(f"{k}={v}" for k, v in cost.items()))
    with open(os.path.join(OUT, "report.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    with open(os.path.join(OUT, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=1)
    return lines


def _wins(dfp):
    """结构 H1 任务级胜负: 谁报警率高谁胜(平局计平)。"""
    import itertools
    w = {}
    for a, b in itertools.combinations(["fc", "cv", "aicc"], 2):
        wa = int((dfp[a] > dfp[b]).sum())
        wb = int((dfp[b] > dfp[a]).sum())
        tie = int(len(dfp) - wa - wb)
        w[f"{a}_vs_{b}"] = {a: wa, b: wb, "tie": tie}
    return w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rerun", action="store_true")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    csv = os.path.join(OUT, "results.csv")
    if args.rerun or not os.path.exists(csv):
        rows = []
        t0 = time.time()
        run_synth(rows)
        run_srsd(rows)
        res = pd.DataFrame(rows)
        res.to_csv(csv, index=False)
        print(f"run done {time.time()-t0:.0f}s, rows={len(res)}")
    res = pd.read_csv(csv)
    for ln in analyze(res):
        print(ln)
    return 0


if __name__ == "__main__":
    sys.exit(main())
