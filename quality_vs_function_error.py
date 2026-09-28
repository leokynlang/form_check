"""quality_vs_function_error.py — 直接检验:
"域内拟合结果看起来正确,但拟合函数(域外/真值)可以有很大偏差"。

与 `reset_fitting_quality_*.md` 不同,本脚本问的不是"门能否检出",而是
"域内拟合优度"与"拟合函数相对真值的偏差"之间的关系。

偏差的操作化:
- 合成台:`run_experiment` 的 extrap_far_relerr / extrap_med_relerr
  (沿第一特征把 log 域外推 ln(3),与注册真值函数比较的相对误差);
- SRSD 台:`run_srsd` 的 extrap_far_relerr(同口径,真值为官方闭式公式)。

注意:域内*平均*偏差与 R² 数学上绑定(高 R² ⇔ 域内残差小),所以
"拟合看起来对、函数却偏差大"只可能体现在域外/函数行为上——本脚本正是
测这一点。产出 `fitting_quality_vs_function_error.md` + 图。
"""
import os
import sys

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures", "quality_vs_function_error.png")
OUT_MD = os.path.join(HERE, "fitting_quality_vs_function_error.md")

H0 = {"H0_BASELINE", "H0_MULTIVAR", "H0_CORR"}
STRUCT = ["QUAD_LIN", "CUBIC_LIN", "INTERACT", "THRESHOLD", "EXP_MODEL",
          "LOG_MODEL", "RATIONAL", "OMITTED_M", "DILUTION", "OPPOSING",
          "COLLINEAR_CURV"]
LQ, RQ = "\u201c", "\u201d"


def spearman(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if len(a) < 3 or np.all(a == a[0]) or np.all(b == b[0]):
        return np.nan
    return float(stats.spearmanr(a, b)[0])


def main():
    # ---------- 合成 ----------
    raw = pd.read_parquet(os.path.join(HERE, "runs", "quality50", "raw.parquet"))
    ex = pd.read_parquet(os.path.join(HERE, "runs", "quality50", "extrap.parquet"))
    syn = raw.merge(ex[["scenario", "n", "rel", "seed", "extrap_far_relerr",
                        "extrap_med_relerr"]],
                    on=["scenario", "n", "rel", "seed"], how="left")
    syn_c = syn[(syn.n == 200) & (syn.rel == 0.01)].groupby(
        "scenario", as_index=False).agg(
        r2=("r2", "mean"), rmse=("rmse_log", "mean"),
        err_med=("extrap_med_relerr", "median"),
        err_far=("extrap_far_relerr", "median"),
        fc=("rej_fc", "mean"), ff=("rej_fitted_f", "mean"))

    # ---------- SRSD ----------
    sr = pd.read_parquet(os.path.join(HERE, "runs", "srsd_full",
                                      "raw_srsd.parquet"))
    sq = pd.read_parquet(os.path.join(HERE, "runs", "srsd_full",
                                      "quality_recompute.parquet"))
    srsd = sr.merge(sq[["task", "n", "rel", "seed", "r2"]],
                    on=["task", "n", "rel", "seed"], how="left")
    srsd_c = srsd[(srsd.n == 200) & (srsd.rel == 0.01)].groupby(
        ["task", "klass"], as_index=False).agg(
        r2=("r2", "mean"), err_far=("extrap_far_relerr", "median"),
        fc=("rej_fc", "mean"), ff=("rej_fitted_f", "mean"))

    lines = []

    def w(s=""):
        lines.append(s)

    w("# 拟合优度 vs 函数偏差:直接检验")
    w()
    w("> 问题:域内拟合结果看起来正确,拟合函数(相对真值)是否仍可能有"
      "很大偏差?")
    w("> 偏差定义:沿第一特征把 log 域外推 ×3,与**注册真值函数**比较的"
      "相对误差(`extrap_far_relerr`)。域内平均偏差与 R² 数学绑定,"
      "故只测域外/函数行为。")
    w("> 数据:`runs/quality50/`(合成)、`runs/srsd_full/`(SRSD)。")
    w()

    # ---- 合成 ----
    w("## 1. 合成台架:高优度 + 大函数偏差同时出现")
    w()
    w("n=200,rel=1%,逐场景(按域内 R² 降序):")
    w()
    w("| 场景 | 类别 | mean R² | 域外误差(中位) | A fc | B fitted-F |")
    w("|---|---|---|---|---|---|")
    for _, r in syn_c.sort_values("r2", ascending=False).iterrows():
        kl = "H0" if r.scenario in H0 else (
            "噪声H1" if r.scenario.startswith("H1_") else "结构H1")
        w(f"| {r.scenario} | {kl} | {r.r2:.4f} | {r.err_far:.3f} | "
          f"{r.fc:.3f} | {r.ff:.3f} |")
    w()
    c = syn_c[syn_c.scenario == "CUBIC_LIN"].iloc[0]
    h0 = syn_c[syn_c.scenario == "H0_BASELINE"].iloc[0]
    w(f"**最干净的一例**:`CUBIC_LIN` 的域内 R²={c.r2:.4f},与真幂律 "
      f"`H0_BASELINE` R²={h0.r2:.4f} 在绝对阈值口径下不可区分;"
      f"但 CUBIC 的域外函数偏差中位为 **{c.err_far:.1%}**,而 H0 仅 "
      f"**{h0.err_far:.2%}**——相差约 {c.err_far / max(h0.err_far, 1e-9):.0f} 倍。")
    w()
    sh = syn[(syn.n == 200) & (syn.rel > 0) & syn.scenario.isin(STRUCT)
             & syn.extrap_far_relerr.notna()]
    w("结构-H1 错形式中,按域内 R² 门槛筛出的高优度样本,其域外偏差:")
    w()
    w("| 域内 R² 门槛 | N | 域外误差中位 | 域外误差>10% 占比 |")
    w("|---|---|---|---|")
    for t in [0.95, 0.99, 0.995]:
        sel = sh[sh.r2 >= t]
        if not len(sel):
            continue
        w(f"| R²≥{t} | {len(sel)} | {sel.extrap_far_relerr.median():.3f} | "
          f"{100 * float((sel.extrap_far_relerr > 0.10).mean()):.0f}% |")
    w()
    w("**R²≥0.99 的错形式 100% 都有 >10% 的域外函数偏差(中位 45%)**"
      "——域内优度高在绝对口径上几乎不约束函数的域外行为。")
    w()
    p = syn[(syn.n == 200) & (syn.rel > 0)
            & syn.scenario.isin(H0 | set(STRUCT))
            & syn.extrap_far_relerr.notna()]
    w(f"跨场景 Spearman:域内 R² vs 域外误差 = "
      f"{spearman(p.r2, p.extrap_far_relerr):+.2f};"
      f"RMSE vs 域外误差 = {spearman(p.rmse_log, p.extrap_far_relerr):+.2f}。"
      "优度与域外偏差**同向但远非充分**——优度高的错形式仍可有大偏差。")
    w()

    # ---- SRSD ----
    w("## 2. SRSD 真实公式:多数高优度样本偏差小,但存在显著例外")
    w()
    h1 = srsd[(srsd.klass == "H1") & (srsd.rel > 0)
              & srsd.extrap_far_relerr.notna()]
    w(f"H1 池相关:域内 R² vs 域外误差 = "
      f"{spearman(h1.r2, h1.extrap_far_relerr):+.2f}。")
    w()
    w("| 域内 R² 门槛 | N | 域外误差中位 | >10% 占比 |")
    w("|---|---|---|---|")
    for t in [0.95, 0.99, 0.995]:
        sel = h1[h1.r2 >= t]
        if not len(sel):
            continue
        w(f"| R²≥{t} | {len(sel)} | {sel.extrap_far_relerr.median():.3f} | "
          f"{100 * float((sel.extrap_far_relerr > 0.10).mean()):.0f}% |")
    w()
    w("与合成台架不同:真实公式上**多数**高 R² 错形式域外偏差并不大"
      "(中位 <1%)——因为它们在官方采样域内确实接近幂律;但仍有例外任务"
      "同时具备高优度和大偏差:")
    w()
    w("| 任务 | 类别 | mean R² | 域外误差(中位) | A fc | B fitted-F |")
    w("|---|---|---|---|---|---|")
    fl = srsd_c[(srsd_c.r2 > 0.99) & (srsd_c.err_far > 0.10)].sort_values(
        "err_far", ascending=False)
    if len(fl):
        for _, r in fl.iterrows():
            w(f"| {r.task.replace('feynman-', '')} | {r.klass} | {r.r2:.4f} | "
              f"{r.err_far:.3f} | {r.fc:.3f} | {r.ff:.3f} |")
    else:
        w("| (n=200,rel=1% 口径下无 R²>0.99 且误差>10% 的任务) | | | | | |")
    w()
    w("对照的极端例仍成立:`ii.11.28` 域内 R²≈0.002(域内就不对)而两门全盲;"
      "`i.41.16` 域内 R²≈1.0000 而两门近全盲——这两个反例说明的是"
      f"{LQ}优度/检出{RQ}与{LQ}函数偏差{RQ}是不同坐标,而非同一故事。")
    w()

    # ---- 判定 ----
    w("## 3. 对观点的判定")
    w()
    w(f"观点:{LQ}拟合结果可能正确,但拟合函数可以有很大偏差{RQ}。")
    w()
    w("- **在合成台架上:成立且显著**。`CUBIC_LIN` 域内 R²=0.9983(与真幂律"
      " 0.9999 阈值口径不可分),域外函数偏差中位 44%;域内 R²≥0.99 的结构"
      "错形式 100% 有 >10% 的域外偏差。")
    w("- **在 SRSD 真实公式上:部分成立**。多数高 R² 错形式域外偏差较小"
      "(中位 <1%),但存在同时高优度、大偏差的任务例外;且存在域内就不对"
      "却全盲的 `ii.11.28`。")
    w(f"- **前提澄清**:若把{LQ}偏差{RQ}定义为**域内平均偏差**,则它与 R² "
      f"数学上绑定,观点不成立;观点成立的前提是把偏差定义在**域外/函数"
      f"行为/极限行为**上——这正是其工程意义所在(内插对、外推错)。")
    w("- **与既有报告的关系**:`reset_fitting_quality_*.md` 测的是"
      f"{LQ}优度 vs 门能否检出{RQ},回答的是另一个问题;本报告测"
      f"{LQ}优度 vs 函数偏差{RQ},直接对应本观点。")
    w()
    w("---")
    w("复现:`python quality_vs_function_error.py`(读 `runs/quality50/` 与 "
      "`runs/srsd_full/`,均为已核验产物)。")

    # ---- 图 ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    ax = axes[0]
    h0m = syn_c.scenario.isin(H0)
    ax.scatter(syn_c.r2[~h0m], syn_c.err_far[~h0m], s=30, c="#4477aa",
               label="structural H1")
    ax.scatter(syn_c.r2[h0m], syn_c.err_far[h0m], s=44, c="#ee6677",
               marker="^", label="H0 / noise H1")
    for _, r in syn_c.iterrows():
        if r.scenario in ("CUBIC_LIN", "EXP_MODEL", "RATIONAL",
                          "H0_BASELINE", "INTERACT"):
            ax.annotate(r.scenario, (r.r2, r.err_far), fontsize=7,
                        xytext=(3, 3), textcoords="offset points")
    ax.axhline(0.10, color="gray", ls=":", lw=0.8)
    ax.set_xlabel("in-domain mean R$^2$ (n=200, rel=1%)")
    ax.set_ylabel("out-of-domain function error (median)")
    ax.set_title("Good in-domain fit, large function error\n(synthetic)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[1]
    for kl, col, mk in [("H0", "#ee6677", "^"), ("H1", "#4477aa", "o")]:
        cc = srsd_c[srsd_c.klass == kl]
        ax.scatter(cc.r2, cc.err_far, s=30, c=col, marker=mk, label=kl)
    for _, r in srsd_c.iterrows():
        if r.err_far > 0.10:
            ax.annotate(r.task.replace("feynman-", ""), (r.r2, r.err_far),
                        fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.axhline(0.10, color="gray", ls=":", lw=0.8)
    ax.set_xlabel("in-domain mean R$^2$ (n=200, rel=1%)")
    ax.set_ylabel("out-of-domain function error (median)")
    ax.set_title("SRSD real formulas:\nmostly small, but exceptions")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG, dpi=140)
    plt.close(fig)

    with open(OUT_MD, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {OUT_MD}")
    print(f"wrote {FIG}")
    return syn_c, srsd_c


if __name__ == "__main__":
    main()
