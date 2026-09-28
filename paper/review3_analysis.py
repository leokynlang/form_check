# -*- coding: utf-8 -*-
"""review3_analysis.py — 第三轮审稿回应的重算脚本(零模拟)。

回应审稿条目(数字回填位置见脚本末尾 NEXT 清单):
  3.2  联合释放的 n 梯度分解:把"危险"定义固定在最大 n 上,再按 n 分层
       ——把(闸门功效)/(分子分母变化)/(最佳 L2 逼近随 n 收缩)三件事分开。
  3.3  危险阈值 sensitivity curve:危险率与三个联合放行率随阈值的曲线。
  3.7  Appendix D 图在 SVD 口径下的重绘(D.1/D.2/D.3/D.4)。
  3.8  Appendix C 里 D 的 SVD 任务级区间与 LOEO;S3 family × method 二项区间;
       Table 6 噪声层在 SVD 路径下的复算(3.8c)。

纪律:
  - 零模拟、零重拟合:只读已核验产物,不训练任何模型。
  - 每个分析块独立 try/except:一块失败不静默,继续跑完其余块,
    退出码非零并逐条列出失败块。
  - 锚点核对:能对上的复现既有报告数字;对不上的打印诊断而非静默通过。
  - 列名自动探测:列不在时跳过该块并明确报告(不猜、不编)。

输出:paper/review3_results.md + figures/review3_fig_D{1..4}.png
"""
import os
import sys
import traceback

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

PAPER_DIR = os.path.dirname(os.path.abspath(__file__))
HERE = os.path.dirname(PAPER_DIR)
# 本脚本可能被复制到仓库之外(例如桌面)运行:此时数据仍从原仓库读,
# 输出(review3_results.md)落在脚本自己旁边。
if not os.path.isdir(os.path.join(HERE, "runs")):
    _fallback = r"D:\form_check"
    if os.path.isdir(os.path.join(_fallback, "runs")):
        HERE = _fallback
FIG_DIR = os.path.join(HERE, "figures")
OUT_MD = os.path.join(PAPER_DIR, "review3_results.md")

STRUCT = ["QUAD_LIN", "CUBIC_LIN", "INTERACT", "THRESHOLD", "EXP_MODEL",
          "LOG_MODEL", "RATIONAL", "OMITTED_M", "DILUTION", "OPPOSING",
          "COLLINEAR_CURV"]
H0_SET = {"H0_BASELINE", "H0_MULTIVAR", "H0_CORR"}
NS = [50, 200, 1000]
THRESHES = [0.02, 0.05, 0.10, 0.20, 0.50]
T_DEFAULT = 0.10
FAILED = []
L = []


def w(s=""):
    L.append(s)


def banner(name):
    print("\n" + "=" * 72 + f"\n{name}\n" + "=" * 72)


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def pct(k, n):
    return f"{100 * k / n:.1f}%" if n else "n/a"


def pick(df, cands):
    """返回第一个存在的列名;都不存在返回 None。"""
    for c in cands:
        if c in df.columns:
            return c
    return None


def load():
    """载入已核验产物。返回 (syn, srsd)。"""
    raw = pd.read_parquet(os.path.join(HERE, "runs", "quality50",
                                       "raw.parquet"))
    ex = pd.read_parquet(os.path.join(HERE, "runs", "quality50",
                                      "extrap.parquet"))
    syn = raw.merge(ex[["scenario", "n", "rel", "seed",
                        "extrap_far_relerr"]],
                    on=["scenario", "n", "rel", "seed"], how="left")
    sr = pd.read_parquet(os.path.join(HERE, "runs", "srsd_full",
                                      "raw_srsd.parquet"))
    sq = pd.read_parquet(os.path.join(HERE, "runs", "srsd_full",
                                      "quality_recompute.parquet"))
    keep = ["task", "n", "rel", "seed"]
    qcols = [c for c in ("r2", "rmse_log", "corr_yhat2_resid")
             if c in sq.columns]
    srsd = sr.merge(sq[keep + qcols], on=keep, how="left")
    cd = pd.read_csv(os.path.join(HERE, "runs", "cd_svd", "s2_cd_svd.csv"),
                     low_memory=False)
    ccols = [c for c in ("rej_C_svd", "rej_D_svd", "fail_C_svd",
                         "fail_D_svd", "p_C_svd", "p_D_svd")
             if c in cd.columns]
    srsd = srsd.merge(cd[keep + ccols], on=keep, how="left")
    print(f"loaded: syn {syn.shape} srsd {srsd.shape}")
    print(f"  cd_svd cols found: {ccols}")
    # 联合沉默(与 no_causal_link.py 同口径:C/D 用 SVD)
    sr_rej = [c for c in ("rej_fc", "rej_fitted_f", "rej_C_svd",
                          "rej_D_svd") if c in srsd.columns]
    sr_fail = [c for c in ("fail_fc", "fail_fitted_f", "fail_C_svd",
                           "fail_D_svd") if c in srsd.columns]
    srsd["_silence"] = (~srsd[sr_rej].fillna(False).any(axis=1)
                        & ~srsd[sr_fail].fillna(False).any(axis=1))
    return syn, srsd


# ────────────────────────────────────────────────── 3.2  n 分解
def analysis_n_decomposition(srsd, t=T_DEFAULT):
    banner("3.2  联合释放的 n 梯度:固定危险定义后的分层")
    base = (srsd.klass == "H1") & (srsd.rel > 0) & srsd.r2.notna() \
        & srsd.extrap_far_relerr.notna()
    w()
    w("### 3.2a 危险定义本身随 n 变化多少(同一 (task, rel, seed) 配对)")
    w()
    piv = srsd[base].pivot_table(index=["task", "rel", "seed"],
                                 columns="n", values="extrap_far_relerr",
                                 aggfunc="first")
    need = [c for c in NS if c in piv.columns]
    if len(need) < 2:
        raise RuntimeError(f"pivot 缺少 n 档:{piv.columns.tolist()}")
    full = piv[need].dropna()
    print(f"配对完整的三元组: {len(full)} / {len(piv)}")
    w(f"配对完整(三档 n 均有域外误差)的 (task, rel, seed): "
      f"{len(full)} / {len(piv)}。")
    w()
    w("| 比较 | 双方都危险 | 仅小 n 危险 | 仅大 n 危险 | 都不危险 |")
    w("|---|---|---|---|---|")
    trans = {}
    for lo, hi in [(50, 200), (200, 1000), (50, 1000)]:
        if lo not in full.columns or hi not in full.columns:
            continue
        a = full[lo] > t
        b = full[hi] > t
        both = int((a & b).sum())
        only_lo = int((a & ~b).sum())
        only_hi = int((~a & b).sum())
        neither = int((~a & ~b).sum())
        trans[(lo, hi)] = (both, only_lo, only_hi, neither)
        w(f"| n={lo} vs n={hi} | {both} | {only_lo} | {only_hi} | {neither} |")
    w()
    w("读法:**仅小 n 危险**一列就是\"估计误差收缩\"直接造成的危险标签翻转"
      "——同一模型、同一真值,只因样本量变大就不再越线;它不经过任何检验。"
      "**仅大 n 危险**一列是反方向(通常极少)。两列之和是危险定义本身对 n"
      "的依赖量,必须从\"检验检出\"的解释里扣除。")
    print(trans)

    w()
    w("### 3.2b 固定危险定义(n=1000 标签)后按 n 分层的联合放行率")
    w()
    if 1000 not in piv.columns:
        raise RuntimeError("缺 n=1000 档,无法固定标签")
    tag = (piv[1000] > t).rename("danger_fixed")
    s = srsd[base].merge(tag, left_on=["task", "rel", "seed"],
                         right_index=True, how="left")
    s = s[s.danger_fixed.notna()]
    w("| n | 危险(固定标签) | 过筛危险 | 联合放行 | 放行/过筛危险 | "
      "放行/全部危险(固定标签) |")
    w("|---|---|---|---|---|---|")
    for n in NS:
        sn = s[s.n == n]
        d = sn[sn.danger_fixed]
        n_d = len(d)
        scr = d[d.r2 >= 0.99]
        n_scr = len(scr)
        rel = int(scr._silence.sum())
        r1 = f"{rel / n_scr:.1%}" if n_scr else "n/a"
        r2 = f"{rel / n_d:.1%}" if n_d else "n/a"
        w(f"| {n} | {n_d} | {n_scr} | {rel} | {r1} | {r2} |")
    w()
    w("这张表才是\"样本量是否关上通道\"的干净读法:危险标签不再随 n 变,"
      "所以剩下的梯度只可能来自筛选与闸门行为。与 Table 15 的原始口径"
      "并读即可分离(2) 分子分母变化与(3) 逼近误差收缩。")


# ────────────────────────────────────────────── 3.3  阈值 sensitivity
def analysis_threshold_sweep(srsd, t=T_DEFAULT):
    banner("3.3  危险阈值 sensitivity curve")
    base = (srsd.klass == "H1") & (srsd.rel > 0) & srsd.r2.notna() \
        & srsd.extrap_far_relerr.notna()
    pool = srsd[base]
    w()
    w("### 3.3a 阈值扫描(同一批模型;每个率都写明分母)")
    w()
    w("| 阈值 | n | 过筛错形式 A | 其中危险 B | 危险/过筛 | 联合放行 C | "
      "C/A | C/B | 全部危险 D | C/D |")
    w("|---|---|---|---|---|---|---|---|---|---|")
    for t_ in THRESHES:
        for n in NS:
            sn = pool[pool.n == n]
            screened = sn[sn.r2 >= 0.99]
            A = len(screened)
            dang = screened[screened.extrap_far_relerr > t_]
            B = len(dang)
            C = int(dang._silence.sum())
            D = int((sn.extrap_far_relerr > t_).sum())
            w(f"| {t_:.0%} | {n} | {A} | {B} | "
              f"{pct(B, A)} | {C} | "
              f"{pct(C, A)} | {pct(C, B)} | {D} | "
              f"{pct(C, D)} |")
    w()
    w(f"论文 Table 15 / §4.10 报的是阈值 = {t:.0%} 那一档的三行;"
      "本表其余行是同一批模型在同一分母定义下的邻域读数。"
      "§4.10 的 headline 对(0.67% / 6.9%)与放行率随阈值的移动幅度,"
      "在回填时按本表逐格引用。")
    w()
    w("### 3.3b 锚点核对(阈值 = 10% 时与 svd_main_results.md §C 对表)")
    w()
    for n in NS:
        sn = pool[pool.n == n]
        scr = sn[sn.r2 >= 0.99]
        B = int((scr.extrap_far_relerr > t).sum())
        C = int(scr[scr.extrap_far_relerr > t]._silence.sum())
        D = int((sn.extrap_far_relerr > t).sum())
        print(f"n={n}: A={len(scr)} B={B} C={C} D={D}")
        w(f"- n = {n}: A = {len(scr)}, B = {B}, C = {C}, D = {D}")
    w()
    w("*期望(svd_main_results.md §C)*:A = 24,060 / 23,932 / 24,000,B = "
      "2,655 / 2,296 / 2,046,C = 350 / 106 / 25,D = 12,369 / 12,348 / "
      "12,046。**若上面数字与期望不符,不要回填差异,先核对分母口径**"
      "(本脚本的 D 用 `klass==H1 且 rel>0 且 extrap 非 NaN`,与报告口径"
      "若差一个筛选条件,差值会全部落在 D 上)。")


# ─────────────────────────────── 3.8a  Appendix C:D 的 SVD 任务级区间
def analysis_task_bootstrap_cd():
    banner("3.8a  Appendix C:D 的 SVD 任务级区间与 LOEO")
    H0 = ["feynman-" + t for t in
          ["i.12.1", "i.29.4", "i.34.27", "i.38.12", "i.39.10", "i.39.22",
           "i.43.31", "i.47.23", "ii.27.16", "ii.27.18", "ii.8.31",
           "iii.12.43", "iii.15.14"]]
    new = pd.read_csv(os.path.join(HERE, "runs", "cd_svd",
                                   "s2_cd_svd.csv"), low_memory=False)
    new = new[new.task.isin(H0) & (new["rel"] > 0)]
    old = pd.concat([pd.read_csv(os.path.join(
        HERE, "runs", f"srsd_full_r{i}", "raw_srsd.csv"), low_memory=False)
        for i in range(5)], ignore_index=True)
    old = old[old.task.isin(H0) & (old.rel > 0)]
    rng = np.random.default_rng(20260920)
    reps = 10000
    tasks = np.array(H0)
    w()
    w("### 3.8a 任务级区间(重抽 13 个 H0 方程,10,000 次)")
    w()
    w("| n | 方法 | 路径 | 池化率 | 任务级 95% | LOEO 范围 |")
    w("|---|---|---|---|---|---|")
    for n in (50, 200, 1000):
        sub = new[new["n"] == n]
        if sub.empty:
            continue
        rates_n = sub.groupby("task").agg(
            c=("rej_C_svd", "mean"), d=("rej_D_svd", "mean"))
        subo = old[old.n == n]
        rates_o = subo.groupby("task").agg(
            c=("rej_fitted_hc3", "mean"), d=("rej_exog_hc3", "mean"))
        weights = sub.groupby("task").size().reindex(tasks).values
        for meth, name in (("c", "C (fitted HC3)"), ("d", "D (exog HC3)")):
            for tag, rates in (("SVD", rates_n), ("pinv", rates_o)):
                vals = rates[meth].reindex(tasks).values
                if np.isnan(vals).any():
                    continue
                bs = np.array([vals[rng.integers(0, len(tasks),
                                                 len(tasks))].mean()
                               for _ in range(reps)])
                lo, hi = np.percentile(bs, [2.5, 97.5])
                loo = [np.average(np.delete(vals, i),
                                  weights=np.delete(weights, i))
                       for i in range(len(tasks))]
                pooled = float(np.average(vals, weights=weights))
                w(f"| {n} | {name} | {tag} | {pooled:.3f} | "
                  f"[{lo:.3f}, {hi:.3f}] | [{min(loo):.3f}, {max(loo):.3f}] |")
                print(f"n={n} {name} {tag}: {pooled:.3f} [{lo:.3f},{hi:.3f}] "
                      f"LOEO [{min(loo):.3f},{max(loo):.3f}]")
    w()
    w("回填目标:§4.1 的 task-level 段落与 Appendix C 的 \"D's SVD-path "
      "task-level interval and leave-one-out range are registered as "
      "outstanding\" 那半句。")


# ─────────────────────────── 3.8b  S3 family × method 二项区间
def analysis_s3_intervals():
    banner("3.8b  S3 family × method 区间")
    sel = pd.read_csv(os.path.join(HERE, "runs", "feynman_families",
                                   "selection_table.csv"))
    fam = dict(zip(sel.eq_name, sel.family))
    klass = dict(zip(sel.eq_name, sel.klass))
    new = pd.read_csv(os.path.join(HERE, "runs", "cd_svd", "s3_cd_svd.csv"),
                      low_memory=False)
    old = pd.read_csv(os.path.join(HERE, "runs", "feynman_families",
                                   "raw_families.csv"), low_memory=False)
    for df in (new, old):
        df["family"] = df.eq_name.map(fam)
        df["klass"] = df.eq_name.map(klass)

    def thr_at_fpr(p, fpr=0.05):
        p = np.sort(np.asarray(p, float))
        p = p[np.isfinite(p)]
        if not len(p):
            return np.nan
        k = max(1, int(round(fpr * len(p))))
        return p[k - 1]

    h0o, h0n = old[old.klass == "H0"], new[new.klass == "H0"]
    h1o, h1n = old[old.klass == "H1"], new[new.klass == "H1"]
    specs = (("A", "p_fc", "p_fc", False),
             ("B", "p_fitted_f", "p_fitted_f", False),
             ("C (SVD)", "p_fitted_hc3", "p_C_svd", True),
             ("D (SVD)", "p_exog_hc3", "p_D_svd", True))
    w()
    w("### 3.8b S3 校准功效 + Wilson 95%(格内独立口径)与方程级 bootstrap")
    w()
    w("| family | method | power | Wilson 95% | 方程级 bootstrap 95% | "
      "方程数 |")
    w("|---|---|---|---|---|---|")
    rng = np.random.default_rng(20260921)
    for f in ("exponential", "rational", "trigonometric"):
        for name, colo, coln, use_new in specs:
            if (colo not in old.columns) or (coln not in new.columns):
                continue
            pool = h0n if use_new else h0o
            col = coln if use_new else colo
            thr = thr_at_fpr(pool[col])
            p1 = (h1n if use_new else h1o)
            p1 = p1[(p1.family == f) & p1[col].notna()]
            pv = p1[col].to_numpy(float)
            hit = (pv < thr).astype(float)
            k, m = int(hit.sum()), len(hit)
            lo, hi = wilson(k, m)
            eqs = p1.eq_name.to_numpy()
            uniq = np.unique(eqs)
            per_eq = np.array([hit[eqs == e].mean() for e in uniq])
            bs = np.array([per_eq[rng.integers(0, len(uniq), len(uniq))].mean()
                           for _ in range(4000)])
            blo, bhi = np.percentile(bs, [2.5, 97.5])
            w(f"| {f} | {name} | {k / m:.3f} | [{lo:.3f}, {hi:.3f}] | "
              f"[{blo:.3f}, {bhi:.3f}] | {len(uniq)} |")
            print(f"{f:14s} {name:8s} {k / m:.3f} "
                  f"[{lo:.3f},{hi:.3f}] boot [{blo:.3f},{bhi:.3f}]")
    w()
    w("回填目标:§4.3 的 \"Family × method binomial intervals are registered "
      "as outstanding\" 与 Table 5 的每格区间。**方程级 bootstrap 一列是"
      "诚实单位**(30 个方程;格内 500 种子不是独立样本)——若两列相差大,"
      "以方程级为准。")


# ────────────────────────────────────── 3.7  Appendix D 重绘(SVD 口径)
def analysis_figures(syn, srsd):
    banner("3.7  Appendix D 重绘(SVD 口径)")
    made = []
    pc = {"A": pick(srsd, ["p_fc", "p_form_check"]),
          "B": pick(srsd, ["p_fitted_f", "p_reset_f"]),
          "C": pick(srsd, ["p_C_svd", "p_fitted_hc3_svd"]),
          "D": pick(srsd, ["p_D_svd", "p_exog_hc3_svd"])}
    print("p 列探测:", pc)

    # D.1 — S1 / S2 的 H0 p 值 CDF(分面板,标注路径)
    try:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
        sc = (pick(syn, ["p_fc", "p_form_check"]),
              pick(syn, ["p_fitted_f", "p_reset_f"]),
              pick(syn, ["p_C_svd", "p_fitted_hc3", "p_reset_hc3"]),
              pick(syn, ["p_D_svd", "p_exog_hc3"]))
        for ax, (df, title, cols, lbl) in zip(
                axes,
                [(syn[syn.scenario.isin(H0_SET)], "S1 synthetic H0 (A/B original path)", sc,
                  ("A form_check", "B RESET fitted F",
                   "C fitted HC3", "D exog HC3")),
                 (srsd[srsd.klass == "H0"], "S2 SRSD-Feynman H0 (C/D SVD)",
                  (pc["A"], pc["B"], pc["C"], pc["D"]),
                  ("A form_check", "B RESET fitted F",
                   "C fitted HC3 (SVD)", "D exog HC3 (SVD)"))]):
            grid = np.linspace(0, 1, 101)
            for col, lab in zip(cols, lbl):
                if col is None or col not in df.columns:
                    continue
                p = df[col].dropna().to_numpy(float)
                if not len(p):
                    continue
                ax.plot(grid, [(p <= g).mean() for g in grid], label=lab)
            ax.plot([0, 1], [0, 1], "k--", lw=0.8)
            ax.set_xlabel("p"); ax.set_ylabel("empirical CDF")
            ax.set_title(title, fontsize=9); ax.legend(fontsize=7)
        fig.suptitle("P-value calibration under strict monomial H0 "
                     "(diagonal = calibrated)", fontsize=10)
        fig.tight_layout()
        pth = os.path.join(FIG_DIR, "review3_fig_D1.png")
        fig.savefig(pth, dpi=170); plt.close(fig)
        made.append(("review3_fig_D1.png", pth))
    except Exception:
        FAILED.append("fig D1"); traceback.print_exc()

    # D.2 — S2 H1 校准后功效热图(真校准版;对照旧图 D.2 的原始率)
    try:
        base = (srsd.klass == "H1") & (srsd.rel > 0) & srsd.r2.notna()
        sub = srsd[base & (srsd.n == 200)]
        tasks = sorted(sub.task.unique())
        methods = [("A", pc["A"]), ("B", pc["B"]),
                   ("C (SVD)", pc["C"]), ("D (SVD)", pc["D"])]
        h0 = srsd[(srsd.klass == "H0") & (srsd.rel > 0) & (srsd.n == 200)]
        mat = np.full((len(tasks), len(methods)), np.nan)
        for j, (nm, col) in enumerate(methods):
            if col is None:
                continue
            t5 = np.nanquantile(h0[col].to_numpy(float), 0.05)
            for i, tk in enumerate(tasks):
                v = sub[sub.task == tk][col].dropna().to_numpy(float)
                if len(v):
                    mat[i, j] = (v < t5).mean()
        fig, ax = plt.subplots(figsize=(5.6, 10))
        im = ax.imshow(mat, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1)
        ax.set_xticks(range(len(methods)))
        ax.set_xticklabels([m[0] for m in methods], rotation=30,
                           ha="right", fontsize=8)
        ax.set_yticks(range(len(tasks)))
        ax.set_yticklabels([t.replace("feynman-", "") for t in tasks],
                           fontsize=6.5)
        for i in range(len(tasks)):
            for j in range(len(methods)):
                if np.isfinite(mat[i, j]):
                    ax.text(j, i, f"{mat[i, j]:.2f}", ha="center",
                            va="center", fontsize=5.5)
        ax.set_title("S2 H1 calibrated power (H0-pool threshold, "
                     "n=200, rel>0)", fontsize=9)
        fig.colorbar(im, ax=ax, fraction=0.04)
        fig.tight_layout()
        pth = os.path.join(FIG_DIR, "review3_fig_D2.png")
        fig.savefig(pth, dpi=170); plt.close(fig)
        made.append(("review3_fig_D2.png", pth))
    except Exception:
        FAILED.append("fig D2"); traceback.print_exc()

    # D.3 — S2 H1 校准功效 vs n(带任务级 bootstrap 区间)
    try:
        fig, ax = plt.subplots(figsize=(7, 4.2))
        rng = np.random.default_rng(20260922)
        for nm, col, mk in (("A form_check", pc["A"], "o"),
                            ("B RESET fitted F", pc["B"], "s"),
                            ("C fitted HC3 (SVD)", pc["C"], "^"),
                            ("D exog HC3 (SVD)", pc["D"], "v")):
            if col is None:
                continue
            xs, ys, los, his = [], [], [], []
            for n in NS:
                h0n = srsd[(srsd.klass == "H0") & (srsd.rel > 0)
                           & (srsd.n == n)]
                subn = srsd[(srsd.klass == "H1") & (srsd.rel > 0)
                            & (srsd.n == n) & srsd.r2.notna()]
                if h0n.empty or subn.empty:
                    continue
                t5 = np.nanquantile(h0n[col].to_numpy(float), 0.05)
                per = subn.groupby("task").apply(
                    lambda g: (g[col].dropna().to_numpy(float) < t5).mean(),
                    include_groups=False).dropna()
                if per.empty:
                    continue
                bs = np.array([per.to_numpy()[rng.integers(0, len(per),
                                                           len(per))].mean()
                               for _ in range(2000)])
                xs.append(n); ys.append(per.mean())
                los.append(np.percentile(bs, 2.5))
                his.append(np.percentile(bs, 97.5))
            if xs:
                ax.errorbar(xs, ys, yerr=[np.array(ys) - los,
                                          np.array(his) - np.array(ys)],
                            marker=mk, capsize=3, label=nm)
        ax.set_xscale("log"); ax.set_xlabel("n")
        ax.set_ylabel("calibrated power (H0-pool threshold)")
        ax.set_title("S2 SRSD-Feynman H1: calibrated power vs n "
                     "(task-level bootstrap 95%)", fontsize=9.5)
        ax.legend(fontsize=7.5); ax.grid(alpha=0.3)
        fig.tight_layout()
        pth = os.path.join(FIG_DIR, "review3_fig_D3.png")
        fig.savefig(pth, dpi=170); plt.close(fig)
        made.append(("review3_fig_D3.png", pth))
    except Exception:
        FAILED.append("fig D3"); traceback.print_exc()

    # D.4 — 报警 vs 外推误差(precision/recall;不含伪影门)
    try:
        t = T_DEFAULT
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
        sub = srsd[(srsd.rel > 0)
                   & srsd.extrap_far_relerr.notna()]
        dang = sub.extrap_far_relerr > t
        names, prec, rec = [], [], []
        for nm, col in (("A", "rej_fc"), ("B", "rej_fitted_f"),
                        ("C (SVD)", "rej_C_svd"), ("D (SVD)", "rej_D_svd")):
            if col not in sub.columns:
                continue
            al = sub[col].fillna(False).astype(bool)
            tp = int((al & dang).sum())
            names.append(nm)
            prec.append(tp / max(1, int(al.sum())))
            rec.append(tp / max(1, int(dang.sum())))
        x = np.arange(len(names))
        axes[0].bar(x - 0.2, prec, 0.4, label="precision")
        axes[0].bar(x + 0.2, rec, 0.4, label="recall")
        axes[0].set_xticks(x); axes[0].set_xticklabels(names)
        axes[0].axhline(dang.mean(), ls="--", lw=0.8,
                        label=f"base rate {dang.mean():.3f}")
        axes[0].set_title(f"alarm precision/recall for extrapolation "
                          f"error > {t:.0%}", fontsize=9)
        axes[0].legend(fontsize=7.5)
        axes[1].hist(sub.extrap_far_relerr, bins=60)
        axes[1].axvline(t, ls="--", lw=0.8)
        axes[1].set_yscale("log")
        axes[1].set_xlabel("far-field relative error")
        axes[1].set_ylabel("count (log)")
        axes[1].set_title("far-field error distribution (S2 tasks)", fontsize=9)
        fig.tight_layout()
        pth = os.path.join(FIG_DIR, "review3_fig_D4.png")
        fig.savefig(pth, dpi=170); plt.close(fig)
        made.append(("review3_fig_D4.png", pth))
    except Exception:
        FAILED.append("fig D4"); traceback.print_exc()

    w()
    w("### 3.7 Appendix D 重绘图(SVD 口径)")
    w()
    if made:
        for nm, pth in made:
            w(f"- `figures/{nm}` — {pth}")
    else:
        w("- 未能生成任何图(见控制台 traceback;通常是 p 值列名探测失败)。")
    w()
    w("回填目标:替换 build_pdf.py 的 `FIG_CAPS` 五图引用为 "
      "review3_fig_D1..D4(如生成成功),并把 Appendix D 的节首"
      "\"archival / 未在 SVD 路径重跑\" 声明改为已重跑。")


# ─────────────────────────── 3.8c  Table 6 噪声层的 SVD 补齐
def analysis_table6_svd():
    banner("3.8c  Table 6 噪声层的 SVD 路径补齐 + 二项区间")
    NOISE = ["H1_HEAVYTAIL", "H1_VAR_GROWING", "H1_AUTOCORR",
             "H1_OUTLIER", "H1_OUTLIER_GRP"]
    LABEL = {"H1_HEAVYTAIL": "Student-t(3) noise",
             "H1_VAR_GROWING": "variance growing with prediction",
             "H1_AUTOCORR": "AR(1) errors, rho = 0.6",
             "H1_OUTLIER": "2% point outliers (+8 sigma)",
             "H1_OUTLIER_GRP": "6% group outliers (+5 sigma)"}
    need = ["rej_fc", "rej_fitted_f", "rej_C_svd", "rej_D_svd"]
    srcs = []
    s1svd = os.path.join(HERE, "runs", "cd_svd", "s1_cd_svd.csv")
    pf = os.path.join(HERE, "runs", "full", "raw.parquet")
    if os.path.exists(pf):
        srcs.append(("runs/full", pf, s1svd))
    srcs.append(("runs/quality50", os.path.join(HERE, "runs", "quality50",
                                                "raw.parquet"),
                 os.path.join(HERE, "runs", "cd_svd", "s1_cd_svd.csv")))
    d, pool = None, None
    for name, path, cdsvd in srcs:
        if not os.path.exists(path):
            continue
        t = pd.read_parquet(path)
        if cdsvd is not None and os.path.exists(cdsvd):
            c1 = pd.read_csv(cdsvd, low_memory=False)
            ks = [k for k in ("scenario", "n", "rel", "seed")
                  if k in c1.columns]
            cc = [c for c in ("rej_C_svd", "rej_D_svd")
                  if c in c1.columns]
            if ks and cc:
                t = t.merge(c1[ks + cc], on=ks, how="left")
                print(f"  merged {cc} from {os.path.basename(cdsvd)}")
        if set(need) <= set(t.columns):
            d, pool = t, name
            break
    if d is None:
        raise RuntimeError(f"找不到同时具备 {need} 的 S1 池")
    sub = d[(d.n == 200) & np.isclose(d.rel, 0.01)
            & d.scenario.isin(NOISE)]
    w()
    w("### 3.8c Table 6(噪声层)在 SVD 路径下的复算 + Wilson 95%")
    w()
    w(f"池:`{pool}`,n = 200,rel = 1%。每格 `率 [Wilson 95%]`。")
    w()
    w("| 噪声场景 | A | B | C (SVD) | D (SVD) | N |")
    w("|---|---|---|---|---|---|")
    for sc in NOISE:
        g = sub[sub.scenario == sc]
        if g.empty:
            w(f"| {LABEL[sc]} | - | - | - | - | 0 |")
            continue
        cells = []
        for c in ("rej_fc", "rej_fitted_f", "rej_C_svd", "rej_D_svd"):
            v = g[c].dropna()
            k, m = int(v.sum()), len(v)
            if m == 0:
                cells.append("-")
                continue
            lo, hi = wilson(k, m)
            cells.append(f"{k / m:.3f} [{lo:.3f},{hi:.3f}]")
        w(f"| {LABEL[sc]} | " + " | ".join(cells) + f" | {len(g)} |")
        print(f"{sc:16s} " + "  ".join(cells) + f"  N={len(g)}")
    w()
    w("回填目标:§4.4 Table 6 的 C/D 两列(论文现值是 pinv 时期行,§6 已标 "
      "outstanding),以及审稿人第 5 节问题 6 的『Table 6 噪声层能否在 SVD "
      "路径补齐 + 二项区间』。")
    w("**两个注意**:①若上面的池是 quality50(200 种子)而非 full(500 "
      "种子),Table 6 题注的『500 seeds』要一并改口径;②审稿人同问的"
      "『S5 的 D 集成』本脚本**不做**——它需要把 D 接进 gate_compare.py "
      "的 harness 重跑,属于重算而非读产物,列为独立待办。")


def main():
    syn, srsd = load()
    w("# review3_results — 第三轮审稿回应的重算结果")
    w()
    w("> 由 `paper/review3_analysis.py` 生成,零模拟,读已核验产物。")
    w("> 数字回填位置见各节末的「回填目标」与脚本末尾 NEXT 清单。")
    for fn, args in ((analysis_n_decomposition, (srsd,)),
                     (analysis_threshold_sweep, (srsd,)),
                     (analysis_task_bootstrap_cd, ()),
                     (analysis_s3_intervals, ()),
                     (analysis_table6_svd, ()),
                     (analysis_figures, (syn, srsd))):
        try:
            fn(*args)
        except Exception:
            FAILED.append(fn.__name__)
            print(f"\n!! {fn.__name__} 失败:")
            traceback.print_exc()
    w()
    w("## NEXT — 回填清单")
    w()
    for s in [
        "§4.10(b): 把「registered as outstanding」换成 3.2a/3.2b 的数字。",
        "§4.10 定义段: 把「registered as outstanding」换成 3.3a 的曲线。",
        "§4.3: 用 3.8b 的区间替换「registered as outstanding」。",
        "Appendix C: 用 3.8a 的 D 的 SVD 区间替换 outstanding 半句。",
        "§4.4 Table 6: 用 3.8c 的 C/D 行替换 pinv 时期行并加区间。",
        "Appendix D: 若 review3_fig_D*.png 生成成功,替换 build_pdf.py 的 "
        "FIG_CAPS 并删掉节首的 archival 声明;否则保留 archival 声明。",
        "§6 Figure provenance / §7: 与上一条保持一致。",
    ]:
        w(f"- {s}")
    w()
    if FAILED:
        w(f"**失败块({len(FAILED)}):** " + ", ".join(FAILED))
    else:
        w("**全部块成功。**")
    with open(OUT_MD, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print(f"\nwrote {OUT_MD}")
    if FAILED:
        print(f"FAILED BLOCKS: {FAILED}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
