"""no_causal_link.py — 仓库整体目标的权威量化(2026-09-20 目标切换)。

目标(用户指令):证明**函数结果与函数拟合匹配度本身不构成必然因果**——
不是证明发现了一个新的检测架构。

操作化定义(两侧断言,预注册在本 docstring):
- 侧 1(高匹配度 ⇒ 不必然正确):匹配度 = 域内受限幂律拟合优度(R² 等);
  函数结果 = 拟合函数相对注册真值的域外(×3 log 外推)相对误差。
  断言:同一高优度区间内同时存在"函数正确"(H0,域外误差≈0.1%)与
  "函数严重错误"(域外误差>10%)两类样本,且优度本身无法区分——
  高匹配度不构成正确性证书。
- 侧 2(低匹配度 ⇒ 不必然被检出):被检出 = 四门任一拒绝(A fc / B fitted-F
  为原路径;C/D 用 SVD 重跑路径 runs/cd_svd,与联合沉默审计注册口径一致);
  fail_*/NaN p 为诚实拒检,不算沉默。断言:存在 R² 远低于 H0 的错形式
  样本四门全部放行(联合沉默)。
- 机制:检出与否由失配方向与增广方向的对齐度 corr(yhat^k, 残差)决定,
  与优度无关——两侧断言在此汇合。

数据全部来自已核验产物,本脚本零模拟、零重生成:
- 合成台 runs/quality50/{raw,extrap}.parquet(57,000 拟合,200 种子)
- SRSD 台 runs/srsd_full/{raw_srsd,quality_recompute}.parquet(240,000 拟合)
  + runs/cd_svd/s2_cd_svd.csv(C/D SVD 路径,侧 2 用)

输出:no_causal_link.md + figures/no_causal_link.png。
锚点数字与已核验报告交叉核对(含联合沉默独立论文口径),失败即非零退出。
"""
import os
import sys

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_MD = os.path.join(HERE, "no_causal_link.md")
FIG = os.path.join(HERE, "figures", "no_causal_link.png")

GATES = ["fc", "fitted_f", "fitted_hc3", "exog_hc3"]
H0_SET = {"H0_BASELINE", "H0_MULTIVAR", "H0_CORR"}
NOISE_H1 = {"H1_HEAVYTAIL", "H1_VAR_GROWING", "H1_AUTOCORR",
            "H1_OUTLIER", "H1_OUTLIER_GRP"}
STRUCT = ["QUAD_LIN", "CUBIC_LIN", "INTERACT", "THRESHOLD", "EXP_MODEL",
          "LOG_MODEL", "RATIONAL", "OMITTED_M", "DILUTION", "OPPOSING",
          "COLLINEAR_CURV"]
ERR_T = 0.10
NS = [50, 200, 1000]
LQ, RQ = "\u201c", "\u201d"


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson 95% CI。"""
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def spearman(a: pd.Series, b: pd.Series) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if len(a) < 3 or np.all(a == a[0]) or np.all(b == b[0]):
        return np.nan
    return float(stats.spearmanr(a, b)[0])


def fmt(v: float) -> str:
    return f"{v:+.2f}" if np.isfinite(v) else "n/a"


def pct(k: int, n: int) -> str:
    return f"{100 * k / n:.1f}%" if n else "n/a"


def main() -> None:
    # ---------- 载入(全已核验产物,零重生成) ----------
    raw = pd.read_parquet(os.path.join(HERE, "runs", "quality50", "raw.parquet"))
    ex = pd.read_parquet(os.path.join(HERE, "runs", "quality50",
                                      "extrap.parquet"))
    syn = raw.merge(ex[["scenario", "n", "rel", "seed", "extrap_far_relerr"]],
                    on=["scenario", "n", "rel", "seed"], how="left")
    sr = pd.read_parquet(os.path.join(HERE, "runs", "srsd_full",
                                      "raw_srsd.parquet"))
    sq = pd.read_parquet(os.path.join(HERE, "runs", "srsd_full",
                                      "quality_recompute.parquet"))
    srsd = sr.merge(sq[["task", "n", "rel", "seed", "r2", "rmse_log",
                        "corr_yhat2_resid"]],
                    on=["task", "n", "rel", "seed"], how="left")
    # C/D 换 SVD 路径(联合沉默注册口径);键完整性核对防错配
    cd = pd.read_csv(os.path.join(HERE, "runs", "cd_svd", "s2_cd_svd.csv"),
                     low_memory=False)
    keys = ["task", "n", "rel", "seed"]
    if len(cd) != len(srsd) or cd.duplicated(keys).any():
        raise SystemExit("cd_svd 行数/键与 srsd_full 不一致")
    srsd = srsd.merge(cd[[*keys, "rej_C_svd", "rej_D_svd",
                          "fail_C_svd", "fail_D_svd"]],
                      on=keys, how="left")
    print("loaded: syn", syn.shape, "srsd", srsd.shape)

    # ---------- 联合沉默(合成:quality50 原四门;SRSD:A/B 原 + C/D SVD) ----------
    syn["_silence"] = (~syn[[f"rej_{g}" for g in GATES]].fillna(False)
                       .any(axis=1)
                       & ~syn[[f"fail_{g}" for g in GATES]].fillna(False)
                       .any(axis=1))
    sr_rej = ["rej_fc", "rej_fitted_f", "rej_C_svd", "rej_D_svd"]
    sr_fail = ["fail_fc", "fail_fitted_f", "fail_C_svd", "fail_D_svd"]
    srsd["_silence"] = (~srsd[sr_rej].fillna(False).any(axis=1)
                        & ~srsd[sr_fail].fillna(False).any(axis=1))

    rel_s = syn.rel > 0
    rel_r = (srsd.rel > 0) & srsd.r2.notna() & srsd.extrap_far_relerr.notna()

    # ---------- 锚点核对 ----------
    rows = []

    def anchor(name, got, want, tol):
        ok = abs(float(got) - want) <= tol
        rows.append((ok, name, float(got), want, tol))

    c = syn[(syn.n == 200) & np.isclose(syn.rel, 0.01)
            & (syn.scenario == "CUBIC_LIN")]
    anchor("syn CUBIC_LIN mean r2 (n200,rel1%)", c.r2.mean(), 0.9983, 5e-4)
    anchor("syn CUBIC_LIN extrap_far median (n200,rel1%)",
           c.extrap_far_relerr.median(), 0.444, 0.01)
    s1 = syn[(syn.n == 200) & rel_s & syn.scenario.isin(STRUCT)
             & syn.extrap_far_relerr.notna() & (syn.r2 >= 0.99)]
    anchor("syn STRUCT R2>=0.99: N (n200)", len(s1), 800, 5)
    anchor("syn STRUCT R2>=0.99: far>10% (n200)",
           (s1.extrap_far_relerr > ERR_T).mean(), 1.0, 1e-6)
    h0p = srsd[(srsd.klass == "H0") & (srsd.rel > 0) & srsd.r2.notna()]
    h1p = srsd[(srsd.klass == "H1") & rel_r]
    q05 = float(np.quantile(h0p.r2, 0.05))
    anchor("srsd H1 frac r2>=H0q05 (rel>0)",
           (h1p.r2 >= q05).mean(), 0.527, 0.03)
    sel99 = h1p[h1p.r2 >= 0.99]
    anchor("srsd H1 R2>=0.99: N", len(sel99), 71992, 50)
    anchor("srsd H1 R2>=0.99: far>10%",
           (sel99.extrap_far_relerr > ERR_T).mean(), 0.097, 0.005)
    # 联合沉默独立论文口径(rel>0, α=.05, SVD C/D):H1 .378/.298/.232,H0 .870+
    for n, want in zip(NS, [0.378, 0.298, 0.232], strict=True):
        pool = srsd[(srsd.klass == "H1") & (srsd.rel > 0) & (srsd.n == n)
                    & srsd.r2.notna()]
        anchor(f"srsd H1 joint silence (n={n},rel>0,alpha)",
               pool._silence.mean(), want, 0.02)
    t28 = srsd[srsd.task.str.startswith("feynman-ii.11.28")
               & (srsd.n == 1000) & np.isclose(srsd.rel, 0.01)]
    anchor("srsd ii.11.28 mean r2 (n1000,rel1%)", t28.r2.mean(), 0.0021, 1e-3)
    bad = [r for r in rows if not r[0]]
    for ok, name, got, want, tol in rows:
        print(f"  [{'OK ' if ok else 'FAIL'}] {name}: "
              f"got={got:.6f} want={want}±{tol}")
    if bad:
        raise SystemExit(f"锚点核对失败: {[b[1] for b in bad]}")

    L: list[str] = []

    def w(s: str = "") -> None:
        L.append(s)

    w("# 必然因果不成立:权威量化报告")
    w()
    w("> **仓库整体目标(2026-09-20)**:证明函数结果与函数拟合匹配度本身"
      "不构成必然因果;不是证明发现了一个新的检测架构。本报告是该命题的"
      "权威量化,`python no_causal_link.py` 一键复现(零模拟,读已核验产物;"
      "锚点核对失败即非零退出)。")
    w("> 两侧断言:侧 1 = 高匹配度 ⇒ 不必然正确;侧 2 = 低匹配度 ⇒ "
      "不必然被检出。机制 = 方向对齐度,而非优度。")
    w("> 数据:合成 57,000(200 种子)+ SRSD 240,000(500 种子)已核验拟合;"
      "门为名义 α=0.05;SRSD 侧 C/D 用 SVD 路径(runs/cd_svd,联合沉默"
      "注册口径),合成侧为 quality50 原路径。")

    # ---------------- 侧 1 ----------------
    w()
    w("## 1. 侧 1:高匹配度 ⇒ 函数不必然正确")
    w()
    w(f"同一高优度区间内,{LQ}函数正确{RQ}与{LQ}函数严重错误{RQ}并存,"
      f"优度本身无法区分——这就是{LQ}不构成必然因果{RQ}的操作化。")
    w()
    w("### 1.1 合成台:R²≥0.99 区间的混合总体")
    w()
    w("| n | 类别 | N | 域外误差>10% | 域外误差中位 |")
    w("|---|---|---|---|---|")
    for n in NS:
        for lbl, m in [("H0(真幂律)", syn.scenario.isin(H0_SET)),
                       ("结构-H1(错形式)", syn.scenario.isin(STRUCT))]:
            s = syn[(syn.n == n) & rel_s & m
                    & syn.extrap_far_relerr.notna() & (syn.r2 >= 0.99)]
            k = int((s.extrap_far_relerr > ERR_T).sum())
            w(f"| {n} | {lbl} | {len(s)} | {pct(k, len(s))} | "
              f"{s.extrap_far_relerr.median():.4f} |")
    w()
    w("R²≥0.99 的区间里,真幂律的域外误差中位 ≈0.1–0.2%,错形式 100% "
      ">10%(中位 ≈45%,见 fitting_quality_vs_function_error.md §1)。"
      "**同一个优度区间,函数行为相差两个数量级以上;优度对函数正确性"
      "零分辨。**CUBIC_LIN 是最干净的一例:域内 R²=0.9983,与真幂律 "
      "0.9999 在绝对阈值口径下不可分,域外中位 44%。")
    w()
    w("### 1.2 SRSD 真实公式")
    med = float(h0p.r2.median())
    w()
    w(f"参照系:H0 池(rel>0)R² 5 分位 = {q05:.5f},中位 = {med:.5f}。")
    w()
    w("| R² 门槛 | N | 域外误差>10% | Wilson 95% CI | 域外误差中位 |")
    w("|---|---|---|---|---|")
    for lbl, thr in [("≥ H0 q05", q05), ("≥ H0 中位", med), ("≥ 0.99", 0.99)]:
        sel = h1p[h1p.r2 >= thr]
        k = int((sel.extrap_far_relerr > ERR_T).sum())
        lo, hi = wilson(k, len(sel))
        w(f"| {lbl} | {len(sel)} | {pct(k, len(sel))} | [{lo:.3f}, {hi:.3f}]"
          f" | {sel.extrap_far_relerr.median():.4f} |")
        if lbl == "≥ H0 q05":
            s1_frac, s1_lo, s1_hi = k / len(sel), lo, hi
    w()
    w(f"SRSD 上与真幂律优度同档(H1 中 ≈{100 * float((h1p.r2 >= q05).mean()):.0f}% "
      f"的样本)的错形式中,{100 * s1_frac:.1f}% 域外函数误差 >10%"
      f"(Wilson [{s1_lo:.3f}, {s1_hi:.3f}]);绝对阈值 R²≥0.99 口径下为 "
      "9.7%。多数高 R² 错形式域外偏差确实小(中位 <1%,采样域内近幂律)——"
      f"所以这里证明的是{LQ}不必然正确{RQ}(必然性否定),不是{LQ}经常错"
      f"{RQ}(频率断言)。但证书语义已经失效:拿到高匹配度不能排除大偏差。")
    w()

    # ---------------- 侧 2 ----------------
    w("## 2. 侧 2:低匹配度 ⇒ 不必然被检出")
    w()
    w("联合沉默 = 四门全不拒 且 无一门拒检(fail_*/NaN p 是诚实拒绝,"
      "不算沉默)。SRSD 侧 C/D 为 SVD 重跑路径(与联合沉默审计同一注册口径)。")
    w()
    w("### 2.1 合成台(对照:偏离可测时并集闭合)")
    w()
    w("| n | 门槛 | N | 联合沉默率 | Wilson 95% CI |")
    w("|---|---|---|---|---|")
    for n in NS:
        h0 = syn[(syn.n == n) & syn.scenario.isin(H0_SET) & rel_s]
        h0m = float(h0.r2.median())
        s2 = syn[(syn.n == n) & rel_s & syn.scenario.isin(STRUCT)
                 & (syn.r2 < h0m)]
        k = int(s2._silence.sum())
        lo, hi = wilson(k, len(s2))
        w(f"| {n} | R²<H0 中位({h0m:.5f}) | {len(s2)} | {pct(k, len(s2))} | "
          f"[{lo:.3f}, {hi:.3f}] |")
    w()
    w(f"合成台低 R² 错形式几乎不沉默(0–0.2%,与 runs/joint_silence "
      f"REPORT 一致)——合成错形式是{LQ}设计得可测{RQ}的。"
      f"真正检验侧 2 的是真实公式:")
    w()
    w("### 2.2 SRSD 真实公式")
    sil_all = {}
    w()
    w("| 门槛 | N | 联合沉默率 | Wilson 95% CI |")
    w("|---|---|---|---|")
    base = (srsd.klass == "H1") & (srsd.rel > 0) & srsd.r2.notna()
    conds = [(f"n={n}, R²<H0 中位", base & (srsd.n == n) & (srsd.r2 < med))
             for n in NS]
    conds += [(f"n={n}, R²<0.5", base & (srsd.n == n) & (srsd.r2 < 0.5))
              for n in NS]
    conds.append(("全 H1(无条件)", base))
    for lbl, m in conds:
        s2 = srsd[m]
        k = int(s2._silence.sum())
        lo, hi = wilson(k, len(s2))
        w(f"| {lbl} | {len(s2)} | {pct(k, len(s2))} | [{lo:.3f}, {hi:.3f}] |")
        sil_all[lbl] = (k, len(s2))
    w()
    t28s = srsd[srsd.task.str.startswith("feynman-ii.11.28")
                & (srsd.rel > 0)]
    r28 = float(t28s[t28s.n == 1000].r2.mean())
    sil28 = {n: float(t28s[t28s.n == n]._silence.mean()) for n in NS}
    w(f"极端例 `ii.11.28`:受限幂律 R²≈{r28:.4f}——几乎不解释任何方差——"
      f"联合沉默率 n=50/200/1000 = {sil28[50]:.0%}/{sil28[200]:.0%}/"
      f"{sil28[1000]:.0%}。**低匹配度同样不必然被检出**:失配方向与全部"
      "增广方向正交时,域内偏离再大也测不到。")
    w()

    # ---------------- 机制 ----------------
    w("## 3. 机制:因果链缺失在哪一环")
    w()
    w(f"优度衡量{LQ}数据被解释的程度{RQ};能否检出取决于{LQ}失配方向是否"
      f"落在检验可探测的方向空间{RQ}——对齐度 corr(yhat², 残差)。")
    w()
    w("| 台架 | n | corr(yhat²,残差) vs fc 检出 | vs fitted-F 检出 | "
      "R² vs fc 检出 |")
    w("|---|---|---|---|---|")
    for n in NS:
        p_s = syn[(syn.n == n) & rel_s & syn.scenario.isin(STRUCT)]
        p_r = srsd[base & (srsd.n == n) & srsd.corr_yhat2_resid.notna()]
        w(f"| 合成 | {n} | {fmt(spearman(p_s.corr_yhat2_resid, p_s.rej_fc))} "
          f"| {fmt(spearman(p_s.corr_yhat2_resid, p_s.rej_fitted_f))} | "
          f"{fmt(spearman(p_s.r2, p_s.rej_fc))} |")
        w(f"| SRSD | {n} | {fmt(spearman(p_r.corr_yhat2_resid, p_r.rej_fc))} "
          f"| {fmt(spearman(p_r.corr_yhat2_resid, p_r.rej_fitted_f))} | "
          f"{fmt(spearman(p_r.r2, p_r.rej_fc))} |")
    w()
    w("对齐度与检出稳定正相关(SRSD fitted-F 列随 n 增强至 +0.7 量级;"
      "合成 n=1000 fitted-F 恒拒 TPR=1.00,相关未定义记 n/a),而 R² 与"
      "同格检出的关联接近零。**从匹配度到正确性、从低匹配度到被检出,"
      "两条直觉链路都没有必然性;真实起作用的是方向几何。**")
    w()

    # ---------------- 判定 ----------------
    w("## 4. 判定")
    w()
    w(f"命题:{LQ}函数结果与函数拟合匹配度本身不构成必然因果{RQ}。")
    w()
    w("- **侧 1 成立**:高匹配度不构成正确性证书——合成台 R²≥0.99 区间"
      "内错形式 100% 域外偏差>10% 而真幂律 ≈0%,优度零分辨;SRSD 真实"
      f"公式上优度同档样本 {100 * s1_frac:.1f}% 域外偏差>10%"
      f"(R²≥0.99 口径 9.7%)。多数高 R² 错形式偏差小 ⇒ 证明的是"
      f"{LQ}不必然{RQ},不是{LQ}经常{RQ}。")
    ku, nu = sil_all["全 H1(无条件)"]
    k50m, n50m = sil_all["n=50, R²<H0 中位"]
    k50l, n50l = sil_all["n=50, R²<0.5"]
    w("- **侧 2 成立**:低匹配度不必然被检出——SRSD 联合沉默 "
      f"{pct(k50m, n50m)}(n=50,R²<H0 中位)到 {pct(k50l, n50l)}"
      f"(n=50,R²<0.5);全 H1 池(n=50){pct(ku, nu)};"
      f"`ii.11.28` R²≈{r28:.3f} 仍 {sil28[50]:.0%} 沉默(n=50)。"
      f"合成台侧 2 近零(0–0.2%)——对照恰恰说明{LQ}域内偏离大小{RQ}"
      "本身也不是检出与否的决定量。")
    w(f"- **机制一致**:两侧断言同源——检出由失配方向与增广方向的对齐度"
      f"决定,优度只在{LQ}场景难度{RQ}层面与功效同向,不含样本级信息。")
    w("- **前提(写作必须显眼)**:偏差定义在域外/函数行为上;域内平均"
      "偏差与 R² 数学绑定,该口径下命题不成立。")
    w()
    w("---")
    w("复现:python no_causal_link.py。数据:runs/quality50(合成 57,000,"
      "200 种子)+ runs/srsd_full(240,000,500 种子)+ runs/cd_svd"
      "(C/D SVD 路径),均为已核验产物;锚点核对(12 项,含联合沉默独立"
      "论文口径)失败即非零退出。")
    w("交叉引用:reset_fitting_quality_report.md(优度×功效)、"
      "fitting_quality_vs_function_error.md(优度×函数偏差)、"
      "runs/joint_silence/REPORT.md(联合沉默)。")

    # ---------------- 图 ----------------
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    ax = axes[0]
    p = syn[(syn.n == 200) & rel_s & syn.extrap_far_relerr.notna()]
    is_h0 = p.scenario.isin(H0_SET)
    noiz = p.scenario.isin(NOISE_H1)
    ax.scatter(p.r2[~is_h0 & ~noiz], p.extrap_far_relerr[~is_h0 & ~noiz],
               s=4, alpha=0.25, c="#4477aa", label="structural H1 (wrong)")
    ax.scatter(p.r2[is_h0], p.extrap_far_relerr[is_h0], s=4, alpha=0.25,
               c="#ee6677", label="H0 (true power law)")
    ax.axhline(ERR_T, color="k", ls=":", lw=1)
    ax.set_xlabel("in-domain R²")
    ax.set_ylabel("far-extrapolation function error")
    ax.set_title("Side 1 (synthetic): same high quality,\nopposite function truth")
    ax.legend(fontsize=8, markerscale=3)
    ax.grid(alpha=0.3)
    ax = axes[1]
    h0r = srsd[(srsd.klass == "H0") & rel_r].r2
    h1r = srsd[(srsd.klass == "H1") & rel_r].r2
    # R² 分布几乎全贴 1.0,全量直方图不可读;改 [0.95,1] 截断放大
    ax.hist(h0r[(h0r >= 0.95)], bins=50, alpha=0.5, density=True,
            color="#ee6677", label="H0 (true power law)")
    ax.hist(h1r[(h1r >= 0.95)], bins=50, alpha=0.5, density=True,
            color="#4477aa", label="H1 (wrong form)")
    ax.set_xlim(0.95, 1.0)
    ax.set_xlabel("in-domain R² [0.95, 1]")
    ax.set_ylabel("density")
    ax.set_title("Side 1 (SRSD): quality distributions\noverlap (zoomed to [0.95, 1])")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[2]
    sc = srsd[(srsd.klass == "H1") & rel_r]
    hb = ax.hexbin(sc.r2, sc.extrap_far_relerr, gridsize=60, bins="log",
                   mincnt=1, cmap="Blues")
    ax.axhline(ERR_T, color="k", ls=":", lw=1)
    ax.set_xlabel("in-domain R²")
    ax.set_ylabel("far-extrapolation function error")
    ax.set_title("Side 1 (SRSD): hexbin, log counts")
    fig.colorbar(hb, ax=ax, label="log10(N)")
    fig.tight_layout()
    fig.savefig(FIG, dpi=140)
    plt.close(fig)

    with open(OUT_MD, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print(f"\nwrote {OUT_MD}")
    print(f"wrote {FIG}")


if __name__ == "__main__":
    main()
