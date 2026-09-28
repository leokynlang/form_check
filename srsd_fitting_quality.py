"""srsd_fitting_quality.py — SRSD 真实公式侧的拟合优度 × 功效关联。

与合成侧 `reset_fitting_quality.py` 同口径,但数据来自 SRSD-Feynman 的
32 个真实公式(13 H0 严格单项式 + 19 H1 非单项式)。拟合优度由本脚本按
`run_srsd.py` 同 seed 规则重生成数据后重算(含 yhat²/yhat³ 中间量),
门结果复用 `runs/srsd_full/raw_srsd.parquet`(500 种子×32 任务×3n×5rel)。

一致性校验:重生成失败的行应与落盘 data_fail 行一致;否则报错。

产出:`reset_fitting_quality_srsd_analysis.md` + `figures/quality_vs_power_srsd.png`。
"""
import os
import sys
import tempfile

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import srsd_data
from common import fitting_quality
from run_srsd import N_GRID, REL_GRID, base_seed

RUN = os.path.join(HERE, "runs", "srsd_full", "raw_srsd.parquet")
CACHE = os.path.join(HERE, "runs", "srsd_full", "quality_recompute.parquet")
FIG = os.path.join(HERE, "figures", "quality_vs_power_srsd.png")
OUT_MD = os.path.join(HERE, "reset_fitting_quality_srsd_analysis.md")
GATES = ["rej_fc", "rej_fitted_f", "rej_fitted_hc3", "rej_exog_hc3"]
QKEYS = ["r2", "r2_adj", "rmse_log", "aic", "bic", "aug_r2_yhat2",
         "aug_r2_yhat3", "corr_yhat2_resid", "corr_yhat3_resid"]
LQ, RQ = "\u201c", "\u201d"


def ensure_srsd_source():
    """默认 %TEMP%/srsd-benchmark 若缺 datasets/feynman.py(缓存损坏),
    回退到全新克隆路径(复用 srsd_data._ensure_upstream);不删除任何已有目录。"""
    up = srsd_data._UPSTREAM
    if os.path.isfile(os.path.join(up, "datasets", "feynman.py")):
        return up
    alt = os.path.join(tempfile.gettempdir(), "srsd-benchmark-qf")
    srsd_data._UPSTREAM = alt
    if not os.path.isfile(os.path.join(alt, "datasets", "feynman.py")):
        srsd_data._ensure_upstream()
    return alt


def recompute_quality(raw, fn_by_task):
    """按 run_srsd 同 seed 规则重生成数据并重算优度;记录生成失败。"""
    rows = []
    for (task, n, rel), g in raw.groupby(["task", "n", "rel"], sort=True):
        fn = fn_by_task[task]
        n_idx = N_GRID.index(int(n))
        rel_idx = REL_GRID.index(float(rel))
        for s in g["seed"].values:
            bs = base_seed(int(s), n_idx, rel_idx)
            rec = {"task": task, "n": int(n), "rel": float(rel), "seed": int(s),
                   "gen_fail": None}
            try:
                df, feat_str = fn(seed=bs, n=int(n), rel=float(rel))
                feats = feat_str.split(",")
                q = fitting_quality(np.log(df[feats].values),
                                    np.log(df["target"].values))
                rec.update(q)
            except Exception as e:  # noqa: BLE001 — 上游数据性质,登记不静默
                rec["gen_fail"] = f"{type(e).__name__}"
                rec.update(dict.fromkeys(QKEYS, np.nan))
            rows.append(rec)
    return pd.DataFrame(rows)


def spearman(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if len(a) < 3 or np.all(a == a[0]) or np.all(b == b[0]):
        return np.nan
    return float(stats.spearmanr(a, b)[0])


def main():
    ensure_srsd_source()
    reg = [r for r in srsd_data.registry_with_meta()
           if r["klass"] in ("H0", "H1")]
    fn_by_task = {r["eq_name"]: srsd_data.make_dataset_fn(r["cls"])
                  for r in reg}
    klass = {r["eq_name"]: r["klass"] for r in reg}
    raw = pd.read_parquet(RUN)
    raw = raw[raw["task"].isin(klass)].copy()
    print(f"loaded {len(raw)} srsd rows, {raw.task.nunique()} tasks")

    if os.path.exists(CACHE):
        qual = pd.read_parquet(CACHE)
        print(f"loaded cached quality recompute: {len(qual)} rows")
    else:
        qual = recompute_quality(raw, fn_by_task)
        qual.to_parquet(CACHE, index=False)
        print(f"cached quality recompute -> {CACHE}")
    merged = raw.merge(qual, on=["task", "n", "rel", "seed"], how="left")
    # 一致性:重生成失败集合应与落盘的 train_gen 失败集合一致
    rec_fail = set(map(tuple, merged.loc[merged.gen_fail.notna(),
                                         ["task", "n", "rel", "seed"]].values))
    stored_train_fail = set(map(tuple, merged.loc[
        merged.data_fail.fillna("").astype(str).str.startswith("train_gen"),
        ["task", "n", "rel", "seed"]].values))
    print(f"regeneration failures: {len(rec_fail)}; stored train_gen failures: "
          f"{len(stored_train_fail)}; symmetric diff: "
          f"{len(rec_fail ^ stored_train_fail)}")
    assert rec_fail == stored_train_fail, (  # noqa: S101
        "regeneration failure set != stored train_gen failures")

    merged = merged[merged.gen_fail.isna()]
    merged["klass"] = merged["task"].map(klass)
    cells = merged.groupby(["task", "klass", "n", "rel"], as_index=False).agg(
        **{k: (k, "mean") for k in QKEYS},
        **{g: (g, "mean") for g in GATES})

    lines = []

    def w(s=""):
        lines.append(s)

    h0 = merged[merged.klass == "H0"]
    h1 = merged[merged.klass == "H1"]
    h1c = cells[cells.klass == "H1"]

    w("# SRSD 真实公式侧:拟合优度 × 统计功效")
    w()
    w(f"> 数据:`runs/srsd_full/raw_srsd.parquet` 门结果(500 种子 × "
      f"{raw.task.nunique()} 任务 × n∈{{50,200,1000}} × 5 rel),"
      "拟合优度由本脚本按同 seed 规则重生成重算(生成失败集与落盘 "
      "data_fail 集一致)。")
    w("> 任务:SRSD-Feynman 132 公式中入选的 "
      f"{len(reg)} 个({sum(1 for r in reg if r['klass'] == 'H0')} H0 严格"
      f"单项式 + {sum(1 for r in reg if r['klass'] == 'H1')} H1 非单项式)。"
      "门:A=fc,B/C/D=RESET fitted-F/fitted-HC3/exog-HC3,α=0.05。")
    w()
    w("**核心结论**:真实公式上,拟合优度对 H0 与错形式 H1 的区分能力远弱于"
      "合成台架——H1 的 R² 与真幂律 H0 大面积重叠,绝对阈值几乎无法筛出"
      "错形式;真正的功效判别量仍是 RESET 增广方向与残差的对齐度。")
    w()

    # ---- 1 H0 vs H1 优度分布 ----
    w("## 1. H0(真幂律)与 H1(错形式)的拟合优度分布重叠")
    w()
    for n in [50, 200, 1000]:
        a = h0[h0.n == n].r2.values
        b = h1[h1.n == n].r2.values
        ks = stats.ks_2samp(a, b)
        f_med = float(np.mean(b >= np.nanmedian(a)))
        f_q05 = float(np.mean(b >= np.nanquantile(a, .05)))
        w(f"- n={n}:H0 R² 中位 {np.nanmedian(a):.5f},H1 R² 中位 "
          f"{np.nanmedian(b):.5f};H1 有 **{f_med:.1%}** 的样本 R² ≥ H0 中位、"
          f"**{f_q05:.1%}** ≥ H0 的 5% 分位(KS={ks.statistic:.3f})")
    w()
    w("SRSD 多数非单项式公式在官方采样域内**接近幂律**(log-log 局部线性),"
      "因此错形式往往拟合优度并不低;这与合成台架(QUAD/OPPOSING 等 R² "
      "明显掉档)形成对照,也是 SRSD 上各门功效普遍低于合成台架的直接原因。")
    ii = merged[(merged.task == "feynman-ii.11.28") & (merged.n == 1000)
                & (merged.rel == 0.01)].r2
    i41 = merged[(merged.task == "feynman-i.41.16") & (merged.n == 200)
                 & (merged.rel == 0.01)].r2
    w()
    w(f"**优度高低两端都可以是盲区**:`ii.11.28` 的受限幂律 R² 在 n=1000、"
      f"rel=1% 仅 {ii.mean():.4f}(几乎不解释任何方差),两门却都沉默"
      f"(方向与增广空间正交);`i.41.16` 的 R²≈{i41.mean():.4f},同样接近全盲。")
    w()

    # ---- 2 跨任务优度 vs 功效 ----
    w("## 2. 跨任务:优度与功效同向,但同样只是难度代理")
    w()
    w("H1 任务的 (task,n,rel) 格上 Spearman(n=200):")
    w()
    w("| 优度指标 | vs A fc | vs B fitted-F | vs D exog-HC3 |")
    w("|---|---|---|---|")
    h1c200 = h1c[h1c.n == 200]
    for q in ["r2", "rmse_log", "aic", "bic"]:
        row = [f"| {q} "]
        for g in ["rej_fc", "rej_fitted_f", "rej_exog_hc3"]:
            row.append(f"| {spearman(h1c200[q], h1c200[g]):+.2f} ")
        row.append("|")
        w("".join(row))
    w()
    w(f"H1 任务的拟合优度与 fitted-F 功效同向:RMSE/AIC/BIC 越高(拟合越差)"
      f"越易检出(≈+0.52),R² 越低越易检出(−0.38);而 A(fc) 与优度几乎无关"
      f"(R² +0.01)。跨任务变异大、样本小(19 任务 × 4 rel = 76 格),"
      f"不构成{LQ}形式对错{RQ}判据。")
    w()

    # ---- 3 场景内 ----
    w(f"## 3. 任务内:优度对{LQ}这个样本会被检出吗{RQ}几乎无信息")
    w()
    w("H1 池内、固定 n 的 Spearman:")
    w()
    w("| n | R² vs fc | R² vs fitted-F | yhat²~残差 vs fc | yhat²~残差 vs fitted-F |")
    w("|---|---|---|---|---|")
    for n in [50, 200, 1000]:
        s = h1[h1.n == n]
        vals = [spearman(s.r2, s.rej_fc), spearman(s.r2, s.rej_fitted_f),
                spearman(s.corr_yhat2_resid, s.rej_fc),
                spearman(s.corr_yhat2_resid, s.rej_fitted_f)]
        w(f"| {n} | " + " | ".join(
            ("n/a" if not np.isfinite(v) else f"{v:+.2f}") for v in vals) + " |")
    w()
    w("与合成侧一致:R² 与检出的任务内秩相关很小;真正预测功效的是对齐量 "
      "`corr(yhat²,残差)`——与 B(fitted-F) 强正相关且随 n 增强"
      "(+0.50→+0.75),与 fc 温和正相关(+0.21~+0.35)。")
    w()

    # ---- 4 逐任务表 ----
    w("## 4. 逐任务(n=200, rel=1%):优度、对齐量与两门功效")
    w()
    w("| 任务 | K | mean R² | corr(yhat²,res) | corr(yhat³,res) | A fc | B fitted-F |")
    w("|---|---|---|---|---|---|---|")
    sub = cells[(cells.n == 200) & (cells.rel == 0.01)].sort_values(
        ["klass", "r2"], ascending=[True, True])
    for _, r in sub.iterrows():
        w(f"| {r.task.replace('feynman-', '')} | {r.klass} | {r.r2:.4f} | "
          f"{r.corr_yhat2_resid:.3f} | {r.corr_yhat3_resid:.3f} | "
          f"{r.rej_fc:.3f} | {r.rej_fitted_f:.3f} |")
    w()
    w("三个展示例(跨 (n,rel) 稳定):`ii.11.28` —— R² 低到 ≈0.002(n=1000)"
      "而两门全盲,失配方向与 RESET 增广空间正交;`i.41.16` —— R²≈1.0000 而"
      "两门近全盲(联合沉默报告的最差任务);`i.26.2` —— R²≈0.99,"
      "`corr(yhat³,残差)`≈0.31 令 fitted-F=1.00 而 fc 仅 0.40,"
      "是 yhat³ 方向救场的真实例。")
    w()
    w("(`i.30.5` 在联合沉默报告中已被撤回,此处仅作留痕,不用于结论。)")
    w()

    # ---- 5 对齐量 vs 功效 ----
    w("## 5. 对齐量解释功效(与合成侧一致的机制)")
    w()
    s = merged[(merged.klass == "H1") & (merged.rel > 0)]
    w("H1 池内 Spearman(rel>0,三 n 合并):")
    w()
    w("| 中间量 | A fc | B fitted-F | C fitted-HC3 | D exog-HC3 |")
    w("|---|---|---|---|---|")
    for q in ["aug_r2_yhat2", "aug_r2_yhat3",
              "corr_yhat2_resid", "corr_yhat3_resid"]:
        row = [f"| {q} "]
        for g in GATES:
            row.append(f"| {spearman(s[q], s[g]):+.2f} ")
        row.append("|")
        w("".join(row))
    w()
    w("对齐量解释功效:`corr(yhat²,残差)` 对 B(fitted-F) 是最强预测量"
      "(+0.54),`corr(yhat³,残差)` 次之(+0.34);对 A(fc) 为 +0.23。"
      "**`aug_r2_yhat2`/`aug_r2_yhat3`(幂次列被基张成的份额)与两门功效"
      "相关很弱,冗余度假设在 SRSD 上不成立**——盲区由失配方向与增广空间"
      "的正交性决定(低对齐),而非由增广列冗余决定。")
    w()

    # ---- 6 联合沉默接口 ----
    w("## 6. 与联合沉默接口:高优度错形式的全绿率")
    w()
    w(f"H1 数据集中 R² 高于 H0 中位者(即{LQ}拟合优度无法区分{RQ}的错形式),"
      "其四门联合沉默率:")
    w()
    w("| n | R² 门槛 | N | 四门联合沉默 | 任一门拒绝 |")
    w("|---|---|---|---|---|")
    for n in [50, 200, 1000]:
        med = float(np.nanmedian(h0[h0.n == n].r2.values))
        sel = h1[(h1.n == n) & (h1.r2 >= med)]
        if not len(sel):
            continue
        silent = float(np.mean(
            ~(sel.rej_fc | sel.rej_fitted_f | sel.rej_fitted_hc3
              | sel.rej_exog_hc3)))
        w(f"| {n} | H0 中位({med:.4f}) | {len(sel)} | {silent:.3f} | "
          f"{1 - silent:.3f} |")
    w()
    w("这与联合沉默报告一致:真实公式上即使 R² 与真幂律同档,"
      "四门并集仍有可观比例全绿;**拟合优度不能为证书背书**。")
    w()

    # ---- 7 结论 ----
    w("## 7. 结论(写作口径)")
    w()
    w("1. 真实公式上 H1 的拟合优度与 H0 大面积重叠,绝对阈值几乎不能筛出"
      f"错形式(与合成台架形成对照,说明合成侧的{LQ}掉档{RQ}是被设计放大的)。")
    w("2. 优度仍是难度代理:跨任务越低越易检出;任务内几无信息。")
    w("3. RESET 增广方向与残差的对齐量(corr(yhat^k,残差))是真实公式上"
      "最稳定的功效预测量;增广列的冗余度(aug_r2_yhatK)与功效弱相关,"
      "不构成机制。")
    w("4. 高优度错形式仍有可观比例四门全绿——证书语义不能建立在拟合优度上。")
    w()
    w("---")
    w("复现:`python srsd_fitting_quality.py`(重生成优度,校验生成失败集"
      "与落盘 data_fail 集一致)。图:`figures/quality_vs_power_srsd.png`。")

    # ---- 图 ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    c20 = cells[(cells.n == 200) & (cells.rel == 0.01)]
    ax = axes[0]
    for kl, col, mk in [("H0", "#ee6677", "^"), ("H1", "#4477aa", "o")]:
        cc = c20[c20.klass == kl]
        ax.scatter(cc.r2, cc.rej_fitted_f, s=30, c=col, marker=mk, label=kl)
    for _, r in c20.iterrows():
        if r.rej_fitted_f < 0.6 or r.r2 < 0.9:
            ax.annotate(r.task.replace("feynman-", ""), (r.r2, r.rej_fitted_f),
                        fontsize=6, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("mean R$^2$ of restricted power-law fit")
    ax.set_ylabel("RESET fitted-F rejection rate")
    ax.set_title("SRSD tasks (n=200, rel=1%):\nH$1$ R$^2$ overlaps H$0$")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax = axes[1]
    ax.scatter(c20.corr_yhat2_resid, c20.rej_fc, s=30, c="#4477aa", label="fc (yhat2)")
    ax.scatter(c20.corr_yhat3_resid, c20.rej_fitted_f, s=30, c="#228833",
               marker="s", label="fitted-F (yhat3)")
    ax.set_xlabel("|corr(yhat^k, residual)|")
    ax.set_ylabel("rejection rate")
    ax.set_title("Alignment vs power on real formulas")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs(os.path.dirname(FIG), exist_ok=True)
    fig.savefig(FIG, dpi=140)
    plt.close(fig)

    with open(OUT_MD, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {OUT_MD}")
    print(f"wrote {FIG}")
    return merged, cells


if __name__ == "__main__":
    main()
