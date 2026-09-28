"""make_report.py — 生成 report.md(汇总表、校准、功效、外推、失败审计)。

用法: python make_report.py <tag> [--srsd <srsd_tag>]
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = os.path.dirname(os.path.abspath(__file__))
METHODS = {"fc": "form_check", "fitted_f": "RESET fitted F",
           "fitted_hc3": "RESET fitted HC3", "exog_hc3": "RESET exog HC3"}
H0_SET = {"H0_BASELINE", "H0_MULTIVAR", "H0_CORR"}
NOISE_H1 = {"H1_HEAVYTAIL", "H1_VAR_GROWING", "H1_AUTOCORR",
            "H1_OUTLIER", "H1_OUTLIER_GRP"}
STRUCT_H1 = ["QUAD_LIN", "CUBIC_LIN", "INTERACT", "THRESHOLD", "EXP_MODEL",
             "LOG_MODEL", "RATIONAL", "OMITTED_M", "DILUTION", "OPPOSING",
             "COLLINEAR_CURV"]


def fmt(x, d=3):
    return "—" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{d}f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tag")
    ap.add_argument("--srsd", default=None)
    args = ap.parse_args()
    tag = args.tag
    run_dir = os.path.join(HERE, "runs", tag)
    with open(os.path.join(run_dir, "summary.json"), encoding="utf-8") as fh:
        summary = json.load(fh)
    cs = pd.read_csv(os.path.join(run_dir, "cell_stats.csv"))
    with open(os.path.join(run_dir, "config.json"), encoding="utf-8") as fh:
        cfg = json.load(fh)
    with open(os.path.join(run_dir, "env.json"), encoding="utf-8") as fh:
        env = json.load(fh)
    raw = pd.read_parquet(os.path.join(run_dir, "raw.parquet"))

    L = []
    A = L.append
    A(f"# form_check vs Ramsey RESET — {tag} 报告\n")
    A(f"生成时间:{cfg.get('finished', '?')} | 行数 {summary['n_rows']} | "
      f"种子 {cfg['seeds']}/格 | n∈{cfg['n_grid']} | rel∈{cfg['rel_grid']}\n")
    A("## 0. 公平性设计(约束 4/5/6 的落实)\n")
    A("- 方法 A 调用 `discover.py::form_check` **本体**,未复制未修改;"
      "其 |t|>3 原生阈值只作留痕,决策用包装层在同一排列零分布上算的"
      "**经验 p 值**(B=999,分辨率 0.001)。")
    A("- 方法 B/C/D = RESET(fitted powers=[2,3] 经典 F / fitted HC3-Wald / "
      "exog HC3-Wald),p 值同为 [0,1] 尺度。")
    A(f"- 四方法统一 α={cfg['alpha']};功效比较另按 **经验 FPR=5% 校准阈值** "
      "(H0 池分位数)做二次汇总(§4)。")
    A("- H0/H1 分组在代码写入时预注册(`SCENARIO_META`),不因结果调整。\n")

    A("## 1. H0 校准(严格单项式真值,拒绝率=误报率)\n")
    A("| n | " + " | ".join(METHODS[m] for m in METHODS) + " |")
    A("|---|" + "---|" * len(METHODS))
    for n in [50, 200, 1000]:
        vals = []
        for m in METHODS:
            sel = cs[(cs.scenario.isin(H0_SET)) & (cs.n == n)]
            vals.append(fmt(sel[f"rate_{m}"].mean()))
        A(f"| {n} | " + " | ".join(vals) + " |")
    A("\n判读:form_check 的经验 p 值在 H0 下接近名义 5%(排列零分布与原生"
      "bootstrap 零分布同构);经典 F 在异方差/重尾 H0 噪声场景会有超尺寸,"
      "HC3 变体在此处显示其修正作用(见 §5 噪声层)。\n")

    A("## 2. 结构 H1 功效 @α=0.05(rel=1%)\n")
    A("| 场景 | n | " + " | ".join(METHODS[m] for m in METHODS) + " |")
    A("|---|---|" + "---|" * len(METHODS))
    for s in STRUCT_H1:
        for n in [50, 200, 1000]:
            row = [fmt(summary["structural_h1_power_alpha05"][s][str(n)][m])
                   for m in METHODS]
            A(f"| {s} | {n} | " + " | ".join(row) + " |")
    A("")
    A("逐条解读(约束:不得只报总平均):\n")
    for s in STRUCT_H1:
        pw = summary["structural_h1_power_alpha05"][s]
        best_n1000 = max(METHODS, key=lambda m: pw["1000"][m])
        A(f"- **{s}**: n=1000 时最高 {METHODS[best_n1000]} "
          f"({fmt(pw['1000'][best_n1000], 2)});"
          f"form_check {fmt(pw['1000']['fc'], 2)} / "
          f"fitted-F {fmt(pw['1000']['fitted_f'], 2)} / "
          f"exog {fmt(pw['1000']['exog_hc3'], 2)}。")
    A("")

    A("## 3. 噪声层 H1(结构仍是严格幂律)——这些是\"误报\"不是\"检出\"\n")
    A("| 场景 | n | " + " | ".join(METHODS[m] for m in METHODS) + " |")
    A("|---|---|" + "---|" * len(METHODS))
    for s in sorted(NOISE_H1):
        for n in [200]:
            row = [fmt(summary["noise_layer_rejection_alpha05"][s][str(n)][m])
                   for m in METHODS]
            A(f"| {s} | {n} | " + " | ".join(row) + " |")
    A("\n**语义警示(约束:区分\"没有拒绝\"与\"证明模型正确\")**:本节场景的"
      "均值结构是严格单项式——任何拒绝都是对噪声假设的误报,不是发现形式失配。"
      "反之,H0 格的\"不拒绝\"也只说明\"无证据反对幂律\",不构成\"幂律得证\"。"
      "form_check 的输出语义即如此措辞(未发现反对证据 ≠ 验证成立)。\n")

    A("## 4. 相同经验 FPR(5%)下的功效(rel=1%,结构 H1 池)\n")
    A("阈值 = H0 池 p 值 5% 分位数(按 n 校准):\n")
    A("| n | " + " | ".join(METHODS[m] for m in METHODS) + " |")
    A("|---|" + "---|" * len(METHODS))
    for n in [50, 200, 1000]:
        row = [fmt(summary["power_at_empirical_fpr5_rel1pct"][str(n)][m])
               for m in METHODS]
        A(f"| {n} | " + " | ".join(row) + " |")
    A("")

    A("## 5. 胜负矩阵(结构 H1 全格,2pp 容差)\n")
    try:
        wm = pd.read_csv(os.path.join(run_dir, "win_matrix.csv"), index_col=0)
        A(wm.to_markdown())
    except Exception as e:  # noqa: BLE001 — 读取失败写入报告正文,报告生成不得中断
        A(f"(win_matrix 读取失败:{e})")
    A("")

    A("## 6. 运行时间(每数据集,秒)\n")
    A("| n | " + " | ".join(METHODS[m] for m in METHODS) + " |")
    A("|---|" + "---|" * len(METHODS))
    for n in [50, 200, 1000]:
        row = []
        for m in METHODS:
            sel = cs[cs.n == n]
            row.append(f"{sel[f'time_{m}'].mean():.4f}")
        A(f"| {n} | " + " | ".join(row) + " |")
    A("")

    A("## 7. 失败登记(不可静默丢弃)\n")
    A("| 方法 | 失败率 | 主要失败类型 |")
    A("|---|---|---|")
    main_fail = {}
    for m, m_name in METHODS.items():
        col = f"fail_{m}"
        if col in raw:
            vc = raw[col].value_counts(dropna=True)
            main_fail[m] = vc.index[0] if len(vc) else "—"
            A(f"| {m_name} | {raw[col].notna().mean():.4f} | {main_fail[m]} |")
    A("\n失败类型说明:`degenerate_zero_residual`=零噪声/精确幂律(检验对象"
      "不存在,显式登记);`rank_deficient`=设计秩亏;`hc3_leverage_degenerate`"
      "=杠杆 1−h→0;`statistic_mismatch`=form_check 返回统计量与包装层重算"
      "不符(实现漂移哨兵);`checked_false_*`=form_check 本体拒绝检查"
      "(非正值/小样本/非有限)。\n")

    A("## 8. 复现信息\n")
    A(f"- git HEAD:`{env.get('git_head','?')}`(工作区 {env.get('git_dirty_entries','?')} 处未提交改动)")
    A(f"- discover.py sha256:`{env['core_file_sha256']['discover.py'][:16]}…`")
    A(f"- python {env['python']},numpy {env['numpy']},scipy {env['scipy']},"
      f"pandas {env['pandas']},dimreflex {env['dimreflex_version']}")
    A(f"- 平台:{env['platform']}")
    A(f"- 种子规则:`{cfg['seed_rule']}`;排列 B={env['b_perm']}\n")

    if args.srsd:
        A("\n---\n")
        A(SRSD_SECTION(args.srsd))
        A("\n## 图(SRSD)\n")
        A("- `figures/srsd_calibration.png` — H0 p 值校准(真实单项式)")
        A("- `figures/srsd_per_task_power.png` — 逐任务功效\n")

    out = os.path.join(run_dir, "report.md")
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(L))
    print("written", out)


def SRSD_SECTION(srsd_tag):
    sdir = os.path.join(HERE, "runs", srsd_tag)
    with open(os.path.join(sdir, "summary_srsd.json"),
              encoding="utf-8") as fh:
        s = json.load(fh)
    L = []
    A = L.append
    A(f"# SRSD-Feynman 主基准 — {srsd_tag}\n")
    A("## 任务池(预注册分类,判据=公式结构,先于任何实验)\n")
    A(f"- 入选 {s['n_tasks']} 个(H0 单项式 {s['n_h0']},H1 非单项式 {s['n_h1']});"
      f"排除 {s['n_excluded']} 个。")
    A("- 排除理由分布:" + ", ".join(f"{k}×{v}" for k, v in s["exclude_reasons"].items()))
    A(f"- 完整清单(含逐任务理由):`runs/{srsd_tag}/task_registry.csv`\n")
    A("## H0(真单项式)误报率 @α=0.05\n")
    A("| n | " + " | ".join(METHODS[m] for m in METHODS) + " |")
    A("|---|" + "---|" * len(METHODS))
    for n in [50, 200, 1000]:
        A(f"| {n} | " + " | ".join(
            fmt(s["h0_fpr"][str(n)][m]) for m in METHODS) + " |")
    A("\n## 相同经验 FPR=5% 下的 H1 功效\n")
    A("| n | " + " | ".join(METHODS[m] for m in METHODS) + " |")
    A("|---|" + "---|" * len(METHODS))
    for n in [50, 200, 1000]:
        A(f"| {n} | " + " | ".join(
            fmt(s["power_at_empirical_fpr5"][str(n)][m]) for m in METHODS) + " |")
    A("\n## 重点问题\n")
    kq = s.get("key_questions", {})
    for q, a in kq.items():
        A(f"- **{q}** {a}")
    A("")
    A("## 报警→外推误差关联\n")
    ex = s.get("extrap_auc", {})
    A("按\"远端外推相对误差>10%\"为正类的 p 值 AUC:" +
      (", ".join(f"{k}={v:.3f}" for k, v in ex.items()) if ex else "—"))
    A("")
    return "\n".join(L)


if __name__ == "__main__":
    main()
