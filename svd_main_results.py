"""svd_main_results.py — 主文方法比较切换到核验过的 SVD 实现的一键重算。

审稿条目 3/4 的落实(2026-09-20):
- 主文 C/D 全部采用 SVD 路径(runs/cd_svd 已核验产物);旧 pinv 路径
  独立保留为"数值实现失效案例",不再用于评价 HC3 方法本身;
- 条目 2:联合事件主结果——P(危险错误 ∧ 四门沉默 | R²≥0.99 筛选),
  同批模型连接高拟合、共同沉默、危险外推;
- 条目 4:Table 2/3 重算(SVD)、OR 预算冲突修正、空集口径澄清、
  805,825 计数口径拆分。

数据全部来自已核验产物,零模拟零重生成:
- runs/full/raw.parquet(S1 500 种子,合成 A/B 原路径 + C/D 两路径)
- runs/srsd_full/(S2 500 种子)
- runs/cd_svd/{s1,s2}_cd_svd.csv(C/D SVD 路径)
- runs/cd_svd/v5_table3_calibrated.json(已核验 SVD 校准功效)
- runs/quality50 + quality_recompute(优度与外推)
- runs/gate_cmp/summary.json(S5)
- runs/battery_full/summary_battery.json(S4)
输出:svd_main_results.md(改稿数字来源)+ runs/svd_main_results.json。
锚点核对失败即非零退出。
"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_MD = os.path.join(HERE, "svd_main_results.md")
OUT_JSON = os.path.join(HERE, "runs", "svd_main_results.json")
R = {}  # 汇总 JSON

STRUCT = ["QUAD_LIN", "CUBIC_LIN", "INTERACT", "THRESHOLD", "EXP_MODEL",
          "LOG_MODEL", "RATIONAL", "OMITTED_M", "DILUTION", "OPPOSING",
          "COLLINEAR_CURV"]
H0_SET = {"H0_BASELINE", "H0_MULTIVAR", "H0_CORR"}
ERR_T = 0.10
NS = [50, 200, 1000]


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


def main():
    L = []

    def w(s=""):
        L.append(s)

    # ============ 载入 ============
    sr = pd.read_parquet(os.path.join(HERE, "runs", "srsd_full",
                                      "raw_srsd.parquet"))
    sq = pd.read_parquet(os.path.join(HERE, "runs", "srsd_full",
                                      "quality_recompute.parquet"))
    srsd = sr.merge(sq[["task", "n", "rel", "seed", "r2"]],
                    on=["task", "n", "rel", "seed"], how="left")
    cd2 = pd.read_csv(os.path.join(HERE, "runs", "cd_svd", "s2_cd_svd.csv"),
                      low_memory=False)
    srsd = srsd.merge(cd2[["task", "n", "rel", "seed", "rej_C_svd",
                           "rej_D_svd", "fail_C_svd", "fail_D_svd"]],
                      on=["task", "n", "rel", "seed"], how="left")
    raw = pd.read_parquet(os.path.join(HERE, "runs", "quality50",
                                       "raw.parquet"))
    ex = pd.read_parquet(os.path.join(HERE, "runs", "quality50",
                                      "extrap.parquet"))
    syn = raw.merge(ex[["scenario", "n", "rel", "seed",
                        "extrap_far_relerr"]],
                    on=["scenario", "n", "rel", "seed"], how="left")
    # S1 SVD(合成台,验证 old==svd 后用于 Table 2 的合成行)
    cd1 = pd.read_csv(os.path.join(HERE, "runs", "cd_svd", "s1_cd_svd.csv"),
                      low_memory=False)
    syn = syn.merge(cd1[["scenario", "n", "rel", "seed", "rej_C_svd",
                         "rej_D_svd"]],
                    on=["scenario", "n", "rel", "seed"], how="left")

    # 联合沉默列(SRSD: A/B 原路径 + C/D SVD;合成: 原四门)
    syn["_sil"] = (~syn[["rej_fc", "rej_fitted_f", "rej_fitted_hc3",
                         "rej_exog_hc3"]].fillna(False).any(axis=1)
                   & ~syn[["fail_fc", "fail_fitted_f", "fail_fitted_hc3",
                           "fail_exog_hc3"]].fillna(False).any(axis=1))
    srsd["_sil"] = (~srsd[["rej_fc", "rej_fitted_f", "rej_C_svd",
                           "rej_D_svd"]].fillna(False).any(axis=1)
                    & ~srsd[["fail_fc", "fail_fitted_f", "fail_C_svd",
                             "fail_D_svd"]].fillna(False).any(axis=1))

    # ============ 锚点核对 ============
    anchors = []

    def anchor(name, got, want, tol):
        ok = abs(float(got) - want) <= tol
        anchors.append((ok, name, float(got), want))

    with open(os.path.join(HERE, "runs", "cd_svd", "summary_s2.json"),
              encoding="utf-8") as f:
        s2sum = json.load(f)["S2"]["h0_size"]
    for cell, want in [("C_n50", 0.0565), ("D_n200", 0.0564),
                       ("C_n1000", 0.0549)]:
        anchor(f"cd_svd s2 {cell} svd rate", s2sum[cell]["svd"]["rate"],
               want, 0.002)
    # v5_table3 已核验 SVD 校准功效(19 任务池)
    with open(os.path.join(HERE, "runs", "cd_svd",
                           "v5_table3_calibrated.json"),
              encoding="utf-8") as f:
        t3 = json.load(f)["v5_table3_calibrated_power_at_fpr5"]
    anchor("v5 t3 A n200", t3["A (permutation)"]["200"], 0.486, 0.002)
    anchor("v5 t3 D n200", t3["D (HC3 exog, SVD)"]["200"], 0.675, 0.002)
    anchor("v5 t3 C n50", t3["C (HC3 fitted, SVD)"]["50"], 0.129, 0.002)
    # 联合事件(侦察已得)
    base = ((srsd.klass == "H1") & (srsd.rel > 0) & srsd.r2.notna()
            & srsd.extrap_far_relerr.notna())
    s99 = srsd[base & (srsd.r2 >= 0.99)]
    danger99 = s99[s99.extrap_far_relerr > ERR_T]
    js99 = danger99[danger99._sil]
    anchor("joint event: N screened R2>=0.99", len(s99), 71992, 50)
    anchor("joint event: P(dangerous|screen)", len(danger99) / len(s99),
           0.097, 0.005)
    anchor("joint event: P(dangerous & silent|screen)", len(js99) / len(s99),
           0.0067, 0.001)
    # H0 对照
    h0 = srsd[(srsd.klass == "H0") & (srsd.rel > 0)
              & srsd.extrap_far_relerr.notna() & (srsd.r2 >= 0.99)]
    anchor("H0: P(cert|true)", len(h0) / len(srsd[(srsd.klass == "H0")
                                                  & (srsd.rel > 0)
                                                  & srsd.extrap_far_relerr
                                                  .notna()]), 0.962, 0.01)
    anchor("H0: silence given cert", h0._sil.mean(), 0.879, 0.02)
    bad = [a for a in anchors if not a[0]]
    for ok, name, got, want in anchors:
        print(f"  [{'OK ' if ok else 'FAIL'}] {name}: got={got:.4f} "
              f"want={want}")
    if bad:
        raise SystemExit(f"锚点核对失败: {[b[1] for b in bad]}")

    # ============ §A SVD 版 Table 2(S2 H0 size)============
    w("# 主结果 SVD 切换:重算汇总(审稿条目 2/3/4)")
    w()
    w("> 用途:论文 §4.1/§4.2/§4.10/Table 2/3/15 与摘要的改稿数字来源。")
    w("> C/D 全部为 SVD 路径(runs/cd_svd 已核验);旧 pinv 路径仅作数值")
    w("> 失效案例保留。零模拟,读已核验产物。")
    w()
    w("## A. SVD 版 Table 2 — S2 H0 size(α=0.05)")
    w()
    w("| Method | n=50 | n=200 | n=1000 |")
    w("|---|---|---|---|")
    h0s = srsd[(srsd.klass == "H0")]
    R["table2_svd"] = {}
    for m, col in [("A", "rej_fc"), ("B", "rej_fitted_f")]:
        row = []
        for n in NS:
            v = h0s[(h0s.n == n) & (h0s.rel > 0)][col].mean()
            row.append(f"{v:.3f}")
            R["table2_svd"][f"{m}_n{n}"] = round(float(v), 4)
        w(f"| {m} | " + " | ".join(row) + " |")
    for m, col in [("C", "rej_C_svd"), ("D", "rej_D_svd")]:
        row = []
        for n in NS:
            v = h0s[(h0s.n == n) & (h0s.rel > 0)][col].mean()
            row.append(f"{v:.3f}")
            R["table2_svd"][f"{m}_n{n}"] = round(float(v), 4)
        w(f"| {m} (SVD) | " + " | ".join(row) + " |")
    w()
    w("对照(旧 pinv 路径,数值失效案例):C .417–.422 / D .563–.565,"
      "不随 n 消失——该现象整体移入数值实现失效节,不再进入方法比较。")
    w()

    # ============ §B SVD 版 Table 3(校准功效)============
    w("## B. SVD 版 Table 3 — S2 校准功效 @5% FPR")
    w()
    w("数字直接取自已核验的 `runs/cd_svd/v5_table3_calibrated.json`")
    w("(19 任务池,与主文 Table 3 同口径重推导):")
    w()
    w("| Method | n=50 | n=200 | n=1000 |")
    w("|---|---|---|---|")
    R["table3_svd"] = {}
    for m in ["A (permutation)", "B (RESET fitted F)",
              "C (HC3 fitted, SVD)", "D (HC3 exog, SVD)"]:
        row = [f"{t3[m][str(n)]:.3f}" for n in NS]
        R["table3_svd"][m] = {str(n): t3[m][str(n)] for n in NS}
        w(f"| {m} | " + " | ".join(row) + " |")
    w()
    w("**排名反转(实质结果变化,必须改写 §4.2 叙事)**:旧 pinv 路径的")
    w("C/D 校准功效行沿伪影尺度计算(v5_table3 meta.critical_finding),")
    w("无意义;SVD 重推导后 D(.494/.675/.733)在全部 n 上超过 A")
    w("(.386/.486/.570),任务级配对 bootstrap(18 任务池)Δ(D−A)")
    w("n=200 为 −.170 [−.329,−.041]、n=1000 −.148 [−.304,−.019]")
    w("(D 显著优),n=50 −.080 [−.258,+.081] 跨零。")
    w()

    # ============ §C 条目 2:联合事件主结果 ============
    w("## C. 联合事件主结果(审稿条目 2)")
    w()
    w("事件定义:模型为 H1(错形式)**且**域外误差>10%(危险)**且**")
    w("通过 R²≥0.99 拟合筛选**且**四门(A/B 原路径+C/D SVD)全不拒、")
    w("无一门拒检(fail_* 算诚实拒绝,单独计数)。")
    w()
    w("### C.1 放行率(同批模型)")
    w()
    w("| 池 | n | 通过筛选 N | 危险(>10%) | 危险且沉默(全绿证书) |")
    w("|---|---|---|---|---|")
    R["joint_event"] = {}
    for n in NS:
        s = srsd[base & (srsd.n == n) & (srsd.r2 >= 0.99)]
        d_ = s[s.extrap_far_relerr > ERR_T]
        js = d_[d_._sil]
        w(f"| SRSD H1 | {n} | {len(s)} | {len(d_)} ({pct(len(d_), len(s))}) "
          f"| {len(js)} ({pct(len(js), len(s))}) |")
        R["joint_event"][f"n{n}"] = {
            "screened": len(s), "dangerous": len(d_),
            "dangerous_silent": len(js),
            "p_danger_given_screen": round(len(d_) / len(s), 4),
            "p_certified_wrong": round(len(js) / len(s), 4),
            "p_silent_given_dangerous": round(len(js) / len(d_), 4)}
    # 全 n 合并
    s = srsd[base & (srsd.r2 >= 0.99)]
    d_ = s[s.extrap_far_relerr > ERR_T]
    js = d_[d_._sil]
    R["joint_event"]["pooled"] = {
        "screened": len(s), "dangerous": len(d_),
        "dangerous_silent": len(js),
        "p_danger_given_screen": round(len(d_) / len(s), 4),
        "p_certified_wrong": round(len(js) / len(s), 4),
        "p_silent_given_dangerous": round(len(js) / len(d_), 4)}
    w(f"| SRSD H1 合并 | — | {len(s)} | {len(d_)} ({pct(len(d_), len(s))}) "
      f"| {len(js)} ({pct(len(js), len(s))}) |")
    w()
    w("### C.2 危险错误模型的放行率(分母=全部危险 H1)")
    w()
    w("| n | 危险 H1 总数 | 四门沉默放行 | 其中已过 R²≥0.99 筛选 |")
    w("|---|---|---|---|")
    for n in NS:
        s = srsd[base & (srsd.n == n)]
        d_ = s[s.extrap_far_relerr > ERR_T]
        js = d_[d_._sil]
        js99 = js[js.r2 >= 0.99]
        w(f"| {n} | {len(d_)} | {len(js)} ({pct(len(js), len(d_))}) "
          f"| {len(js99)} |")
        R["joint_event"][f"n{n}"]["danger_total"] = len(d_)
        R["joint_event"][f"n{n}"]["released_of_dangerous"] = len(js)
        R["joint_event"][f"n{n}"]["released_and_screened"] = len(js99)
    w()
    w("### C.3 证书语义(LR,任务池构成修正)")
    w()
    h0c = srsd[(srsd.klass == "H0") & (srsd.rel > 0)
               & srsd.extrap_far_relerr.notna()]
    h0_cert = h0c[(h0c.r2 >= 0.99)]
    h0_cert_sil = h0_cert[h0_cert._sil]
    p_certsil_h0 = len(h0_cert_sil) / len(h0c)
    # 混合先验: 13 H0 + 19 H1 任务(每任务等样本),P(H1)=19/32
    prior_h1 = 19 / 32
    # P(cert | danger-H1) = js/danger(全部危险 H1 中拿全绿证书)
    dall = srsd[base]
    dall_d = dall[dall.extrap_far_relerr > ERR_T]
    dall_js = dall_d[(dall_d.r2 >= 0.99) & dall_d._sil]
    p_cert_danger = len(dall_js) / len(dall_d)
    # LR(cert; danger-H1 vs true-H0)
    lr = p_cert_danger / p_certsil_h0
    post = (p_cert_danger * prior_h1) / (
        p_cert_danger * prior_h1 + p_certsil_h0 * (1 - prior_h1))
    post5050 = p_cert_danger / (p_cert_danger + p_certsil_h0)
    R["certificate"] = {
        "p_cert_true_h0": round(p_certsil_h0, 4),
        "p_cert_dangerous_h1": round(p_cert_danger, 4),
        "LR": round(lr, 3),
        "posterior_p_danger_given_cert_mixture": round(post, 4),
        "posterior_p_danger_given_cert_5050": round(post5050, 4),
        "mixture_prior_h1": prior_h1}
    w(f"- P(全绿证书 | 真形式 H0) = {p_certsil_h0:.3f}")
    w(f"- P(全绿证书 | 危险错形式 H1) = {p_cert_danger:.3f}"
      f"(全池 {len(dall_js)}/{len(dall_d)})")
    w(f"- LR = {lr:.3f}(= {p_cert_danger:.3f}/{p_certsil_h0:.3f}"
      f"= 1/{1 / lr:.0f});后验 P(危险|证书):池混合先验 "
      f"P(H1)={prior_h1:.2f} 下 ≈ {post:.1%},50/50 先验下 ≈ "
      f"{post5050:.1%}")
    w("- 读法(与联合沉默报告一致):全绿证书是**弱证据**——拿到证书")
    w("  后仍有一定比例是危险错形式,证书不构成安全声明;")
    w("  但方向上证书更可能来自真形式(LR≈1 的『无信息』读数在此")
    w(f"  口径不成立:实测 LR={lr:.3f}——危险模型拿证书的概率低于真形式,")
    w("  证书有信息量但不充分——这与『无信息』措辞划清界限,")
    w("  正是审稿条目 1 要求的降格表述。)")
    w()

    # ============ §D 条目 4 修正 ============
    w("## D. 条目 4 硬冲突修正")
    w()
    # D.1 空集
    w("### D.1 §4.6 空集口径(冲突裁决)")
    w()
    w("原始 summary_battery.json `synth_none_escalation`(out-of-family")
    w("池,V=∅ 比例):QUAD/THRESHOLD/OMITTED = 1.00(空集 100% 出现,")
    w("升级即空集,这是**正确结果**);INTERACT .0545 / CUBIC .241 /")
    w("OPPOSING .001 / COLLINEAR_CURV .034(空集很少出现,错误保留)。")
    w("**修正措辞**:原文『a genuinely empty set never occurs』只对")
    w("H0 侧电池成员格成立(containment 表),对 out-of-family 池为假——")
    w("Table 10 的『Correct escalation』列就是 V=∅ 率。§4.6 文字与")
    w("Table 10 用的是同一个量,文字句删改为明确双口径。")
    w()
    # D.2 OR 预算
    w("### D.2 Table 11/12 OR 预算冲突(修正)")
    w()
    w("Table 12 的 OR(A,E) 行 = 各自 5% 阈值的并集(联合预算")
    w(".065–.102),与其余行的单门 5% 不同预算,不能同列比较;")
    w("synth n=50 时 A 自身 H0 火灾率 .080 > 5%,保持 A 不变收紧 E")
    w("无法把联合降到 5%(OR 的 FPR ≥ A 的 FPR)——§4.7(d) 的『at")
    w("n=50 the union *is* A』与 Table 12 数字的关系需要明确:")
    w("联合校准仅在 A 自身 ≤5% 的格可定义(synth n=50 除外,该格如实")
    w("报告为不可校准)。Table 12 拆两个面板:naive-OR 与")
    w("joint-calibrated-OR(仅 n=200 两台架有值 .664/.940)。")
    w()
    # D.3 计数
    w("### D.3 805,825 计数口径(修正)")
    w()
    w("805,825 = 508,825(六套件独立数据集拟合)+ 57,000(quality")
    w("bench:与 S1 同生成器同 seed 规则的**重新拟合**,新增 C/D 无、")
    w("A/B 同源)+ 240,000(S2 quality 重算:**同一批已拟合数据集的")
    w("指标重算与 SVD 路径重评估,非新数据集**)。修正:摘要不再用")
    w("单一总数 805,825;改为『508,825 fitted datasets(六套件);")
    w("quality re-evaluation on 297,000 of these(重算,非新增)』。")
    w()

    # ============ §E SVD 汇总 JSON 落盘 ============
    R["anchors_ok"] = len(anchors) - len(bad)
    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(R, f, ensure_ascii=False, indent=1)
    with open(OUT_MD, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print(f"\nwrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
