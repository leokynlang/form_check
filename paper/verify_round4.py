"""verify_round4.py — 审稿第四轮修订的一致性核验(逐条对应意见 3.1–3.10)。

与前几轮的 verify_round*.py 同规矩:纯文本断言,全过 exit 0,任一失败
exit 1 并列出条目。**这是文本一致性检查,不重算数字**——数字重算见
`review3_analysis.py`。

检查:
  3.1  Table 15 三个分母列齐备;正文与摘要口径一致;无 "of scoped"
  3.2  §4.10(b) 列三个混淆源;§6 同款 limitation;n 梯度措辞不再无条件
  3.3  10% 阈值标为约定;合成侧 100% 标为设计产物;9.7% 标为过筛危险占比
  3.4  证书 LR 为 65×;后验先验标签正确(池混合 vs 50/50);标题含 Alone
  3.5  §6 声明未运行任何 SR 方法
  3.6  S6 效力在摘要内限定
  3.7  build_pdf 的 Appendix D 图注与图一致 + 节首代际声明;§6/§7 呼应
  3.8  §4.3 敏感度降级;Appendix C 走 SVD 口径;pinv 数字仅作路径对照
  3.9  Table 12 题注去 "matched per-gate FPR 5%";A 的预算不对称在场
  3.10 "real formulas" 在主张层清除;Table 15 题注澄清 synthetic draws
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DRAFT = os.path.join(HERE, "draft_v3_complete.md")
BUILD = os.path.join(HERE, "build_pdf.py")
ANALYSIS = os.path.join(HERE, "review3_analysis.py")

with open(DRAFT, encoding="utf-8") as f:
    draft = f.read()
with open(BUILD, encoding="utf-8") as f:
    build = f.read()

lines = draft.splitlines()
head = "\n".join(lines[:140])
concl_start = draft.index("## 8 Conclusion")
concl_end = draft.index("## References")
concl = draft[concl_start:concl_end]
claim_layer = head + concl
# 摘要 = 第一个 "## Abstract" 到 "## 1 Introduction"
abs_start = draft.index("## Abstract")
abs_end = draft.index("## 1 Introduction")
abstract = draft[abs_start:abs_end]

checks = []


def ok(name, cond):
    checks.append((bool(cond), name))


# ── 3.1 分母 ────────────────────────────────────────────────
ok("Table 15 有 A/B/C 分母定义列",
   "Screened wrong forms (A)" in draft
   and "Dangerous among them (B)" in draft
   and "Released dangerous (C)" in draft)
ok("Table 15 三个比率列齐备",
   "| C / A | C / B | C / D |" in draft)
ok("Table 15 表注定义 D 的分母",
   "D = all dangerous H1 fits" in draft
   and "12,369 / 12,348 / 12,046" in draft)
ok("Table 15 表注写明 n=50 的 616 缺口",
   "616 further dangerous fits" in draft)
ok("Table 15 表注承认旧版口径混用",
   'under the same label, "release rate of' in draft
   or "under the same label" in draft)
ok("旧 'of scoped' 措辞已清除", "of scoped" not in draft)
ok("摘要写明两个分母",
   "counting all dangerous fits" in abstract
   and "counting screened dangerous models" in abstract)
ok("§4.10(b) 报两个分母",
   "With all dangerous H1 fits as" in draft
   and "with screened dangerous models as denominator" in draft)

# ── 3.2 n 梯度 ─────────────────────────────────────────────
ok("§4.10(b) 列三个混淆源",
   "(i) gate power rises with n" in draft
   and "(ii) numerator and denominator both change with n" in draft
   and "(iii) the best" in draft)
ok("§4.10(b) 否认 'sample size closes the channel' 的强读法",
   "is not, by itself, evidence that sample size" in draft)
ok("§4.10(b) 承认 11.0→9.6→8.5 与 (iii) 一致",
   "11.0% → 9.6% → 8.5%" in draft and "with (iii) alone" in draft)
ok("§4.10(b) 指出筛选在大 n 是绑定约束",
   "the binding constraint is the fit screen" in draft)
ok("§6 有同款 limitation",
   "Threshold\nconventions and the n-gradient" in draft
   or "Threshold" in draft and "n-gradient" in draft
   and "contracting best L" in draft)
ok("旧无条件措辞已清除",
   "sample size, not any single gate, is what closes" not in draft)

# ── 3.3 阈值 ────────────────────────────────────────────────
ok("10% 阈值标为约定",
   "Both thresholds are conventions" in draft)
ok("10% 无下游误差预算",
   "carries no\ndownstream error budget" in draft
   or "downstream error budget or decision loss" in draft)
ok("headline 对标为 threshold-specific",
   "is threshold-specific" in draft)
ok("合成侧 100% 标为设计产物",
   "a design property of that bench" in draft)
ok("9.7% 标为过筛危险占比",
   "9.7% of the wrong forms" in draft
   or "9.7% of screened wrong forms are dangerous" in draft)
ok("§5 的 100%/9.7% 并置已拆开",
   "9.7% on real formulas" not in draft)

# ── 3.4 证书语义与标题 ──────────────────────────────────────
ok("证书 LR 为 65×", "65×" in draft)
ok("旧 '50× less likely' 已清除", "50× less likely" not in draft)
ok("标题含 'Alone Cannot Certify'",
   "Goodness-of-Fit Alone Cannot Certify" in lines[0])
ok("后验先验标签正确",
   "under the pool's own" in draft and "50/50 prior" in draft)
ok("旧 '2.2% at a 50/50' 标签已清除",
   "2.2% at a 50/50" not in draft)

# ── 3.5 SR 方法空缺 ─────────────────────────────────────────
ok("§6 声明未运行任何 SR 方法",
   "No symbolic-regression method is run anywhere" in draft)
ok("§6 指明这是缺口而非对照",
   "not a claim of parity with" in draft)

# ── 3.6 S6 效力 ─────────────────────────────────────────────
ok("摘要限定 S6 只给签名一致性",
   "signature-consistency claims only" in abstract)

# ── 3.7 Appendix D ─────────────────────────────────────────
ok("build_pdf: 旧 42–56% 图注清除",
   "HC3 variants sit at 42" not in build)
ok("build_pdf: 节首声明 SVD 重绘",
   "regenerated on the verified SVD path" in build)
ok("build_pdf: D.1-D.4 已换 SVD 重绘图",
   all("review3_fig_D%d.png" % i in build for i in (1, 2, 3, 4)))
ok("build_pdf: 旧四图已退役",
   "full_500_calibration.png" not in build
   and "full_500_power_heat.png" not in build
   and "full_500_n_scaling.png" not in build
   and "full_500_extrap.png" not in build)
ok("build_pdf: 封面改 budget-audited",
   "budget-audited gate comparison" in build
   and "matched-basis gate comparison" not in build)
ok("§6 有 Figure provenance 条",
   "Figure provenance (Appendix D)" in draft)
ok("§7 bit-identical 与图已同步",
   "the Appendix D panels are regenerated from" in draft)

# ── 3.8 统计单元 ────────────────────────────────────────────
ok("§4.3 敏感度降级到 30 seeds 噪声",
   "30 seeds\nper cell and cannot be separated" in draft
   or "cannot be separated from sampling noise" in draft)
ok("§4.3 family × method 区间已回填",
   "Family × method intervals are now attached" in draft)
ok("Appendix C 走 SVD 口径",
   "C .053 [.051–.056]; D .056" in draft)
ok("Appendix C 任务级 C 为 SVD 区间",
   "C [.045, .068]" in draft)
ok("Appendix C 保留 pinv 作路径对照",
   "under the retracted pinv path read C .417" in draft)
ok("§4.1 自相矛盾句已修",
   "[.045, .068], against [.197, .707] under the retracted pinv" in draft)
ok("§4.1 旧矛盾串已清除",
   "SVD-path C a 95%\ninterval of\n[.197, .707]" not in draft)

# ── 3.9 matched FPR ────────────────────────────────────────
ok("Table 12 题注去 'matched per-gate FPR 5%'",
   "at matched per-gate FPR 5%" not in draft)
ok("Table 12 后标明范围",
   '"Matched FPR" applies to three of the four columns' in draft)
ok("synth n=50 的 A 预算被点出",
   "A fires at .080" in draft)

# ── 3.10 措辞 ──────────────────────────────────────────────
ok("Table 15 题注澄清 synthetic draws",
   "synthetic draws from the" in draft
   and "registered Feynman equations, not measurements" in draft)
ok("主张层无 'real formulas'", "real formulas" not in claim_layer)

# ── 交付物 ────────────────────────────────────────────────
ok("Appendix G 有第三轮条目",
   "Third-round revision (reviewer response)" in draft)
ok("review3_analysis.py 已交付", os.path.exists(ANALYSIS))
ok("outstanding 指向 review3_analysis.py",
   draft.count("paper/review3_analysis.py") >= 4)

fails = [n for c, n in checks if not c]
for c, n in checks:
    print(f"  [{'OK ' if c else 'FAIL'}] {n}")
print(f"\n{len(checks) - len(fails)}/{len(checks)} passed")
sys.exit(1 if fails else 0)
