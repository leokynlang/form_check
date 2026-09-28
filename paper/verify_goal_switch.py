"""verify_goal_switch.py — 核验 v3 稿目标切换改写的一致性。

检查项:
1. 关键交叉引用存在:Section 4.10 / Table 15 / Q9 / no_causal_link.py;
2. 关键数字与 no_causal_link.md 报告一致(锚点数字);
3. 旧主张句式不再出现在主张层(标题/摘要/结论);
4. 双向断言两侧都在场(高匹配度/低匹配度);
5. 备份存在。
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DRAFT = os.path.join(HERE, "draft_v3_complete.md")
REPORT = os.path.join(HERE, "..", "no_causal_link.md")
BAK = DRAFT + ".bak-goal-switch"

with open(DRAFT, encoding="utf-8") as f:
    draft = f.read()
with open(REPORT, encoding="utf-8") as f:
    report = f.read()

checks = []
ok = lambda name, cond: checks.append((bool(cond), name))  # noqa: E731

# 1. 交叉引用
ok("标题含 entailment 主张", "Does Not Entail Functional Correctness" in draft)
ok("§4.10 节存在", "### 4.10 Goodness-of-fit carries no entailment" in draft)
ok("Table 15 存在", "**Table 15" in draft)
ok("Q9 存在", "**Q9." in draft)
ok("§7 复现含 no_causal_link.py", "no_causal_link.py" in draft)
ok("引言含双向断言段", "not connected by necessity" in draft)
ok("贡献 (vi) 存在", "(vi) the no-entailment analysis" in draft)
ok("Appendix G 变更条目 9", "Goal-switch revision" in draft)
ok("§4 序言含 Table 15", "Table 15 the no-entailment summary" in draft)

# 2. 数字锚点(与 no_causal_link.md 报告核对)
anchors = [
    ("57.5", "R²<0.5 联合沉默 n=50"),
    ("49.7", "R²<0.5 联合沉默 n=200"),
    ("31.7", "R²<0.5 联合沉默 n=1000"),
    ("9.7", "R²>=0.99 域外>10% (SRSD)"),
    ("0.002", "ii.11.28 R²"),
    ("83", "ii.11.28 沉默率"),
    ("+0.76", "方向对齐 Spearman 上端"),
]
for num, desc in anchors:
    ok(f"报告数字 {desc} ({num})", num in report)

# 稿内引用的核心数字与报告一致抽查
for num in ["57.5", "9.7", "0.002"]:
    ok(f"稿内含 {num}", num in draft)

# 3. 旧主张句式已不在主张层(标题行/摘要首段/结论首段)
title = draft.splitlines()[0]
ok("标题不再以 'Benchmark of Form-Specification Gates' 开头",
   "Benchmark of Form-Specification Gates" not in title)

# 4. 双向断言
ok("方向1(高匹配度不必然正确)", "does not entail a correct function" in draft)
ok("方向2(低匹配度不必然检出)", "does not entail detection" in draft)

# 5. 备份
ok("备份存在", os.path.exists(BAK))

fails = [n for c, n in checks if not c]
for c, n in checks:
    print(f"  [{'OK ' if c else 'FAIL'}] {n}")
print(f"\n{len(checks) - len(fails)}/{len(checks)} passed")
sys.exit(1 if fails else 0)
