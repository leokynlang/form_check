"""verify_round3.py — 审稿第三轮(六条)修订的一致性核验。

检查:
1. 过强措辞已清除(主张层无 entail/independence/no information);
2. 联合事件主结果在场(Table 15 joint release + 0.67%/6.9%/7.8→0.2%);
3. SVD 替换完成(Tables 2/3/5 C/D 标注 SVD;pinv 只在 Q3 案例节);
4. 硬冲突修复(空集双口径、OR 移出 Table 12、805,825 拆分);
5. 数字与 svd_main_results.md / no_causal_link.md 锚点一致。
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DRAFT = os.path.join(HERE, "draft_v3_complete.md")
SVD_MD = os.path.join(HERE, "..", "svd_main_results.md")
NCL_MD = os.path.join(HERE, "..", "no_causal_link.md")

with open(DRAFT, encoding="utf-8") as f:
    draft = f.read()
with open(SVD_MD, encoding="utf-8") as f:
    svd_md = f.read()
with open(NCL_MD, encoding="utf-8") as f:
    ncl_md = f.read()

# 主张层 = 标题+摘要(前 120 行)+结论
head = "\n".join(draft.splitlines()[:130])
concl_start = draft.index("## 8 Conclusion")
concl_end = draft.index("## References")
concl = draft[concl_start:concl_end]
claim_layer = head + concl

checks = []


def ok(name, cond):
    checks.append((bool(cond), name))


# 1. 过强措辞清除(主张层)
ok("标题为 bounded 版", "Cannot Certify Functional Form" in draft.splitlines()[0])
ok("主张层无 'Does Not Entail'", "Does Not Entail" not in claim_layer)
ok("主张层无 'logically independent'", "logically independent" not in claim_layer)
ok("主张层无 'carry, by necessity'", "by necessity" not in claim_layer)
ok("主张层无 805,825", "805,825" not in claim_layer)
ok("摘要声明不主张独立性",
   "not statistical independence" in head or "We do not claim" in head)

# 2. 联合事件主结果
ok("Table 15 = joint release", "Joint release on the same models" in draft)
ok("0.67% 在场", "0.67%" in draft)
ok("6.9% of dangerous 在场", "6.9%" in draft)
ok("7.8% → 0.2% 在场", "7.8% → 0.2%" in draft)
ok("50x LR 在场", "50×" in draft or "50x" in draft)
ok("released dangerous 定义", "released dangerous model" in draft.lower()
   or "Released dangerous" in draft)
ok("低拟合重定位为机制", "not\nepidemiological" in draft or "not epidemiological"
   in draft.replace("\r\n", "\n").replace("\n", " "))

# 3. SVD 替换
ok("Table 2 标注 SVD", "C/D on the verified SVD path" in draft)
ok("Table 2 C 行 SVD 数字 .057", "| C (SVD) | .057" in draft)
ok("Table 3 D 行 .494/.675/.733", "| D (SVD) | **.494** | **.675** | **.733** |"
   in draft)
ok("Table 5 SVD 列", "C (SVD) | D (SVD)" in draft)
ok("§4.2 排名反转叙事", "reorders the field" in draft)
ok("pinv 仅作案例(Q3)", "numerical-implementation case" in draft)
ok("§4.1 冲突句已修", "structural miscalibration of the robustness" not in draft)

# 4. 硬冲突
ok("空集双口径", "never\r\nempty *on this pool*" in draft
   or "empty *on this pool*" in draft)
ok("Table 10 题注 V=∅ 率", "column 1 is the V = ∅ rate" in draft)
ok("OR 移出 Table 12", "excluded from this table" in draft)
ok("OR uncalibratable n50", "uncalibratable" in draft)
ok("508,825 + 297,000 拆分", "297,000" in head and "508,825" in head)

# 5. 锚点数字(与重算报告一致)
for token in ["0.67", "7.8", "0.2%", ".494", ".675", ".733", ".057"]:
    ok(f"报告含 {token}", token in svd_md or token in ncl_md)
with open(os.path.join(HERE, "..", "runs", "svd_main_results.json"),
          encoding="utf-8") as f:
    svd_json = f.read()
with open(os.path.join(HERE, "..", "runs", "cd_svd",
                       "s3_family_table_v5_output.txt"),
          encoding="utf-8") as f:
    s3_txt = f.read()
ok("报告含 6.9 (json 0.0687)", "0.0687" in svd_json)
ok("S3 表含 .860", ".860" in s3_txt or "0.860" in s3_txt)

fails = [n for c, n in checks if not c]
for c, n in checks:
    print(f"  [{'OK ' if c else 'FAIL'}] {n}")
print(f"\n{len(checks) - len(fails)}/{len(checks)} passed")
sys.exit(1 if fails else 0)
