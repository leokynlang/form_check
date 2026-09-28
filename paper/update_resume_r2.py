# -*- coding: utf-8 -*-
"""update_resume_r2.py — RESUME.md 追加第二轮修订的交付记录。"""
p = 'RESUME.md'
md = open(p, encoding='utf-8').read()
old = """- 交付：桌面 `Silent_Acceptance_paper_v4_revision.pdf`（31 页，覆盖更新）；
  judge 终验 4 页全 pass；工作稿 `draft_v3_complete.pdf` 同步"""
new = """- 交付：桌面 `Silent_Acceptance_paper_v4_revision.pdf`（32 页，覆盖更新）；
  judge 终验 4 页全 pass；工作稿 `draft_v3_complete.pdf` 同步
- 追加（同日）：按用户要求把修改记录收进论文——markdown 变更记录更名
  Appendix G — Change notes (v2 → v5) 并移到文末；build_pdf 停止排除该节、
  解析时 E/G 分离投递（G 在 F 代码清单之后，正文附录顺序 A–G），
  正式 PDF 现含完整两轮修订记录（第 7 条 v4 撤回与重归因、第 8 条二轮
  扩展验证）；同时修复编号条目折行成段导致的句中断开（_join_wrapped 把
  编号行视为段首吸收续行）；第 8(d) 条措辞如实化（字形审计零缺字、
  未做字体回退，审稿报告的"方框"不可复现）"""
assert old in md
open(p, 'w', encoding='utf-8').write(md.replace(old, new))
print('RESUME updated')
