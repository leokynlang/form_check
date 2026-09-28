# -*- coding: utf-8 -*-
"""update_resume_r2b.py — RESUME.md 追加:修改目标切换到桌面 v5 谱系。"""
p = 'RESUME.md'
md = open(p, encoding='utf-8').read()
old = """- 未做（审稿意见中明示可后置）：主文三条故事线结构性删减（电池/附录 E
  移补充材料）、Appendix D 图幅与代码字号排版、期刊版式（去封面/目录）；
  HC3 全套件重实现扫描（论文已标注 outstanding）"""
new = """- 未做（审稿意见中明示可后置）：主文三条故事线结构性删减（电池/附录 E
  移补充材料）、Appendix D 图幅与代码字号排版、期刊版式（去封面/目录）；
  HC3 全套件重实现扫描（论文已标注 outstanding）

### 重大更正与修改目标切换到桌面 v5（2026-09-08 晚，本节）
- **更正**：第二轮审稿意见审的是**桌面 v5 谱系**（`Silent_Acceptance_v5.md`，
  另一会话维护的平行稿，含 §3.7 SVD 重跑/App B 撤题/App D 案例/§3.8 双
  口径），不是我的 `paper/draft_v3_complete.md`。我此前"审稿侧幻觉"的
  F8 归因**是错的**——grep 证据为真但核对对象错版；8 条意见对 v5 有 7 条
  成立。F8 已改写为"版本错配引发的错误归因"，教训：幻觉归因前必须先做
  工件版本对齐
- v5 谱系的数据/脚本在本仓库：`runs/cd_svd/`（24 万行 SVD 重跑 +
  i.30.5 补丁 + v5_table3_calibrated.json）、`paper/revision_loc/`
- 第二轮意见在 v5 上的落实（`paper/v5_round2_analysis.py` →
  `paper/v5_round2_results.json`；编辑脚本 `apply_v5_round2.py`，14 处
  全命中）：①三处滚动检查声明改为如实口径 + D 验证扩展至 13 方程/176
  fits/三路径（SVD≡statsmodels |Δp|≤.0023）；②.067→8/120 Wilson
  [.034,.126]"名义相容"；③H1 配对 bootstrap：Δ(D−A) n50 +.080 跨零、
  n200 +.170/n1000 +.148 显著、D 胜 44–72% 任务；Δ(F−D) n200 +.003
  平局——"D 每个样本量最强"→"当前任务池最高平均功效"；④数值口径：
  NumPy 2.5.2 rcond=10⁻¹⁵ 相对值（SciPy 公式未用）、cond=10^k 损失 k 位
  语义、cond 统一为实测 10^21.8；⑤双口径 508,825/501,325 + C/D
  32-vs-31 敏感性（D 移 ≤.028）；⑥匿名链接改如实；⑦表号导航修正
  （Tables 7–9 battery 等）；⑧修订范围段移出 + 摘要压缩
- 新文件：桌面 `Silent_Acceptance_v5_response_notes.md`（修订范围承接 +
  逐条回应表）、`paper/render_desktop_v5.py`（字形修复渲染器：∅→Ø、
  ∈/∉→SimHei、逐位上标兜底）、备份 `paper/v5_backup/`
- 交付：桌面 `Silent_Acceptance_v5.pdf` 重渲染（21 页），字形审计零缺失、
  ∈ 300dpi 实拍正常、12 项内容断言全过（"rolling/rcond"两项初报 FAIL 均
  为误报：前者是否定句自身，后者是上标提取形态）
- 平行会话警告：`paper/draft_v5_complete.md`（= 我的 v3 副本）、
  `paper/render_v5_pdf.py` 属另一会话遗留，未触碰；两谱系（paper/draft*
  与桌面 v5）此后独立维护，修改前先确认目标谱系"""
assert old in md
open(p, 'w', encoding='utf-8').write(md.replace(old, new))
print('RESUME updated')
