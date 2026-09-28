# 第八部分 论文与文档构建层(`build_hx_paper.py`、`artifacts/`、`experiments/form_check_vs_reset/paper/`)

> 本部分覆盖 4 个构建脚本(约 1,300 行)与它们的产物。这一层不产生新的统计证据,而是把前七部分的证据**排版成可提交的手稿**,并且同样体现项目的纪律:数据直接从实验产物 CSV 读取(不手抄数字)、图表程序化绘制、诚实的证据分层声明写进正文。

---

## 8.1 `build_hx_paper.py` — HX 论文 Word 构建器(747 行)

**产物**:`artifacts/HX_Hybrid_Directional_Specification_Test.docx`(匿名稿 "HX: A Hybrid Directional Specification Test for Log-Linear Physical Model Discovery")。

### 8.1.1 排版基础设施(build_hx_paper.py:39-293)

- **OOXML 级精确控制**:`set_cell_shading/set_cell_margins/set_table_borders`(直接操作 `w:shd/w:tcMar/w:tblBorders`)、表头跨页重复(`w:tblHeader`)、行防拆分(`w:cantSplit`)。
- **真 OMML 数学公式**:`add_equation` 构造 `m:oMathPara/m:oMath/m:r/m:t` 元素(build_hx_paper.py:147-160)——Word 原生公式对象而非图片或斜体文本,8 个编号公式(H0 定义/退化判据/方向池/正交化/maxT/联合 F/置换 p/Bonferroni 合并)。
- **`add_field`**:`PAGE` 域代码动态页码;Title 样式的默认蓝色底边框被显式移除(注释:"学术标题页只用排版")。
- **样式系统**:`setup_styles` 全套 Aptos/等线字体、1.5 倍行距、两端对齐、widow control;`setup_page` 信纸尺寸+0.72-0.78in 边距;页脚页码、页眉留空。
- **程序化插图**:Pillow 手绘两张 1500px 图——`draw_performance_chart`(FPR/TPR 双面板条形图,α=0.05 红线)与 `draw_workflow`(HX 六步流水线框图:Fit base→Generate→Residualize→Select→maxT/Joint F 分叉→Bonferroni Combine);`font()` 在 aptos/arial 字体文件间回退。

### 8.1.2 数据直读与摘要数字(295-320)

- `load_random_data()` 直接读 `experiments/fc_multidirection/runs/random_1000.csv`——**论文数字不是手抄的**,1000 组泛化结果改动后重跑构建器即自动更新。
- `metric_summary` 按 klass∈{H0,H1} 且失败字段为空过滤,逐方法计 FPR/TPR;`wilson`(z=1.959963984540054 精确值)算 95% CI。
- **摘要关键数字**(build_hx_paper.py:438):1000 组(393 H0/607 H1)上 HX FPR=0.0509(Wilson 0.0332-0.0773)、power=0.9522(0.9322-0.9665),对照单方向 form_check power 0.4942、结构-only MD4 0.8913、RESET 0.7512;并诚实声明"随机套件是 stress test 而非独立外部验证;冻结 fitted-value 方向使残差置换校准只是近似"。

### 8.1.3 正文结构与诚实条款

- **结构**:摘要→Plain Language Summary→目录→11 章(引言/相关工作/问题形式化/HX 核心算法/统计解释与守卫/评估设计/结果/机理/局限/可复现协议/结论)+ 附录 A 算法规格(10 步)、B 失败分类学(六类命名失败码)、C 最小报告清单(8 条)、13 条参考文献(Ramsey 1969/MacKinnon-White/Freedman-Lane/Winkler/wild bootstrap/Westfall-Young/Buckingham-Pi/SRBench/AI Feynman/PySR/SRSD)。
- **方法论上值得记录的三处**:
  1. **§4.6 诚实声明**:开发实现用"冻结方向 + 受限残差置换",fitted-value 方向是响应自适应的,"其有限样本精确性不保证,p 值应读作校准近似";确认性实现应改 refit mode(每次置换内重建)或 split mode;异方差换 wild bootstrap。
  2. **§5.1 错误控制边界**:"maxT 控制的是方向族内搜索;Bonferroni 控制的是两分支并集;两者都不修复错定的重采样律"。
  3. **§10 证据分层表**(build_hx_paper.py:652-657):四档证据各配"允许的主张"——Development benchmark / Post-development randomized stress test / **Locked synthetic confirmation** / External scientific datasets——把"这个数字能说什么"制度化。
- **失败可见性进正文**(§7.5):F10 上 HX 秩亏失败而 RESET 告警("为什么失败状态必须保持可见");L1(长利)n=16、残差 df=3、条件数 4.9e9,"HX 拒绝 RESET 不拒绝,但数值域使其不可归因而非干净的方法胜利"。

---

## 8.2 `artifacts/build_unified_paper.py` — 统一手稿构建器(产 `Falsification_Side_of_Physical_Model_Discovery.docx`)

- **【计算】** 把两篇工作合成一篇 Word:"Silent Acceptance(证伪侧基准)"+"HX(混合方向检验)"统一为 *Falsification Side of Physical Model Discovery*。中文排版(宋体/黑体、eastAsia rFonts 显式设置、首行缩进 24pt、1.5 倍行距、A4)。
- **【原因】** 两篇手稿共享同一个证据链(Silent Acceptance 的"指标集体静默"现象学 → HX 的"修复方案"),统一稿面向同一读者;单独 docx 保留在 artifacts/ 供分别投稿。
- 注意:`OUT = r"D:\dimreflex\artifacts\..."` 硬编码 D 盘路径(与 build_hx_paper.py 的 `ROOT = Path(r"D:\dimreflex")` 一致)——构建器曾在另一工作副本(D:\dimreflex)上运行,产物复制回本仓库的 artifacts/。

## 8.3 `experiments/form_check_vs_reset/paper/` — Silent Acceptance 手稿

### 8.3.1 论文草稿(draft_v0_annotated → v1_complete → v2_complete)

- **`draft_v2_complete.md`**(当前版)标题:*Silent Acceptance: A Falsification-Side Benchmark of Form-Specification Tests for Physical Model Discovery*。开头即声明 "Complete draft v2 — working manuscript, not submitted",**变更日志内嵌**(v2 相对 v1:电池章节合并、外推章节两处错误修正、表编号修复)。
- **四个规模数字**:三套件 1,157 网格格 / 394,260 拟合数据集 / 472,500 留出评估(240,000 同分布测试 + 232,500 外探针);第四套件(电池)188 格 / 94,000 数据集。
- **四个主要发现**(摘要):①聚合排名跨任务族不迁移(fitted-F 在合成二次上 power≥0.99,对有理真值近盲 0.13-0.19);②HC3-Wald 在真单项式物理定律上误报 42-56%、**不随 n 变平**(n=50→1000),匹配经验 FPR 后优势消失;③无测试处处最优——排列曲率检验两池校准最好(5%±1%)但漏 exponential 族且有可构造盲区;各法失败画像互补。报警提高严重外推误差的概率(精确率 0.39-0.56 对 0.20 基率)但不预测其幅度;④电池:同一统计量可在族间路由——**但必须逐族重建零分布**(排列零分布在 log-log 坐标精确,在 rational 的 (1/x,1/y) 坐标产生 30-34% 假升级,wild bootstrap 修到 4.5-9%;族内真值包含 93-95%,四族最外结构 1.00 正确升级)。
- **协议三承诺**(§1):统一 α;共享排列零分布;只在匹配经验 FPR(真单项式任务上估)下比功效——"我们第一个 pilot 恰恰犯了直接比较的错误,产生了共同基线下消失的倒置排名"(pilot 的失败被写进论文动机)。
- **论文标题现象**(§1 开头):幂律拟合 PK 数据 R²=0.99、过 k-fold、赢 IC 比较、标准误极小——域外一个十进位 −20% 误差;"这不是假设:失败作为受控场景在实验室笔记中复现,并作为场景 S1 出现在本基准中(四个检验 power≥0.99 全检出)"。

### 8.3.2 构建脚本

- **`build_pdf.py`**(ReportLab 路线):读 draft_v2_complete.md(排除过程性内容)+ mindmap.png + figures/full_500_*.png;正文 body.pdf 自动 TOC → 合并封面 → `Silent_Acceptance_paper.pdf`。字体注册(Times/SimHei/Consolas 族);设计引擎生成的米色调色板(PAGE_BG #f5f4f3/HEADER_FILL #756b4c/ACCENT #297189);qa_render/ 下 24 张逐页 PNG 是 PDF 的 QA 渲染(视觉验收流程的产物)。
- **`md_to_docx.py`**:Markdown→Word 转换器(标题/粗斜体/行内码/表格/注解块/分隔线/列表),历史版(v1)已被删除,单文件收敛。
- **`make_mindmap.py`**:中文思维导图(中心节点+左右各 3 分支,每分支叶≤4,叶子文案经宽度预算;字号按"图宽 8.6in→PDF 内容宽 6.3in 缩放 0.73"换算)——研究问题/四种受检检验/关键发现(HC3 超尺寸 42-56% 不随 n)/协议承诺/电池发现/结论分支;`mindmap.png` 嵌入 PDF 首页。
- **`battery_extension_section.md`**:v2 新增电池章节的独立 Markdown 片段(合并前的增量草稿)。
- **`figures/full_500_*.png`**:analyze.py 五张图(full tag)+ full_500 五张图(500 种子正式跑)——校准/功效热图/n 扩缩/外推/运行时间。

---

## 8.4 文档层小结

三条产出线的分工与共同纪律:

| 产物 | 讲的故事 | 数字来源 |
|---|---|---|
| Silent_Acceptance_paper.pdf | 现象学基准:指标静默 + 四检验失败画像 + 电池路由 | runs/full(500 种子)+ battery_full |
| HX_...docx | 修复方案:混合方向检验的完整数学规格与随机套件证据 | fc_multidirection/runs(直读 CSV) |
| Falsification_Side_...docx | 两篇合一的统一稿 | 同上两者 |

共同纪律:①论文数字程序化提取(不手抄);②"允许的主张"随证据分层写死;③变更日志/失败案例进正文;④QA 渲染逐页留档。
