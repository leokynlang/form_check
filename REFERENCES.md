# REFERENCES.md — dimreflex 文献对照地图

> 建立于 2026-09-03。目的：标清哪些想法已被先行者覆盖（不可再当卖点），
> 哪些是 dimreflex 的真实增量（论文主张只能落在后者）。
> 每篇按三段写：**他们做了什么 / 与你的设计对照（已覆盖·你的增量·引用方式）**。

---

## 1. Udrescu & Tegmark — AI Feynman（Science Advances 6:eaay2631, 2020；arXiv:1905.11481）

**他们做了什么**：递归符号回归：神经拟合 + 一组"物理启发"分解步骤
——**量纲分析（Buckingham π）是第一步**，其后是广义可分性、乘法可分性、
多项式拟合、暴力搜索。100 道 Feynman Lectures 方程全部发现（当时公开
软件只能解 71 道）；更难测试集从 15% 提到 90%。附带发布
FSReD 基准数据库（6GB，120 题）。

**对照**：
- **已覆盖**：「先用 π 定理降维再在无量纲空间做发现」这条中心路线。
  你独立重建了它（无 AI 背景、未读此文献，独立收敛本身是想法站得住的证据，
  但**不能**再作为新颖性主张）。
- **他们没做（你的增量）**：① π 是他们的一个可选预处理步骤，量纲分析是
  装饰；你的架构里 π 是地基，且配目标隔离、共享变量检测、有理零空间。
  ② 他们输出单一表达式、无降级阶梯、无"数据不足/形式失配显式拒绝"
  的任何机制——论文通篇没有"诚实失败"概念。③ 他们不做工程落地
  （无模板路径、无监控、无不确定度）。
- **引用方式**："AI Feynman applies Buckingham π analysis as one of several
  physics-inspired decomposition steps [Udrescu & Tegmark 2020]. dimreflex
  instead places dimensional analysis at the architectural foundation and
  treats honest failure handling — explicit degradation, refusal on
  insufficient data or mismatched functional form — as a first-class output."

## 2. Bakarji, Callaham, Brunton & Kutz — Dimensionally Consistent Learning with Buckingham Pi（Nature Computational Science 2(12):834–844, 2022；arXiv:2202.04643）

**他们做了什么**：**数据驱动地选 π 基**——利用测量数据的对称/自相似结构，
找"最能坍缩数据"的那组无量纲群（注意：不是任意基，是择优基）。三种方法：
受约束非参拟合、BuckiNet（深度学习第一层投影到 π 空间）、SINDy 版。
案例：旋转环上小球、层流边界层、Rayleigh–Bénard 对流。

**对照**：
- **已覆盖（重要校准）**：你待办 6 里的「π 群基选择」问题，他们已经从
  数据驱动方向做过一版。你的"共享变量检测 + 最小重叠选基"设想不再是
  无人区——**但方向不同**：他们让数据选基，你走"领域知识显式注入"
  （模板路径）+ 检测留痕。这个哲学分歧本身是引用时的定位点。
- **他们没做（你的增量）**：① 无降级阶梯/审计链/监控。② 不做
  "形式合适性"检查（他们只关心发现得出什么，不关心该不该拒）。
  ③ 无测量不确定度传播。
- **引用方式**：在相关工作里引为"自动 π 基选择的代表工作"，然后声明
  设计哲学差异："Where Bakarji et al. let the data select the optimal
  dimensionless basis, dimreflex treats basis selection as domain knowledge
  that must be declared explicitly (template path), with automatic bases
  only flagged, never silently preferred."

## 3. arXiv:2411.15919（2024, NSERC）— Enhancing Symbolic Regression and UPINNs with Dimensional Analysis

**他们做了什么**：离你最近的工程组合：Buckingham π + PySR / AI Feynman
跑代数方程；Ipsen 方法 + UPINNs（universal physics-informed neural nets）
挖微分方程隐藏项。结论与你的体验一致：无量纲化降维→收敛快、精度高、
稀疏数据下更稳。声称首次把 universal differential equations 与
Buckingham π 结合。

**对照**：
- **已覆盖**：「π + 符号回归提升精度/效率」的实证主张——他们已经发表。
  你的正对照（泵）证明的同一类结论不能再当贡献写。
- **他们没做（你的增量）**：全篇是"怎么发现得更好"，没有一行关于
  "发现错了/发现不了时怎么办"：无形式失配检查（你的 form_check 在这一
  光谱上是空的）、无降级、无审计、无监控、无不确定度分离。他们也是
  单赢家输出。
- **引用方式**：作为"π+SR 组合的最近实证"引用，然后转向："These works
  optimize what is discovered; none addresses what the system must do when
  the discovered form is wrong or the data is insufficient — the question
  dimreflex is built around."

## 4. SRBench（La Cava et al., NeurIPS 2021 D&B）与 SRBench++（IEEE Trans. Evol. Comput. 29(4), 2024）

**他们做了什么**：符号回归的活基准。SRBench：25 种 SR + 7 种 ML 方法，
252+ 数据集，测准确率/简洁性/可重发现性。SRBench++：加入领域专家
可解释性评估 + 性质化子任务（噪声鲁棒性、冗余特征、外推、特征选择、
局部极小规避）。2025 年另有"Call for Action"更新版（arXiv:2505.03977）。

**对照**：
- **已覆盖**：命中率/精度/复杂度/外推 这些基准维度。
- **他们没做（你的空位，已两次确认）**：**没有任何一个版本测量
  "系统在该拒答时是否拒答"**——诚实失败率、降级行为、审计链质量
  不在任何指标里。你的负对照方法论（Prandtl–Kármán 台架：判据预注册、
  特异性对照、形式失配证据带）就是这个缺失维度的原型。
- **引用方式**：引为基准标准，然后："SRBench and SRBench++ evaluate what
  algorithms find; no existing benchmark evaluates what algorithms honestly
  report when they should not find — abstention rate, calibrated degradation,
  and audited refusal are unmeasured dimensions."

---

## 位置总结：三个空位（论文主张只能落在这里）

1. **工程诚实作为架构**（不是正则化项）：降级阶梯、双噪声分离
   （回归残差/测量不确定度）、审计链、form_check 逐设计标定的阈值。
   对标物：无。
2. **诚实失败基准**：测拒答率而非命中率。SRBench++ 证明了这一维度
   确实没人测。你的负对照就是它的第一个实例。
3. **Text2π**：量纲检查作为对话式科学推理中最便宜的验证层
   （线性代数即判，无需求解器）。邻近有单位检查器（VerityMath/UnitMath），
   无同形工作。

**写作纪律**（写给未来的自己和任何接手的会话）：
- 绝不再说"我们提出用 π 定理辅助模型发现"——Udrescu/Bakarji/2411.15919 在先。
- 所有新颖性句子必须通过检验：这个具体机制在上面四篇里有没有？
- 引用先行者要大方（他们验证了路线可行），分歧点要精确（他们不回答的问题）。
