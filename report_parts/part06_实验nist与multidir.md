# 第六部分 实验二与实验三:NIST 数值认证与多方向扩展(`experiments/nist_strd/`、`experiments/fc_multidirection/`)

> 两个实验目录是 form_check vs RESET 主基准(P1-P3)的**外部效度与盲区修复**延伸:
> - `nist_strd/`(5 个 Python 文件,约 1,300 行):用 NIST 统计参考数据集(StRD)的**认证值**验证整个比较框架的数值正确性,并把"检测能力问题"与"数值计算问题"二分归因。
> - `fc_multidirection/`(6 个文件,约 1,500 行):把 form_check 的"单一中心二次方向"机制扩展为**多方向 + RESET 混合**检验(hx),攻击已知盲区(OPPOSING/CUBIC/THRESHOLD),并用**随机真值泛化测试**验证设计视野外的行为。
>
> 两个目录都有强烈的"列名隔离"纪律:`p_fc` 是产品成绩,`p_qda`/`p_md*`/`p_hx` 是实验 analogue,docstring 反复声明"任何报告中不得把 analogue 结果表述为产品 form_check 的成绩"。

---

## 6.1 `experiments/nist_strd/` — NIST StRD 认证验证

### 6.1.1 `parse_nist.py` — 数据解析器与双通道自检(204 行)

**纪律**(docstring):只读 data/ 目录;解析结果与 EXPECTED(计划阶段从 NIST 网页独立摘录的第二通道)逐项断言一致才放行;断言失败立即抛错,"绝不带病放行"。

- `parse_dat(path)`(parse_nist.py:40-99):正则抽取数据集名/观测数/参数数、认证系数与标准误(`B(\d+)` 行,Fortran 记数 `D`→`E` 转换)、残差标准差、R²、ANOVA 表、数据块(形状断言 `(n_obs, 1+n_pred)`)。
- **双通道裁决记录**(parse_nist.py:23-32 注释,实验方法学的亮点):Longley B5 认证标准误在 NIST 网页摘录通道为 0.0226,.dat 文件为 0.226,两通道不一致。裁决**以 .dat 文件为准**,证据三条:(a) 三条独立 SE 路径(QR/normal-eq/scaled-cov)从真实数据复算一致得 0.2260;(b) R/statsmodels 公认输出 t=−0.226 与 0.226 组出自洽;(c) t 一致性:取 0.0226 会得 t=−2.26 与公认输出不符。**原因**:数据 provenance 冲突不掩盖、不留含糊,裁决与证据留痕。
- `self_check()`(parse_nist.py:153-185):sha256 对账(provenance.json)+ 全部认证值双通道断言(relerr<1e-13)+ **认证值内部一致性**(`s² == MS_resid`、`SST == SSR+SSE`、`R² == SSR/SST`、`F == MS_reg/MS_res`,1e-12)。**原因**:连 NIST 自己的数都做代数自洽校验——解析 bug 在进入任何实验之前就死掉。
- EXPECTED 字典(parse_nist.py:105-144):四数据集(Norris 36×2、Pontius 40×3、Filip 82×11、Longley 16×7)完整认证值硬编码。

### 6.1.2 `analogue.py` — 原空间"二次方向" analogue(112 行)

- **声明**(docstring):"本模块不是产品代码,不是 dimreflex.form_check 的替身、移植版或成绩"——form_check 定义在 log-log 空间(要求全正值),Filip 等数据 x 全为负、产品按结构化拒绝。本模块把"中心二次方向列 + 残差置换零分布"机制**镜像到任意原空间回归**,只回答"检验机制本身在认证回归世界里行为是否正常"。列名 `p_qda/rej_qda/t_qda` 与产品 `p_fc` 严格分开。
- `quad_direction_analogue(X_raw, y_raw)`(analogue.py:33-97):z 标准化特征平方的行均值作 quad 列(注释解释:"log 在 form_check 里同时承担'幂律基准'与'尺度压缩'两个角色;原空间没有幂律语义,尺度压缩角色由 z 标准化承担");设计 `[1, X_raw, quad]`;pinv 末行向量化置换零分布(与 common.py 同技巧);失败四类写 fail 字段不抛异常。另有 `__main__` 冒烟:线性真值应不告警、二次真值应告警。
- **导向**:`run_nist.detection_layer` 的 qda 列;`analogue_calibration`(H0 校准抽查,n=36/60 各 200 次线性真值,报告留痕)。

### 6.1.3 `run_nist.py` — 三层检验主入口(537 行)

**八规格设计**(run_nist.py:45-62):

| 规格 | 数据集 | 模型 | 角色 |
|---|---|---|---|
| N1 | Norris | 1 次多项式 | certified_correct(误报考察主战场) |
| P1 | Pontius | 1 次 | misspecified_omit_x²(应检出) |
| P2 | Pontius | 2 次 | certified_correct |
| F1-F3 | Filip | 1/2/3 次 | misspecified(漏 2-10 阶) |
| F10 | Filip | 10 次 | certified_correct(数值压力关键点) |
| L1 | Longley | 6 预测元线性 | certified_correct(共线性考察) |

- **`fit_paths(X, y)`**(run_nist.py:115-165)——四求解路径互检:lstsq/pinv/qr/scaled(中心化标准化后求解再仿射映射回原系数,协方差经雅可比 `J = [[1,−μ/σ],[0,diag(1/σ)]]` 变换)。三条 SE 路径对照:qr(数值较稳)、normal_eq(条件数平方——"经典翻车路径,预期在 Pontius/Filip raw 上退化,如实展示")、scaled_cov。诊断记录 `_cond/_rank/_rank_default(numpy 默认容差判定)/_warn_all`。**原因**:`_stats` 注释写明"NIST 认证规格数学上满秩;numpy 默认 matrix_rank 在病态设计上会降秩,那是数值现象,单独记入诊断,不用来定 df/s/F"——自由度按设计列数算,不按数值秩。
- **`certified_check`**(run_nist.py:182-210):逐系数相对误差(近零系数以自身认证 se 为尺度——"有统计意义的下限")、s/SSE/SSR/F/R² 相对认证值偏差;**框架正确性以 scaled 路径为准**(阈值 Norris/Pontius/Longley<1e-8,Filip<1e-6——"认证值 500 位精度算出,双精度可达位数按 1e-6 判");raw 路径偏差用于展示设计尺度数值退化(NIST 本意),不计为框架错误。
- **`detection_layer`**(run_nist.py:222-255):RESET 三变体(raw 设计)+ qda;**raw → scaled(Z) → full(Z,yz) 三级翻转对照**——"拟合值参数化不变 ⇒ 统计量理论不变 ⇒ 翻转即数值伪影"。full 级的必要性注释:"Longley y~6e4 时 yhat³~3e14 直接引爆增广设计秩判定,Z-scaled 仍救不回,full 才能运行"。
- **`effective_decision` + `classify` + `numerical_reasons`**(run_nist.py:283-356)——二分归因规则:
  - scaled 认证核对超阈 ⇒ `FRAMEWORK_BUG(立即停下报告)`;
  - 数值红旗仅在 raw 尺度而 scaled/full 干净 ⇒ 数值计算问题(设计尺度),不影响框架判定;
  - 数值全绿仍告警/不告警 ⇒ 检测层现象(能力或语义)。
  - 数值红旗收集器:raw/scaled 决策翻转、exog 失败两分类("结构性(幂列已在基设计中,变体定义如此)" vs "数值性(原始尺度幂列引爆秩判定)")、HC3 scaled/full p 相差 >6 个数量级("pinv(X'X) 条件数平方路径失效"哨兵)、numpy 默认秩降秩、求解警告。
- **`fc_applicability`**(run_nist.py:259-279):产品 form_check 在四数据集上的适用性矩阵,解读逐条写死:Norris/Pontius "适用(全正),认证真值是原空间线性/二次非幂律——log-log 告警属'形式层面检出',不是同规格 RESET 意义的误报";Filip "不适用:x 全为负,产品按结构化拒绝(预期行为)";Longley "适用,6 特征走多列平均二次方向分支"。
- **`analogue_calibration`**(run_nist.py:359-370):固定种子线性真值 200 次重复,H0 拒绝率写入报告供复现。
- **`write_report`**(run_nist.py:388-479):报告开头声明 qda 列不得表述为产品成绩;五节(适用性矩阵/规格表/统计量与 p 值/认证偏差/归因结论)+ 双通道不一致裁决留痕 + analogue 校准表("轻度自由倾向与产品 form_check 在 log 域的实测校准 5.7/3.2/2.4% 同型,属机制属性")。
- **【导向】** `runs/nist_strd_results.csv / certified_check.csv / fc_applicability.csv / env_snapshot.json / report.md`。**实验结论**(报告与 README):框架 scaled 路径全部通过认证核对(无 FRAMEWORK_BUG);Pontius/Filip 遗漏规格被 RESET/qda 检出;认证正确规格上的少量告警被归因为"检验零分布假设与数据噪声结构的语义现象"而非错误;Longley 是 n=16 边界情形("任何单一结论都不稳")。

### 6.1.4 `compare_reset_vs_fc.py` — 同语义/跨语义两层对照(171 行)

- **【计算】** 两层比较:
  - **同语义层**(log-log 世界):`form_check(p_fc)` vs **log 空间 RESET**(`loglog_reset`:log y ~ [1, log X] + yhat²,yhat³)——同一基准模型(幂律)、两种检验机制,回答"机制是否一致"。
  - **跨语义层**:raw 空间认证规格 RESET vs form_check——回答"语义不同时结论为何分歧"。
  - **Norris 噪声掩盖诊断**(`norris_noise_diagnostic`):认证噪声 s≈0.885 为绝对量,y<1 区间(仅 6 个点)相对噪声 332%,信号被淹没——fc(0.789)与 log-RESET(0.334)双静默是"检测功效受限,非误报,两机制一致"的机理解释。
- **【产出】** `runs/compare_reset_vs_fc.md`,含"机制属性对照"总表(统计量/基准模型/H0 校准/结构功效盲区/反向盲区/HC3 变体/数值稳健/决策阈值/NIST 边界行为 九行对照)。关键结论:fc 与 log-RESET 经典 F 在 Norris/Pontius 决策一致;Longley 同语义下三零分布分裂(fc 置换 0.004 告警/经典 F 0.396 静默/HC3 4.5e-29 告警)——"边界情形,归为检测层方向/零分布敏感性,非数值错误"。

---

## 6.2 `experiments/fc_multidirection/` — 多方向扩展与混合检验

### 6.2.1 `multidir.py` — 多方向 analogue 探索台(219 行)

- **【计算】** 方向池(`build_directions`,每列中心化,z 为标准化特征):`quad_j = z_j²`、`cubic_j = z_j³`、`inter_ij = z_i·z_j`、`osc_j = sin(π·z_j)`,可选 `yhat^k`(RESET 自引用)。算法四步(multidir.py:102-192):
  1. **贪心枢轴选择**:每列须对 [基∪已选方向] 添加新维度才保留(按最大残差范数;注释:"方向列可'相对基共线',逐列阈值不够");预算 = min(max_dirs, n − k_base − df_floor) 保证高维短样本可运行;
  2. **修正 Gram-Schmidt 正交化**:保留列对基与彼此做 GS → W(单位正交,⊥基)。动机:"近平行的高次代理(cubic_x2∝x⁶ 与 yhat² 等)在 min-norm 求解下互相分摊信号,逐方向 t 被稀释;正交坐标下归因稳定";
  3. **双统计量**(同一批置换零分布,受限残差置换;W⊥基 ⇒ RSS_r 置换不变、RSS_u = RSS_r − Σd²,全向量化):maxT = max_j|t_j|(抓集中信号)+ joint F = (Σd²/m)/((RSS_r−Σd²)/(n−k))(聚合分散信号);
  4. **合并 p = min(1, 2·min(p_maxT, p_F))**(Bonferroni,H0 下保守有效)。
- 返回携带逐方向 `t_dir/p_dir/argmax_dir/kept_names/n_dir_dropped`——**归因**(哪个方向抓住信号)是本扩展的独有能力。
- 失败七类写 fail 不抛(含 `all_directions_in_span`)。
- `__main__` 冒烟:四案例(纯幂律静默/二次/三次/震荡应告警)。

### 6.2.2 `hybrid.py` — RESET × 多方向混合检验(226 行)

- **出身**(docstring):"md 方向族抓反向/稀释/交互/震荡曲率(RESET 的盲区),yhat 幂列自动可达任意高阶代理(md 固定阶梯的盲区,如 Filip F3 的 4-10 阶残差)"。`hybrid_test(Xlog, ly, kinds=(quad,cubic,interact,osc), extra_fitted=(2,3))` 与 multidir 的算法骨架相同,默认池含 yhat²,yhat³。
- 同样有"非产品"声明与独立失败分类。`__main__` 冒烟含"四次真值(应告警,归因 yhat)"案例。

### 6.2.3 `run_bench.py` — 场景基准(261 行)

- **【计算】** 问题:"form_check 能否检测高次方程/震荡曲率?多方向扩展能否补盲区而不破坏 H0 校准?" 五方法并列:`p_fc`(产品本体)/`p_md3`(quad+cubic+interact)/`p_md4`(+osc)/`p_hx`(混合)/`p_reset`(参照)。场景**只读复用** form_check_vs_reset/scenarios.py("预注册分类不改"):H0 校准 3 场景 + H1 功效 6 场景(QUAD/CUBIC/INTERACT/DILUTION/OPPOSING/THRESHOLD)。网格 n∈{50,200},rel=0.01,50 seeds;**CUBIC_LIN/OPPOSING 加 n=1000**("文档记载 fc 功效随 n 反降,验证 md 是否免疫")。种子规则与主基准同一族 `100_000 + s·97 + n_idx·13`。
- 汇总:逐场景×n 的拒绝率+Wilson CI(bench_summary.csv);**argmax 方向分布诊断**(argmax_directions.csv——CUBIC_LIN 的信号被 cubic_x1 抓住等归因证据);NIST 真实数据 log 空间抽查(Norris 告警的注释:"信号源自噪声主导区,稳健性存疑,属形式层面检出";Pontius 与 fc/raw RESET 一致检出);报告五节含逐场景数字结论。
- **【结论】**(runs/report.md):多方向在 H0 校准保持下补掉 OPPOSING/CUBIC 盲区,THRESHOLD 上所有方法都弱(方向阶梯不吸收折角)。

### 6.2.4 `random_ten.py` — 随机真值泛化测试(288 行)

- **目的**(docstring):"预注册场景之外的随机漂移泛化检验(CONTEXT.md 建议的'设计视野外漂移')"。
- **【计算】** 主 RNG 独立于场景库:每组随机抽类别(H0:H1=4:6)、真值家族(12 种 H1 家族:quad/cubic/quartic/interact/cross_cubic/sin1/sin2/exp_deform/threshold/logistic_bump/opposing_quad/mixed_cubic_osc)、参数、n∈[60,400]、rel∈[0.5%,3%]、特征数 1-2;**非幂律成分按随机目标信噪比缩放 `max|g|/rel ∈ [3,12]`**("非平凡可检");四方法盲跑(fc/md4/hx/reset),期望"H0 全体静默、H1 全体告警"。N≤10 走原始口径(存档可复现),N>10 通用路径(种子 MASTER_SEED+n_groups 保证与 10 组版流不同)。
- 汇总:FPR/TPR(Wilson CI)、**H1 逐家族检出率**、配对分歧(检出集合差)、总符合率、hx 归因分布(top-8 argmax)。FAIL 计入不符合;注释:"单组二值方差大,本测试为泛化检验,预注册场景功效见 run_bench"。
- **【结论】**(runs/random_1000.md):1000 组泛化上 hx 总符合率最高、fc 在 opposing_quad/高次家族漏检率可见——支持"多方向扩展修复盲区"主结论。

### 6.2.5 `compare_vs_reset.py` — 三块对比(200 行)

- **A 合成基准逐种子配对一致性**(不重跑 bench_raw.csv):每场景×n 统计双拒绝/md4 独拒/reset 独拒/双接受——"比边际率更锐利"(配对设计消掉样本间变异)。
- **B NIST 认证数据交叉对比**:md4/hx 原空间运行 Filip(fc 结构上不可)与 RESET 同台;Longley n=16 全方向集(33+7 列)超 df → **quad-only 降级**(维度扩展性对比点)。RESET p 值取 raw 级、失败回退 full 级的映射逻辑。结果解读内嵌:F3 已由混合修复("md4 单独时方向阶梯对基共线被剔,够不到 4-10 阶残差;hx 的 yhat²,yhat³ 自引用列补上,argmax=cubic_x3 即 x⁹ 代理,p=0.002,与 RESET 8.7e-26 同判");F10 仍不可运行("需正交多项式基,未来工作");L1 边界("残差 df=3、cond 4.9e9,单数据集不可归因")。
- **C 单次调用成本微基准**:n=200/p=2/30 次中位毫秒,md4 vs hx vs RESET(置换 999 次 vs 渐近 F,慢 ~N×,均毫秒级)。
- **【导向】** `runs/compare_vs_reset.md + compare_nist_md4.csv`;连同 bench/random 报告构成 fc_multidirection/runs/report.md 的完整证据链。

---

## 6.3 本部分小结

两个延伸实验共同回答三个元问题:

1. **"比较框架本身可信吗?"** —— NIST 认证核对层:解析器双通道断言 + 四路径互检 + scaled 路径全绿,把"检测结论"与"数值伪影"用三级翻转对照(raw→scaled→full)严格分离。
2. **"已知盲区可修吗?"** —— 多方向 + 正交化 + maxT/joint-F 双统计量 + Bonferroni 合并,在 H0 校准保持下补掉 OPPOSING/CUBIC;混合 yhat 列再补高阶代理盲区(Filip F3)。
3. **"修复泛化吗?"** —— 1000 组随机真值(12 家族、随机信噪比)盲跑验证,期望符合率 + 配对分歧 + 归因分布三视角。

同时两个目录都以 docstring 级纪律守护"产品成绩"与"实验 analogue"的边界——这是论文 claim 写作不越界(只把 p_fc 写成 dimreflex 成绩)的代码级保障。
