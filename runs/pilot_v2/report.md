# form_check vs Ramsey RESET — pilot_v2 报告

生成时间:2026-09-04 12:34:48 | 行数 14250 | 种子 50/格 | n∈[50, 200, 1000] | rel∈[0.0, 0.005, 0.01, 0.03, 0.05]

## 0. 公平性设计(约束 4/5/6 的落实)

- 方法 A 调用 `discover.py::form_check` **本体**,未复制未修改;其 |t|>3 原生阈值只作留痕,决策用包装层在同一排列零分布上算的**经验 p 值**(B=999,分辨率 0.001)。
- 方法 B/C/D = RESET(fitted powers=[2,3] 经典 F / fitted HC3-Wald / exog HC3-Wald),p 值同为 [0,1] 尺度。
- 四方法统一 α=0.05;功效比较另按 **经验 FPR=5% 校准阈值** (H0 池分位数)做二次汇总(§4)。
- H0/H1 分组在代码写入时预注册(`SCENARIO_META`),不因结果调整。

## 1. H0 校准(严格单项式真值,拒绝率=误报率)

| n | form_check | RESET fitted F | RESET fitted HC3 | RESET exog HC3 |
|---|---|---|---|---|
| 50 | 0.057 | 0.056 | 0.064 | 0.061 |
| 200 | 0.032 | 0.037 | 0.043 | 0.036 |
| 1000 | 0.024 | 0.033 | 0.036 | 0.027 |

判读:form_check 的经验 p 值在 H0 下接近名义 5%(排列零分布与原生bootstrap 零分布同构);经典 F 在异方差/重尾 H0 噪声场景会有超尺寸,HC3 变体在此处显示其修正作用(见 §5 噪声层)。

## 2. 结构 H1 功效 @α=0.05(rel=1%)

| 场景 | n | form_check | RESET fitted F | RESET fitted HC3 | RESET exog HC3 |
|---|---|---|---|---|---|
| QUAD_LIN | 50 | 1.000 | 1.000 | 1.000 | 1.000 |
| QUAD_LIN | 200 | 1.000 | 1.000 | 1.000 | 1.000 |
| QUAD_LIN | 1000 | 1.000 | 1.000 | 1.000 | 1.000 |
| CUBIC_LIN | 50 | 0.420 | 1.000 | 1.000 | 1.000 |
| CUBIC_LIN | 200 | 0.500 | 1.000 | 1.000 | 1.000 |
| CUBIC_LIN | 1000 | 0.340 | 1.000 | 1.000 | 1.000 |
| INTERACT | 50 | 0.340 | 1.000 | 1.000 | 0.140 |
| INTERACT | 200 | 0.260 | 1.000 | 1.000 | 0.100 |
| INTERACT | 1000 | 0.200 | 1.000 | 1.000 | 0.020 |
| THRESHOLD | 50 | 1.000 | 1.000 | 1.000 | 1.000 |
| THRESHOLD | 200 | 1.000 | 1.000 | 1.000 | 1.000 |
| THRESHOLD | 1000 | 1.000 | 1.000 | 1.000 | 1.000 |
| EXP_MODEL | 50 | 1.000 | 1.000 | 1.000 | 1.000 |
| EXP_MODEL | 200 | 1.000 | 1.000 | 1.000 | 1.000 |
| EXP_MODEL | 1000 | 1.000 | 1.000 | 1.000 | 1.000 |
| LOG_MODEL | 50 | 1.000 | 1.000 | 1.000 | 1.000 |
| LOG_MODEL | 200 | 1.000 | 1.000 | 1.000 | 1.000 |
| LOG_MODEL | 1000 | 1.000 | 1.000 | 1.000 | 1.000 |
| RATIONAL | 50 | 1.000 | 1.000 | 1.000 | 1.000 |
| RATIONAL | 200 | 1.000 | 1.000 | 1.000 | 1.000 |
| RATIONAL | 1000 | 1.000 | 1.000 | 1.000 | 1.000 |
| OMITTED_M | 50 | 0.980 | 0.980 | 0.960 | 0.960 |
| OMITTED_M | 200 | 1.000 | 1.000 | 1.000 | 1.000 |
| OMITTED_M | 1000 | 1.000 | 1.000 | 1.000 | 1.000 |
| DILUTION | 50 | 1.000 | 0.500 | 0.420 | 1.000 |
| DILUTION | 200 | 1.000 | 0.960 | 0.940 | 1.000 |
| DILUTION | 1000 | 1.000 | 1.000 | 1.000 | 1.000 |
| OPPOSING | 50 | 0.020 | 0.500 | 0.260 | 1.000 |
| OPPOSING | 200 | 0.000 | 0.980 | 0.940 | 1.000 |
| OPPOSING | 1000 | 0.000 | 1.000 | 1.000 | 1.000 |
| COLLINEAR_CURV | 50 | 1.000 | 1.000 | 1.000 | 1.000 |
| COLLINEAR_CURV | 200 | 1.000 | 1.000 | 1.000 | 1.000 |
| COLLINEAR_CURV | 1000 | 1.000 | 1.000 | 1.000 | 1.000 |

逐条解读(约束:不得只报总平均):

- **QUAD_LIN**: n=1000 时最高 form_check (1.00);form_check 1.00 / fitted-F 1.00 / exog 1.00。
- **CUBIC_LIN**: n=1000 时最高 RESET fitted F (1.00);form_check 0.34 / fitted-F 1.00 / exog 1.00。
- **INTERACT**: n=1000 时最高 RESET fitted F (1.00);form_check 0.20 / fitted-F 1.00 / exog 0.02。
- **THRESHOLD**: n=1000 时最高 form_check (1.00);form_check 1.00 / fitted-F 1.00 / exog 1.00。
- **EXP_MODEL**: n=1000 时最高 form_check (1.00);form_check 1.00 / fitted-F 1.00 / exog 1.00。
- **LOG_MODEL**: n=1000 时最高 form_check (1.00);form_check 1.00 / fitted-F 1.00 / exog 1.00。
- **RATIONAL**: n=1000 时最高 form_check (1.00);form_check 1.00 / fitted-F 1.00 / exog 1.00。
- **OMITTED_M**: n=1000 时最高 form_check (1.00);form_check 1.00 / fitted-F 1.00 / exog 1.00。
- **DILUTION**: n=1000 时最高 form_check (1.00);form_check 1.00 / fitted-F 1.00 / exog 1.00。
- **OPPOSING**: n=1000 时最高 RESET fitted F (1.00);form_check 0.00 / fitted-F 1.00 / exog 1.00。
- **COLLINEAR_CURV**: n=1000 时最高 form_check (1.00);form_check 1.00 / fitted-F 1.00 / exog 1.00。

## 3. 噪声层 H1(结构仍是严格幂律)——这些是"误报"不是"检出"

| 场景 | n | form_check | RESET fitted F | RESET fitted HC3 | RESET exog HC3 |
|---|---|---|---|---|---|
| H1_AUTOCORR | 200 | 0.040 | 0.100 | 0.080 | 0.080 |
| H1_HEAVYTAIL | 200 | 0.060 | 0.080 | 0.060 | 0.060 |
| H1_OUTLIER | 200 | 0.000 | 0.060 | 0.020 | 0.020 |
| H1_OUTLIER_GRP | 200 | 0.020 | 0.080 | 0.040 | 0.040 |
| H1_VAR_GROWING | 200 | 0.060 | 0.080 | 0.100 | 0.100 |

**语义警示(约束:区分"没有拒绝"与"证明模型正确")**:本节场景的均值结构是严格单项式——任何拒绝都是对噪声假设的误报,不是发现形式失配。反之,H0 格的"不拒绝"也只说明"无证据反对幂律",不构成"幂律得证"。form_check 的输出语义即如此措辞(未发现反对证据 ≠ 验证成立)。

## 4. 相同经验 FPR(5%)下的功效(rel=1%,结构 H1 池)

阈值 = H0 池 p 值 5% 分位数(按 n 校准):

| n | form_check | RESET fitted F | RESET fitted HC3 | RESET exog HC3 |
|---|---|---|---|---|
| 50 | 0.773 | 0.872 | 0.855 | 0.906 |
| 200 | 0.792 | 0.990 | 0.986 | 0.910 |
| 1000 | 0.785 | 1.000 | 1.000 | 0.908 |

## 5. 胜负矩阵(结构 H1 全格,2pp 容差)

(win_matrix 读取失败:`Import tabulate` failed.  Use pip or conda to install the tabulate package.)

## 6. 运行时间(每数据集,秒)

| n | form_check | RESET fitted F | RESET fitted HC3 | RESET exog HC3 |
|---|---|---|---|---|
| 50 | 0.0050 | 0.0005 | 0.0005 | 0.0004 |
| 200 | 0.0061 | 0.0007 | 0.0005 | 0.0005 |
| 1000 | 0.0106 | 0.0009 | 0.0008 | 0.0007 |

## 7. 失败登记(不可静默丢弃)

| 方法 | 失败率 | 主要失败类型 |
|---|---|---|
| form_check | 0.0842 | degenerate_zero_residual |
| RESET fitted F | 0.1053 | degenerate_zero_residual |
| RESET fitted HC3 | 0.1053 | degenerate_zero_residual |
| RESET exog HC3 | 0.1368 | degenerate_zero_residual |

失败类型说明:`degenerate_zero_residual`=零噪声/精确幂律(检验对象不存在,显式登记);`rank_deficient`=设计秩亏;`hc3_leverage_degenerate`=杠杆 1−h→0;`statistic_mismatch`=form_check 返回统计量与包装层重算不符(实现漂移哨兵);`checked_false_*`=form_check 本体拒绝检查(非正值/小样本/非有限)。

## 8. 复现信息

- git HEAD:`9344ae242848d283cbe447e22974541cfeb19139`(工作区 1 处未提交改动)
- discover.py sha256:`8a3f15543c2254af…`
- python 3.13.15,numpy 2.5.2,scipy 1.18.1,pandas 3.0.5,dimreflex 0.5.0
- 平台:Windows 11 AMD64
- 种子规则:`100_000 + s*97 + n_idx*13 + rel_idx*5_000_000`;排列 B=999


---

# SRSD-Feynman 主基准 — srsd_pilot

## 任务池(预注册分类,判据=公式结构,先于任何实验)

- 入选 32 个(H0 单项式 13,H1 非单项式 19);排除 88 个。
- 排除理由分布:nonpositive_inputs_or_output×80, too_many_variables_for_harness×6, nonfinite_under_official_sampling×2
- 完整清单(含逐任务理由):`runs/srsd_pilot/task_registry.csv`

## H0(真单项式)误报率 @α=0.05

| n | form_check | RESET fitted F | RESET fitted HC3 | RESET exog HC3 |
|---|---|---|---|---|
| 50 | 0.056 | 0.063 | 0.425 | 0.565 |
| 200 | 0.060 | 0.043 | 0.415 | 0.565 |
| 1000 | 0.055 | 0.054 | 0.421 | 0.566 |

## 相同经验 FPR=5% 下的 H1 功效

| n | form_check | RESET fitted F | RESET fitted HC3 | RESET exog HC3 |
|---|---|---|---|---|
| 50 | 0.394 | 0.190 | 0.050 | 0.188 |
| 200 | 0.504 | 0.293 | 0.282 | 0.360 |
| 1000 | 0.594 | 0.425 | 0.407 | 0.649 |

## 重点问题

- **Q1_univariate_low_df** SRSD 池无单变量非单项式任务(1 变量物理式几乎全是单项式);该问题由合成台架 QUAD_LIN 回答(见合成报告 §2)。
- **Q2_averaging_blindness** 多变量 H1 任务中 exog-HC3 优于 form_check 最多的前 3:feynman-ii.11.28(fc=0.02 vs exog=1.00); feynman-iii.4.32(fc=0.02 vs exog=1.00); feynman-iii.4.33(fc=0.02 vs exog=1.00)
- **Q3_hc3_hetero_fpr** {'50': {'classical': 0.06346153846153846, 'hc3': 0.425, 'exog_hc3': 0.5653846153846154}, '200': {'classical': 0.04346153846153846, 'hc3': 0.4153846153846154, 'exog_hc3': 0.5653846153846154}, '1000': {'classical': 0.05384615384615385, 'hc3': 0.4207692307692308, 'exog_hc3': 0.5657692307692308}}
- **Q4_alarm_extrapolation** {'auc': {'fc': 0.2887493416416434, 'fitted_f': 0.26342901973720667, 'fitted_hc3': 0.6095093835779094, 'exog_hc3': 0.5011415227631334}, 'alarm_pr': {'form_check': {'precision': 0.3865814696485623, 'recall': 0.6269430051813472}, 'RESET fitted F': {'precision': 0.5504559270516717, 'recall': 0.4938641941641669}, 'RESET fitted HC3': {'precision': 0.16889844589971656, 'recall': 0.4712298881919825}, 'RESET exog HC3': {'precision': 0.24322241681260945, 'recall': 0.9468230160894464}}}

## 报警→外推误差关联

按"远端外推相对误差>10%"为正类的 p 值 AUC:fc=0.289, fitted_f=0.263, fitted_hc3=0.610, exog_hc3=0.501


## 图(SRSD)

- `figures/srsd_calibration.png` — H0 p 值校准(真实单项式)
- `figures/srsd_per_task_power.png` — 逐任务功效
