# form_check vs Ramsey RESET — full 报告

生成时间:2026-09-04 20:18:44 | 行数 142500 | 种子 500/格 | n∈[50, 200, 1000] | rel∈[0.0, 0.005, 0.01, 0.03, 0.05]

## 0. 公平性设计(约束 4/5/6 的落实)

- 方法 A 调用 `discover.py::form_check` **本体**,未复制未修改;其 |t|>3 原生阈值只作留痕,决策用包装层在同一排列零分布上算的**经验 p 值**(B=999,分辨率 0.001)。
- 方法 B/C/D = RESET(fitted powers=[2,3] 经典 F / fitted HC3-Wald / exog HC3-Wald),p 值同为 [0,1] 尺度。
- 四方法统一 α=0.05;功效比较另按 **经验 FPR=5% 校准阈值** (H0 池分位数)做二次汇总(§4)。
- H0/H1 分组在代码写入时预注册(`SCENARIO_META`),不因结果调整。

## 1. H0 校准(严格单项式真值,拒绝率=误报率)

| n | form_check | RESET fitted F | RESET fitted HC3 | RESET exog HC3 |
|---|---|---|---|---|
| 50 | 0.047 | 0.045 | 0.048 | 0.054 |
| 200 | 0.037 | 0.041 | 0.045 | 0.042 |
| 1000 | 0.041 | 0.038 | 0.038 | 0.040 |

判读:form_check 的经验 p 值在 H0 下接近名义 5%(排列零分布与原生bootstrap 零分布同构);经典 F 在异方差/重尾 H0 噪声场景会有超尺寸,HC3 变体在此处显示其修正作用(见 §5 噪声层)。

## 2. 结构 H1 功效 @α=0.05(rel=1%)

| 场景 | n | form_check | RESET fitted F | RESET fitted HC3 | RESET exog HC3 |
|---|---|---|---|---|---|
| QUAD_LIN | 50 | 1.000 | 1.000 | 1.000 | 1.000 |
| QUAD_LIN | 200 | 1.000 | 1.000 | 1.000 | 1.000 |
| QUAD_LIN | 1000 | 1.000 | 1.000 | 1.000 | 1.000 |
| CUBIC_LIN | 50 | 0.450 | 1.000 | 1.000 | 1.000 |
| CUBIC_LIN | 200 | 0.404 | 1.000 | 1.000 | 1.000 |
| CUBIC_LIN | 1000 | 0.470 | 1.000 | 1.000 | 1.000 |
| INTERACT | 50 | 0.228 | 1.000 | 1.000 | 0.172 |
| INTERACT | 200 | 0.212 | 1.000 | 1.000 | 0.098 |
| INTERACT | 1000 | 0.204 | 1.000 | 1.000 | 0.074 |
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
| OMITTED_M | 50 | 0.942 | 0.918 | 0.958 | 0.958 |
| OMITTED_M | 200 | 1.000 | 1.000 | 1.000 | 1.000 |
| OMITTED_M | 1000 | 1.000 | 1.000 | 1.000 | 1.000 |
| DILUTION | 50 | 1.000 | 0.330 | 0.294 | 1.000 |
| DILUTION | 200 | 1.000 | 0.880 | 0.866 | 1.000 |
| DILUTION | 1000 | 1.000 | 1.000 | 1.000 | 1.000 |
| OPPOSING | 50 | 0.016 | 0.510 | 0.310 | 1.000 |
| OPPOSING | 200 | 0.018 | 0.980 | 0.976 | 1.000 |
| OPPOSING | 1000 | 0.012 | 1.000 | 1.000 | 1.000 |
| COLLINEAR_CURV | 50 | 1.000 | 1.000 | 1.000 | 1.000 |
| COLLINEAR_CURV | 200 | 1.000 | 1.000 | 1.000 | 1.000 |
| COLLINEAR_CURV | 1000 | 1.000 | 1.000 | 1.000 | 1.000 |

逐条解读(约束:不得只报总平均):

- **QUAD_LIN**: n=1000 时最高 form_check (1.00);form_check 1.00 / fitted-F 1.00 / exog 1.00。
- **CUBIC_LIN**: n=1000 时最高 RESET fitted F (1.00);form_check 0.47 / fitted-F 1.00 / exog 1.00。
- **INTERACT**: n=1000 时最高 RESET fitted F (1.00);form_check 0.20 / fitted-F 1.00 / exog 0.07。
- **THRESHOLD**: n=1000 时最高 form_check (1.00);form_check 1.00 / fitted-F 1.00 / exog 1.00。
- **EXP_MODEL**: n=1000 时最高 form_check (1.00);form_check 1.00 / fitted-F 1.00 / exog 1.00。
- **LOG_MODEL**: n=1000 时最高 form_check (1.00);form_check 1.00 / fitted-F 1.00 / exog 1.00。
- **RATIONAL**: n=1000 时最高 form_check (1.00);form_check 1.00 / fitted-F 1.00 / exog 1.00。
- **OMITTED_M**: n=1000 时最高 form_check (1.00);form_check 1.00 / fitted-F 1.00 / exog 1.00。
- **DILUTION**: n=1000 时最高 form_check (1.00);form_check 1.00 / fitted-F 1.00 / exog 1.00。
- **OPPOSING**: n=1000 时最高 RESET fitted F (1.00);form_check 0.01 / fitted-F 1.00 / exog 1.00。
- **COLLINEAR_CURV**: n=1000 时最高 form_check (1.00);form_check 1.00 / fitted-F 1.00 / exog 1.00。

## 3. 噪声层 H1(结构仍是严格幂律)——这些是"误报"不是"检出"

| 场景 | n | form_check | RESET fitted F | RESET fitted HC3 | RESET exog HC3 |
|---|---|---|---|---|---|
| H1_AUTOCORR | 200 | 0.038 | 0.040 | 0.042 | 0.042 |
| H1_HEAVYTAIL | 200 | 0.034 | 0.036 | 0.030 | 0.030 |
| H1_OUTLIER | 200 | 0.040 | 0.066 | 0.036 | 0.036 |
| H1_OUTLIER_GRP | 200 | 0.038 | 0.060 | 0.054 | 0.054 |
| H1_VAR_GROWING | 200 | 0.048 | 0.058 | 0.050 | 0.050 |

**语义警示(约束:区分"没有拒绝"与"证明模型正确")**:本节场景的均值结构是严格单项式——任何拒绝都是对噪声假设的误报,不是发现形式失配。反之,H0 格的"不拒绝"也只说明"无证据反对幂律",不构成"幂律得证"。form_check 的输出语义即如此措辞(未发现反对证据 ≠ 验证成立)。

## 4. 相同经验 FPR(5%)下的功效(rel=1%,结构 H1 池)

阈值 = H0 池 p 值 5% 分位数(按 n 校准):

| n | form_check | RESET fitted F | RESET fitted HC3 | RESET exog HC3 |
|---|---|---|---|---|
| 50 | 0.777 | 0.879 | 0.862 | 0.908 |
| 200 | 0.786 | 0.987 | 0.985 | 0.908 |
| 1000 | 0.782 | 1.000 | 1.000 | 0.906 |

## 5. 胜负矩阵(结构 H1 全格,2pp 容差)

|                  |   form_check |   RESET fitted F |   RESET fitted HC3 |   RESET exog HC3 |
|:-----------------|-------------:|-----------------:|-------------------:|-----------------:|
| form_check       |            0 |               10 |                  8 |               12 |
| RESET fitted F   |           36 |                0 |                  7 |               12 |
| RESET fitted HC3 |           37 |                4 |                  0 |               12 |
| RESET exog HC3   |           25 |               17 |                 14 |                0 |

## 6. 运行时间(每数据集,秒)

| n | form_check | RESET fitted F | RESET fitted HC3 | RESET exog HC3 |
|---|---|---|---|---|
| 50 | 0.0120 | 0.0010 | 0.0011 | 0.0010 |
| 200 | 0.0144 | 0.0014 | 0.0013 | 0.0012 |
| 1000 | 0.0253 | 0.0018 | 0.0017 | 0.0017 |

## 7. 失败登记(不可静默丢弃)

| 方法 | 失败率 | 主要失败类型 |
|---|---|---|
| form_check | 0.0842 | degenerate_zero_residual |
| RESET fitted F | 0.1053 | degenerate_zero_residual |
| RESET fitted HC3 | 0.1053 | degenerate_zero_residual |
| RESET exog HC3 | 0.1368 | degenerate_zero_residual |

失败类型说明:`degenerate_zero_residual`=零噪声/精确幂律(检验对象不存在,显式登记);`rank_deficient`=设计秩亏;`hc3_leverage_degenerate`=杠杆 1−h→0;`statistic_mismatch`=form_check 返回统计量与包装层重算不符(实现漂移哨兵);`checked_false_*`=form_check 本体拒绝检查(非正值/小样本/非有限)。

## 8. 复现信息

- git HEAD:`9344ae242848d283cbe447e22974541cfeb19139`(工作区 46 处未提交改动)
- discover.py sha256:`bdf7b06978f2e1d3…`
- python 3.13.15,numpy 2.5.2,scipy 1.18.1,pandas 3.0.5,dimreflex 0.5.0
- 平台:Windows 11 AMD64
- 种子规则:`100_000 + s*97 + n_idx*13 + rel_idx*5_000_000`;排列 B=999


---

# SRSD-Feynman 主基准 — srsd_full

## 任务池(预注册分类,判据=公式结构,先于任何实验)

- 入选 32 个(H0 单项式 13,H1 非单项式 19);排除 88 个。
- 排除理由分布:nonpositive_inputs_or_output×80, too_many_variables_for_harness×6, nonfinite_under_official_sampling×2
- 完整清单(含逐任务理由):`runs/srsd_full/task_registry.csv`

## H0(真单项式)误报率 @α=0.05

| n | form_check | RESET fitted F | RESET fitted HC3 | RESET exog HC3 |
|---|---|---|---|---|
| 50 | 0.049 | 0.050 | 0.422 | 0.564 |
| 200 | 0.053 | 0.050 | 0.417 | 0.563 |
| 1000 | 0.052 | 0.053 | 0.418 | 0.565 |

## 相同经验 FPR=5% 下的 H1 功效

| n | form_check | RESET fitted F | RESET fitted HC3 | RESET exog HC3 |
|---|---|---|---|---|
| 50 | 0.409 | 0.194 | 0.050 | 0.186 |
| 200 | 0.513 | 0.281 | 0.282 | 0.367 |
| 1000 | 0.593 | 0.422 | 0.408 | 0.653 |

## 重点问题

- **Q1_univariate_low_df** SRSD 池无单变量非单项式任务(1 变量物理式几乎全是单项式);该问题由合成台架 QUAD_LIN 回答(见合成报告 §2)。
- **Q2_averaging_blindness** 多变量 H1 任务中 exog-HC3 优于 form_check 最多的前 3:feynman-ii.11.27(fc=0.04 vs exog=1.00); feynman-ii.11.28(fc=0.05 vs exog=1.00); feynman-iii.4.32(fc=0.05 vs exog=0.99)
- **Q3_hc3_hetero_fpr** {'50': {'classical': 0.04988461538461538, 'hc3': 0.42196153846153844, 'exog_hc3': 0.5637307692307693}, '200': {'classical': 0.049692307692307695, 'hc3': 0.41673076923076924, 'exog_hc3': 0.5632307692307692}, '1000': {'classical': 0.05269230769230769, 'hc3': 0.4178846153846154, 'exog_hc3': 0.5646153846153846}}
- **Q4_alarm_extrapolation** {'auc': {'fc': 0.2857667969205327, 'fitted_f': 0.25942463070068456, 'fitted_hc3': 0.6078050588594601, 'exog_hc3': 0.5006655751872304}, 'alarm_pr': {'form_check': {'precision': 0.3912631703563395, 'recall': 0.6286684782608696}, 'RESET fitted F': {'precision': 0.5638041608008759, 'recall': 0.4897282608695652}, 'RESET fitted HC3': {'precision': 0.1710393788564779, 'recall': 0.47529891304347827}, 'RESET exog HC3': {'precision': 0.24415435174090658, 'recall': 0.9468478260869565}}}

## 报警→外推误差关联

按"远端外推相对误差>10%"为正类的 p 值 AUC:fc=0.286, fitted_f=0.259, fitted_hc3=0.608, exog_hc3=0.501


## 图(SRSD)

- `figures/srsd_calibration.png` — H0 p 值校准(真实单项式)
- `figures/srsd_per_task_power.png` — 逐任务功效
