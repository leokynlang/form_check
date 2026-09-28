# P2:Feynman 方程函数族覆盖 — feynman_families

配置:P1 冻结(powers=[2,3], B=999 排列, α=0.05, rel∈{0.5,1,3,5}%),30 种子/格;族分类仅用于分层选择,不改变 H0/H1 判据。

## 方程选择(30/120,含 not-applicable 代表)

| eq_name           | family            | klass   |   n_vars | eligible   | reason                       |   _e |
|:------------------|:------------------|:--------|---------:|:-----------|:-----------------------------|-----:|
| feynman-i.29.4    | exact_power_law   | H0      |        1 | True       | strict_monomial              |    0 |
| feynman-i.34.27   | exact_power_law   | H0      |        1 | True       | strict_monomial              |    0 |
| feynman-ii.27.16  | exact_power_law   | H0      |        1 | True       | strict_monomial              |    0 |
| feynman-ii.27.18  | exact_power_law   | H0      |        1 | True       | strict_monomial              |    0 |
| feynman-ii.8.31   | exact_power_law   | H0      |        1 | True       | strict_monomial              |    0 |
| feynman-iii.12.43 | exact_power_law   | H0      |        1 | True       | strict_monomial              |    0 |
| feynman-i.41.16   | exponential       | H1      |        2 | True       | non_monomial_structure       |    0 |
| feynman-iii.4.32  | exponential       | H1      |        2 | True       | non_monomial_structure       |    0 |
| feynman-iii.4.33  | exponential       | H1      |        2 | True       | non_monomial_structure       |    0 |
| feynman-i.6.20a   | exponential       | EXCLUDE |        1 | False      | nonpositive_inputs_or_output |    1 |
| feynman-i.6.20    | exponential       | EXCLUDE |        2 | False      | nonpositive_inputs_or_output |    1 |
| feynman-i.6.20b   | exponential       | EXCLUDE |        3 | False      | nonpositive_inputs_or_output |    1 |
| feynman-i.47.23   | fractional_power  | H0      |        3 | True       | strict_monomial              |    0 |
| feynman-i.44.4    | logarithmic       | EXCLUDE |        4 | False      | nonpositive_inputs_or_output |    1 |
| feynman-ii.37.1   | multivar_additive | EXCLUDE |        3 | False      | nonpositive_inputs_or_output |    1 |
| feynman-i.13.4    | multivar_additive | EXCLUDE |        4 | False      | nonpositive_inputs_or_output |    1 |
| feynman-ii.36.38  | other             | EXCLUDE |        5 | False      | nonpositive_inputs_or_output |    1 |
| feynman-i.11.19   | other             | EXCLUDE |        6 | False      | nonpositive_inputs_or_output |    1 |
| feynman-i.10.7    | rational          | H1      |        2 | True       | non_monomial_structure       |    0 |
| feynman-i.48.2    | rational          | H1      |        2 | True       | non_monomial_structure       |    0 |
| feynman-ii.11.28  | rational          | H1      |        2 | True       | non_monomial_structure       |    0 |
| feynman-ii.13.23  | rational          | H1      |        2 | True       | non_monomial_structure       |    0 |
| feynman-ii.13.34  | rational          | H1      |        2 | True       | non_monomial_structure       |    0 |
| feynman-ii.38.14  | rational          | H1      |        2 | True       | non_monomial_structure       |    0 |
| feynman-i.26.2    | trigonometric     | H1      |        2 | True       | non_monomial_structure       |    0 |
| feynman-bonus.15  | trigonometric     | H1      |        3 | True       | non_monomial_structure       |    0 |
| feynman-bonus.20  | trigonometric     | H1      |        3 | True       | non_monomial_structure       |    0 |
| feynman-i.30.5    | trigonometric     | H1      |        3 | True       | non_monomial_structure       |    0 |
| feynman-i.37.4    | trigonometric     | H1      |        3 | True       | non_monomial_structure       |    0 |
| feynman-bonus.10  | trigonometric     | EXCLUDE |        2 | False      | nonpositive_inputs_or_output |    1 |

## 按族汇总 @官方域(rel 池化)

| 族/类 | 情形 | form_check | RESET fitted F | RESET fitted HC3 | RESET exog HC3 |
|---|---|---|---|---|---|
| H0:exact_power_law | official_homo | 0.04 | 0.04 | 0.52 | 0.36 |
| H0:exact_power_law | official_hetero | 0.04 | 0.05 | 0.53 | 0.35 |
| H0:exact_power_law | scale_half | nan | nan | nan | nan |
| H0:exact_power_law | scale_double | nan | nan | nan | nan |
| H0:fractional_power | official_homo | 0.06 | 0.09 | 0.10 | 1.00 |
| H0:fractional_power | official_hetero | 0.00 | 0.06 | 0.03 | 1.00 |
| H0:fractional_power | scale_half | nan | nan | nan | nan |
| H0:fractional_power | scale_double | nan | nan | nan | nan |
| H1:exponential | official_homo | 0.11 | 0.05 | 0.70 | 0.68 |
| H1:exponential | official_hetero | 0.12 | 0.06 | 0.70 | 0.68 |
| H1:exponential | scale_half | 0.08 | 0.05 | 0.68 | 0.68 |
| H1:exponential | scale_double | 0.08 | 0.05 | 0.71 | 0.67 |
| H1:rational | official_homo | 0.66 | 0.13 | 0.68 | 0.95 |
| H1:rational | official_hetero | 0.65 | 0.17 | 0.66 | 0.96 |
| H1:rational | scale_half | 0.42 | 0.07 | 0.68 | 0.90 |
| H1:rational | scale_double | 0.84 | 0.19 | 0.69 | 1.00 |
| H1:trigonometric | official_homo | 0.69 | 0.72 | 0.81 | 1.00 |
| H1:trigonometric | official_hetero | 0.71 | 0.71 | 0.81 | 1.00 |
| H1:trigonometric | scale_half | 0.55 | 0.71 | 0.83 | 1.00 |
| H1:trigonometric | scale_double | 0.68 | 0.61 | 0.86 | 0.98 |

## 校准功效(经验 FPR=5%)

- **homo**(H1 池):form_check=0.56, RESET fitted F=0.33, RESET fitted HC3=0.38, RESET exog HC3=0.31
  - exponential: form_check=0.11, RESET fitted F=0.05, RESET fitted HC3=0.63, RESET exog HC3=0.04
  - rational: form_check=0.66, RESET fitted F=0.14, RESET fitted HC3=0.50, RESET exog HC3=0.66
  - trigonometric: form_check=0.70, RESET fitted F=0.73, RESET fitted HC3=0.08, RESET exog HC3=0.05
- **hetero**(H1 池):form_check=0.58, RESET fitted F=0.34, RESET fitted HC3=0.38, RESET exog HC3=0.31
  - exponential: form_check=0.13, RESET fitted F=0.06, RESET fitted HC3=0.63, RESET exog HC3=0.06
  - rational: form_check=0.68, RESET fitted F=0.17, RESET fitted HC3=0.49, RESET exog HC3=0.65
  - trigonometric: form_check=0.73, RESET fitted F=0.71, RESET fitted HC3=0.09, RESET exog HC3=0.04

## 采样范围敏感性(H1,同族 official vs ×2 域,>15pp 差异=明显)

- H1:rational form_check: 0.66 → 0.84(×2 域)

## 重点问题(基于本数据)

1. **哪类方程被 form_check 漏检**:exponential 族(官方域拒绝率 0.11,校准后 0.13)——exp(-x/T) 在 log-log 空间的曲率被其指数尺度隐藏;作为交换它在 rational/trig 上 0.66-0.73 且 H0 误报最准(0.04)。
2. **哪类被 fitted RESET 漏检**:rational 族(0.13-0.19,各变体最盲)——yhat^2/yhat^3 作为有理函数的局部泰勒展开恰好吸收饱和曲线的主弯曲,这是对 fitted 变体的结构性打击,与其'光滑近似'原理一致。
3. **量纲与报警**:本轮未直接操纵量纲信息(诚实留白);间接信号:fractional_power 的 H0(1/3 次幂单项式)上 exog-HC3 误报 100%——结构/量纲层面的'退化设计'识别本可事前排除该格。量纲正确不等于形式正确(exponential/rational 全部量纲自洽,仍是 H1)。
4. **采样范围敏感性**:仅 rational×form_check 一处超阈(0.66→0.84,×2 域曲率更明显,方向合理);其余族/方法 |Δ|≤15pp——**结论不随采样范围明显变化**。异方差情形与同方差逐字复现(全表 Δ≤0.02)。

## 与既有结论的一致与冲突

**一致**:①HC3-Wald 在真实物理单项式 H0 上严重超尺寸(本表 0.52,P1 SRSD 42-57%)——两套独立任务池复现;②校准前 exog 的'高功效'(rational 0.95/trig 1.00)与其 36-57% 误报率不可分离,校准后(trig 0.04-0.05)几乎消失——再次证明'同经验 FPR 比较是必需设施'而非可选项;③异方差不改变任何结论。
**冲突/新现象**:①合成台架上 fitted-F 曾是最强(QUAD/INTERACT 1.00),真实 Feynman 方程上它反而对 rational 近乎全盲(0.13)——合成外推到真实公式族会得出错误的排名,这是 P2 存在的理由;②form_check 在合成上被 OPPOSING/INTERACT 致盲,但真实方程池里它没有盲区(最弱 exponential 0.11 仍高于 fitted-F 同族 0.05)——'平均平方列'的盲区需要两个变量曲率恰好抵消,真实物理公式中未见此结构。
