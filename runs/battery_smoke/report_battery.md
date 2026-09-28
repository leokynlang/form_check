# 形式合格性电池 — battery_smoke

行数 2280,数据层失败 0(显式登记)。
wild bootstrap 零分布(Rademacher);α=0.05 每族独立;viable=未拒绝族。

## 层1 合成校准(power 族尺寸,结构 H0+噪声层)

| 场景 | 尺寸 |
|---|---|
| H0_BASELINE | 0.067 |
| H0_CORR | 0.050 |
| H0_MULTIVAR | 0.075 |
| H1_AUTOCORR | 0.067 |
| H1_HEAVYTAIL | 0.058 |
| H1_OUTLIER | 0.067 |
| H1_OUTLIER_GRP | 0.058 |
| H1_VAR_GROWING | 0.050 |

池化:**0.061**(目标 0.05)

## 层2 合成判别(EXP/LOG/RAT)

### COLLINEAR_CURV(truth=none)
- rel=0.005: containment=nan 假升级=0.033 假包含={}
- rel=0.01: containment=nan 假升级=0.000 假包含={}
- rel=0.03: containment=nan 假升级=0.033 假包含={}
- rel=0.05: containment=nan 假升级=0.033 假包含={}
### CUBIC_LIN(truth=none)
- rel=0.005: containment=nan 假升级=0.167 假包含={}
- rel=0.01: containment=nan 假升级=0.300 假包含={}
- rel=0.03: containment=nan 假升级=0.333 假包含={}
- rel=0.05: containment=nan 假升级=0.400 假包含={}
### DILUTION(truth=none)
- rel=0.005: containment=nan 假升级=0.900 假包含={}
- rel=0.01: containment=nan 假升级=0.933 假包含={}
- rel=0.03: containment=nan 假升级=0.967 假包含={}
- rel=0.05: containment=nan 假升级=0.967 假包含={}
### EXP_MODEL(truth=exponential)
- rel=0.005: containment=1.000 假升级=0.000 假包含={'(none)': 30}
- rel=0.01: containment=1.000 假升级=0.000 假包含={'(none)': 30}
- rel=0.03: containment=1.000 假升级=0.000 假包含={'(none)': 30}
- rel=0.05: containment=0.967 假升级=0.033 假包含={'(none)': 30}
### INTERACT(truth=none)
- rel=0.005: containment=nan 假升级=0.067 假包含={}
- rel=0.01: containment=nan 假升级=0.167 假包含={}
- rel=0.03: containment=nan 假升级=0.033 假包含={}
- rel=0.05: containment=nan 假升级=0.067 假包含={}
### LOG_MODEL(truth=logarithmic)
- rel=0.005: containment=0.967 假升级=0.033 假包含={'rational': 26, '(none)': 4}
- rel=0.01: containment=0.933 假升级=0.000 假包含={'rational': 27, '(none)': 3}
- rel=0.03: containment=1.000 假升级=0.000 假包含={'rational': 29, '(none)': 1}
- rel=0.05: containment=1.000 假升级=0.000 假包含={'rational': 29, '(none)': 1}
### OMITTED_M(truth=none)
- rel=0.005: containment=nan 假升级=1.000 假包含={}
- rel=0.01: containment=nan 假升级=1.000 假包含={}
- rel=0.03: containment=nan 假升级=1.000 假包含={}
- rel=0.05: containment=nan 假升级=1.000 假包含={}
### OPPOSING(truth=none)
- rel=0.005: containment=nan 假升级=0.000 假包含={}
- rel=0.01: containment=nan 假升级=0.000 假包含={}
- rel=0.03: containment=nan 假升级=0.000 假包含={}
- rel=0.05: containment=nan 假升级=0.000 假包含={}
### QUAD_LIN(truth=none)
- rel=0.005: containment=nan 假升级=1.000 假包含={}
- rel=0.01: containment=nan 假升级=1.000 假包含={}
- rel=0.03: containment=nan 假升级=1.000 假包含={}
- rel=0.05: containment=nan 假升级=1.000 假包含={}
### RATIONAL(truth=rational)
- rel=0.005: containment=0.933 假升级=0.067 假包含={'(none)': 30}
- rel=0.01: containment=0.900 假升级=0.100 假包含={'(none)': 30}
- rel=0.03: containment=1.000 假升级=0.000 假包含={'(none)': 30}
- rel=0.05: containment=0.900 假升级=0.100 假包含={'(none)': 30}
### THRESHOLD(truth=none)
- rel=0.005: containment=nan 假升级=1.000 假包含={}
- rel=0.01: containment=nan 假升级=1.000 假包含={}
- rel=0.03: containment=nan 假升级=1.000 假包含={}
- rel=0.05: containment=nan 假升级=1.000 假包含={}

## 层2b truth=none 场景(升级正确结局)

| 场景 | 正确升级率 | 假包含率 | top 假包含组合 |
|---|---|---|---|
| COLLINEAR_CURV | 0.025 | 0.975 | rational(75) |
| CUBIC_LIN | 0.300 | 0.700 | power(75) |
| DILUTION | 0.942 | 0.058 | rational(7) |
| INTERACT | 0.083 | 0.917 | power(110) |
| OMITTED_M | 1.000 | 0.000 | -(0) |
| OPPOSING | 0.000 | 1.000 | power+rational(96) |
| QUAD_LIN | 1.000 | 0.000 | -(0) |
| THRESHOLD | 1.000 | 0.000 | -(0) |