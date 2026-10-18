# 形式合格性电池 — battery_full

行数 94000,数据层失败 0(显式登记)。
wild bootstrap 零分布(Rademacher);α=0.05 每族独立;viable=未拒绝族。

## 层1 合成校准(power 族尺寸,结构 H0+噪声层)

| 场景 | 尺寸 |
|---|---|
| H0_BASELINE | 0.053 |
| H0_CORR | 0.049 |
| H0_MULTIVAR | 0.049 |
| H1_AUTOCORR | 0.049 |
| H1_HEAVYTAIL | 0.054 |
| H1_OUTLIER | 0.043 |
| H1_OUTLIER_GRP | 0.045 |
| H1_VAR_GROWING | 0.054 |

池化:**0.050**(目标 0.05)

## 层2 合成判别(EXP/LOG/RAT)

### COLLINEAR_CURV(truth=none)
- rel=0.005: containment=nan 假升级=0.036 假包含={}
- rel=0.01: containment=nan 假升级=0.034 假包含={}
- rel=0.03: containment=nan 假升级=0.042 假包含={}
- rel=0.05: containment=nan 假升级=0.024 假包含={}
### CUBIC_LIN(truth=none)
- rel=0.005: containment=nan 假升级=0.264 假包含={}
- rel=0.01: containment=nan 假升级=0.292 假包含={}
- rel=0.03: containment=nan 假升级=0.230 假包含={}
- rel=0.05: containment=nan 假升级=0.178 假包含={}
### DILUTION(truth=none)
- rel=0.005: containment=nan 假升级=0.924 假包含={}
- rel=0.01: containment=nan 假升级=0.918 假包含={}
- rel=0.03: containment=nan 假升级=0.924 假包含={}
- rel=0.05: containment=nan 假升级=0.910 假包含={}
### EXP_MODEL(truth=exponential)
- rel=0.005: containment=0.956 假升级=0.044 假包含={'(none)': 500}
- rel=0.01: containment=0.934 假升级=0.066 假包含={'(none)': 500}
- rel=0.03: containment=0.956 假升级=0.044 假包含={'(none)': 500}
- rel=0.05: containment=0.956 假升级=0.044 假包含={'(none)': 500}
### INTERACT(truth=none)
- rel=0.005: containment=nan 假升级=0.058 假包含={}
- rel=0.01: containment=nan 假升级=0.070 假包含={}
- rel=0.03: containment=nan 假升级=0.046 假包含={}
- rel=0.05: containment=nan 假升级=0.044 假包含={}
### LOG_MODEL(truth=logarithmic)
- rel=0.005: containment=0.928 假升级=0.004 假包含={'rational': 465, '(none)': 35}
- rel=0.01: containment=0.944 假升级=0.006 假包含={'rational': 454, '(none)': 46}
- rel=0.03: containment=0.934 假升级=0.008 假包含={'rational': 445, '(none)': 55}
- rel=0.05: containment=0.956 假升级=0.002 假包含={'rational': 465, '(none)': 35}
### OMITTED_M(truth=none)
- rel=0.005: containment=nan 假升级=1.000 假包含={}
- rel=0.01: containment=nan 假升级=1.000 假包含={}
- rel=0.03: containment=nan 假升级=1.000 假包含={}
- rel=0.05: containment=nan 假升级=1.000 假包含={}
### OPPOSING(truth=none)
- rel=0.005: containment=nan 假升级=0.002 假包含={}
- rel=0.01: containment=nan 假升级=0.000 假包含={}
- rel=0.03: containment=nan 假升级=0.002 假包含={}
- rel=0.05: containment=nan 假升级=0.000 假包含={}
### QUAD_LIN(truth=none)
- rel=0.005: containment=nan 假升级=1.000 假包含={}
- rel=0.01: containment=nan 假升级=1.000 假包含={}
- rel=0.03: containment=nan 假升级=1.000 假包含={}
- rel=0.05: containment=nan 假升级=1.000 假包含={}
### RATIONAL(truth=rational)
- rel=0.005: containment=0.938 假升级=0.062 假包含={'(none)': 500}
- rel=0.01: containment=0.918 假升级=0.082 假包含={'(none)': 500}
- rel=0.03: containment=0.948 假升级=0.052 假包含={'(none)': 500}
- rel=0.05: containment=0.932 假升级=0.068 假包含={'(none)': 500}
### THRESHOLD(truth=none)
- rel=0.005: containment=nan 假升级=1.000 假包含={}
- rel=0.01: containment=nan 假升级=1.000 假包含={}
- rel=0.03: containment=nan 假升级=1.000 假包含={}
- rel=0.05: containment=nan 假升级=1.000 假包含={}

## 层2b truth=none 场景(升级正确结局)

| 场景 | 正确升级率 | 假包含率 | top 假包含组合 |
|---|---|---|---|
| COLLINEAR_CURV | 0.034 | 0.966 | rational(1220) |
| CUBIC_LIN | 0.241 | 0.759 | power(1358) |
| DILUTION | 0.919 | 0.081 | rational(162) |
| INTERACT | 0.054 | 0.946 | power(1891) |
| OMITTED_M | 1.000 | 0.000 | -(0) |
| OPPOSING | 0.001 | 0.999 | power+rational(1681) |
| QUAD_LIN | 1.000 | 0.000 | -(0) |
| THRESHOLD | 1.000 | 0.000 | -(0) |

## 层3 SRSD(逐方程,H1 语义=必升级)

| 方程 | 升级率 | 中位|V| | top viable |
|---|---|---|---|
| feynman-bonus.15 | 0.029 | 2.0 | power,rational(0.677) |
| feynman-bonus.20 | 0.100 | 1.0 | exponential,rational(0.359) |
| feynman-i.10.7 | 0.096 | 1.0 | rational(0.558) |
| feynman-i.26.2 | 0.001 | 3.0 | power,exponential,rational(0.533) |
| feynman-i.30.5 | 0.169 | 1.0 | exponential(0.775) |
| feynman-i.37.4 | 0.025 | 1.0 | rational(0.49) |
| feynman-i.41.16 | 0.013 | 2.0 | power,rational(0.748) |
| feynman-i.48.2 | 0.233 | 1.0 | rational(0.463) |
| feynman-ii.11.28 | 0.001 | 4.0 | power,exponential,logarithmic,rational(0.87) |
| feynman-ii.13.23 | 0.240 | 1.0 | rational(0.46) |
| feynman-ii.13.34 | 0.086 | 1.0 | rational(0.615) |
| feynman-ii.38.14 | 0.669 | 0.0 | (0.669) |
| feynman-iii.4.32 | 0.003 | 2.0 | power,exponential(0.627) |
| feynman-iii.4.33 | 0.011 | 2.0 | power,rational(0.867) |