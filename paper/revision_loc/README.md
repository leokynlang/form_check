# revision_loc — HC3 假阳性故障定位(2026-09-07)

回应审稿意见第 4 条的定位实验。**结论:v4 论文 Q3 的"pinv 实现路径"表述仍不够准——
失效载体是 `pinv(X'X)` 的默认截断容差作用在条件数平方的 Gram 矩阵上;不是 pinv 与
SVD 两个实现族的本质差异(np.linalg.pinv 本身就是 SVD,gelsd)。**

## 证据文件
- `locate_pinv_failure.py` — E1 复现 / E2 载体分解 / E3 条件数扫描
- `results_localization.json` — 全部原始数字

## E1 复现(n=200, rel=1%, 60 seeds,同一数据喂四条路径)
| 任务 | harness HC3 | classical F | statsmodels | SVD-native |
|---|---|---|---|---|
| i.34.27 | **1.000** | .050 | .033 | .033 |
| iii.12.43 | **1.000** | .033 | .033 | .033 |
| i.12.1(对照) | .050 | .067 | .050 | .050 |

## E2 载体分解(i.34.27, 20 seeds,hat × 协方差 2×2)
| 路径 | 中位 p | 拒绝率 |
|---|---|---|
| gram_hat + **gram_cov** | 0 | **1.00** |
| svd_hat + **gram_cov** | 0 | **1.00** |
| gram_hat + svd_cov | .62 | .00 |
| svd_hat + svd_cov | .63 | .00 |
| statsmodels | .63 | .00 |

→ **失效在夹心协方差的 `pinv(X'X)` 因子,不在杠杆 hat。** 换成 SVD 的协方差立即恢复,
hat 用哪条路径无影响。cond(X)=10^10.7,cond(X'X)=10^21.05;Gram 截断最大偏差 1
(完全奇异),SVD 截断偏差 1.1e-06。

## E2b 三种候选解释裁决(i.34.27, 10 seeds)
| 修复方式 | 中位 p | 拒绝率 |
|---|---|---|
| pinv(G) 默认 rcond | 0 | 1.00 |
| pinv(G), rcond=0(不截断) | .241 | .30 |
| 列均衡化后 pinv(G) | 0 | 1.00 |
| SVD-on-X(statsmodels 同路) | .628 | .00 |

- **容差是直接死因**:seed 0 上同一矩阵,默认容差 p=0,rcond=0 p=0.433。
  机制:pinv 默认截断阈 = max(n,k)·eps·σ_max(G) ≈ 5.1e-3,而 σ_min(G) ≈ 8.5e-10
  ——被截方向比阈值低 7 个数量级,rank 4 被砍成 3。
- **尺度假设被排除**:列均衡化不恢复(σ 方向在 Gram 里已经死到底,缩放救不了)。
- **参考分布假设被排除**:各路径同用 F(q, n−k),良态任务上各实现 ≥3 位一致。
- **rcond=0 只部分修复**(拒绝率 .30,余留膨胀):cond(X'X)≈10^21.8,在 float64
  里反解 Gram 必然丢 ~6 位精度——**正确修复是"在高条件数下根本不组 X'X"**,
  对 X 直接做 SVD/QR(statsmodels 路线),或列均衡化后对 X 做 QR。
  均衡化本身有效的前提是不经过 Gram(与 Bellamy 2021 结论一致)。

## E3 条件数扫描(13 H0 任务 × 10 seeds)
拒绝率 1.00 的 5 个任务 cond(X) ∈ [10^7.6, 10^11.5];≤0.10 的任务除 ii.27.18
(10^7.6)外全部 < 10^7。点二列 r(任务中位 log10 cond, 拒绝率) = **0.844**
(论文 Q3 报 .78,聚合层级不同,方向一致)。

## 对论文的修正建议(第 4 条的准确表述)
原文"silent numerical pathology of pinv-squared-conditioning inference paths"
应改为类似的:**"在近共线增广设计上先组 Gram 矩阵 X'X 再求逆,条件数平方后
σ_min(X'X) 跌破 float64 噪声地板;`pinv` 的默认相对截断容差将该方向判为精确零,
Wald 统计量随之爆表。任何实现(pinv 或显式求逆)在此路径上同样失效;对 X 本身
做 SVD/QR 分解则完全恢复。"**

## 给 harness 的最小修复(common.py)
`_ols` 第 91 行 `XtX_inv = np.linalg.pinv(X.T @ X)` 与 `reset_test` 第 234 行
`V = XtX_inv @ meat @ XtX_inv` 是唯一需要改的两处:
1. hat 改由 SVD 因子取:`U, s, Vt = svd(X); hat = einsum("ij,ij->i", U, U)`;
   `XtX_inv = (Vt.T / s**2) @ Vt`(几何上等于 pinv(X)·pinv(X)^T,不组 Gram);
2. 新增登记类型 `gram_conditioning_hazard`:当 cond(X_aug) > 1e8 时拒绝走
   Gram 路径并登记(或将 C/D 的 wrapper 与 statsmodels 交叉验证,复用 A 已有的
   statistic_mismatch 哨兵机制——该哨兵目前只覆盖 A,不覆盖 C/D,这正是本次
   伪影"静默"通过的原因)。

## 复跑
```
"$LOCALAPPDATA/Programs/Python/Python313/python.exe" \
  "D:/form_check/paper/revision_loc/locate_pinv_failure.py"
```
依赖:3.13 环境(numpy 2.5.2 / statsmodels 0.15.0 / sympy / yaml);
SRSD 上游浅克隆在 %TEMP%/srsd-benchmark(已存在)。

## v5 全量重跑结果(runs/cd_svd/,rerun_cd_svd.py)
- S1(142,500 fits):C/D 旧 vs 新**逐格一致**(H0 尺寸、全部 H1 功效差 ≤.02)。
  Gram 伪影不攻击合成宽域(cond ~10¹),Tables 1/4/7 不变——伪影边界被精确划定。
- S2(240,000 fits,12/13 任务哨兵逐位过):H0 尺寸 C .417→.053、D .563→.056;
  H1 raw 功效(n=200) C .623→.227、D .932→.661。v2/v3 的 headline 表是伪影。
- S3(11,760):H0 池 C .453→.046、D .450→.044(1/3 幂律 D=1.00 亦伪影)。
- S6(120 格)与工程台架 E(147 格,wb/smoke 缓存重放):A/B 零漂移;
  E 的 C 拒绝率 .735→.633、D .639→.592。
- i.30.5 一个任务(7,500 行)数据流漂移、哨兵真失配,已剔除待披露
  (上游克隆与环境均同一,漂移限于该 3 变量 exp 方程内部;其 C/D 差异
  与 S1 EXP_MODEL twin 一致为零,不改变任何结论)。
- 汇总:runs/cd_svd/v5_master_summary.json;逐行数据 s{1,2,3}_cd_svd.csv。
  校准功效重推导:runs/cd_svd/v5_table3_calibrated.json(D .494/.675/.733,
  v4 的 C/D 校准行因阈值取自伪影 p 分布而作废);任务级 bootstrap:
  task_bootstrap_v5.py 输出(13 方程全部落 .052–.062,v4 的双峰消失);
  S3 家族表:s3_family_table_v5.py 输出;E 台架指纹:v5_e_bench_findings.json
  (31 格 p<1e-20 下溢,Gram 独有,SVD 全部回到正常范围)。

## 代码调试收尾(2026-09-07)
- real_data/smoke_test.py:输出路径参数化(--out,默认行为不变)。
  事故根因修复:重放曾硬编码回写 real_data/smoke_results.csv 覆盖已审计
  的 Gram 原版;原版已从缓存确定性恢复并另存
  runs/cd_svd/e_smoke_results_GRAM_original.csv。
- replay_s6_smoke_svd.py:replay_smoke 改走 --out 直写独立目录,不再回写
  原档案;sentinel_check 清除死分支(两表统一 utf-8-sig 读取)。
- rerun_cd_svd.py:删除 S6/E 死桩函数(重放职责归 replay 脚本,注释说明);
  去掉未用的 json 导入;S1 断点续写 partial 逻辑保留。
- 回归验证:8 脚本 py_compile 通过;replay 重跑后原档案逐字节未动、
  SVD 输出与修复前 150 行逐位一致、patch 对 Gram 路径的差异仍生效;
  三个套件行数(142500/240000/11760)与关键 JSON 值抽检通过。
