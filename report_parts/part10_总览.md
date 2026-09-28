# 第十部分 总览:代码地图、流水线、统计方法索引与工程纪律

---

## 10.1 全仓代码地图

```
dimreflex/ (仓库根 = 包本身,扁平布局)
├─ 核心包(13 文件,≈3,900 行,v0.5.0)
│   core.py        量纲引擎:pint 单位 → 量纲矩阵 → 精确有理零空间 → π 群
│   _ols.py        唯一 OLS 引擎:全框架 log 空间拟合的统计同源点
│   discover.py    发现层总装:_prepare_pi 公共前半段 + form_check + 三入口
│   autodiscovery.py  决策层:k-fold CV + 一标准误差 + AICc/BIC + 后端调整
│   families.py    函数族体系:CandidateModel schema + 模板方法基类 + 4 族 + 注册表
│   templates.py   工程知识注入:用户声明的 Group + fit_with_templates
│   anomaly.py     守卫层:三层检测 + 五种调整策略 + 全程审计
│   uncertainty.py 测量不确定度蒙特卡洛传播(与回归 se 分开报告)
│   scenarios.py   多场景结构比较(pooled/shared_exponents/independent)
│   monitoring.py  上线后监视:包络/PSI/残差/场景构成 + 连续升级
│   cli.py + launcher.py  命令行门面 + PyInstaller 薄启动器
│   __init__.py    31 个符号的公共 API 面
├─ validation/     正/负对照验证脚本(预注册判据 + 退出码承载结论)
├─ experiments/
│   form_check_vs_reset/  主基准 P1-P3(合成+SRSD+函数族分层+电池),15 文件
│   nist_strd/            NIST StRD 认证验证(数值正确性 + 二分归因)
│   fc_multidirection/    多方向扩展 + HX 混合检验 + 1000 组泛化
├─ tests/          116 个单元/集成测试(每类 docstring 携带事故档案)
├─ artifacts/ + */paper/ + build_hx_paper.py  三条论文产出线
├─ pyproject.toml  扁平打包 + 依赖钉界 + ruff 14 规则族全开
└─ heat_exchanger_data.csv / pump_bep_coefficients.csv / airfoil_self_noise.csv
```

## 10.2 端到端数据流水线

```
CSV + schema.json(单位/角色/不确定度)
  → System.check_data ────────────── Level 0-2 门控(缺列/样本量)
  → dimensionless_groups ─────────── 精确有理零空间 → 目标隔离 → 化简
  → compute_groups ───────────────── π 数值表(零点合法,NaN 是事故)
  → 逐特征降级 / 共享变量检测 ────── dropped_features / shared_vars 留痕
  → form_check ───────────────────── 偏残差曲率 + 本数据 bootstrap 零分布
  → 三入口分流
       power_law → fit_power_law(共享 _ols)
       multi     → discover_multi(R² 排序,失败候选带 reason)
       auto      → k-fold CV → 一标准误差竞争集 → AICc/BIC/参数量推荐
                     → 幂律推荐 + adjust_backend → fit_with_guard(auto)
  → _promote_level_3 ────────────── CV R²>0.9 才升级,矛盾写 warnings
  →(可选)propagate_uncertainty ── 测量噪声 → 指数 MC 区间(vs 回归 se)
  →(可选)fit_scenarios ────────── 场景分层 CV 三结构择一
  →(部署)EngineeringMonitor ───── 基线标定 → 批次五段检查 → 告警留痕
```

## 10.3 统计方法索引(全仓统计技术 → 落点)

| 统计方法 | 实现位置 | 服务的判断 |
|---|---|---|
| 精确有理高斯消元/零空间 | core.py | 可读的整数 π 指数 |
| 目标隔离 + 整数坐标下降 | core.py | 物理可解释的群基 |
| OLS + 经典标准误 | _ols.py | 幂律指数与置信度(全框架同源) |
| MAD 稳健 z 分数 | anomaly.py | 残差离群(50% breakdown point) |
| hat 杠杆 / 标准化条件数 | anomaly.py | 设计极端 / 共线性 |
| Huber IRLS | anomaly.py | 降权不剔除的稳健重拟合 |
| k-fold 交叉验证岭 | anomaly.py | 共线性下指数稳定化 |
| 残差 bootstrap 零分布 | discover.form_check | 形式失配的 per-data 阈值 |
| 残差置换零分布(pinv 末行向量化) | 实验层 common.py/multidir/hybrid | 检验的有限样本校准 |
| wild bootstrap(Rademacher) | battery.py | 异方差坐标下的合法零分布 |
| Ramsey RESET(fitted/exog × F/HC3) | 实验层 reset_test | 统计学对照系 |
| maxT(max|t|)+ 联合 F + Bonferroni 合并 | multidir/hybrid | 集中/分散失配双分支 |
| 贪心枢轴选择 + 修正 Gram-Schmidt | multidir/hybrid | 方向冗余消除与归因稳定 |
| k-fold CV R² + 一标准误差规则 | autodiscovery/scenarios | 复杂度惩罚的模型选择 |
| AIC/AICc/BIC(log 空间残差) | autodiscovery/scenarios | 复杂度诊断 |
| 变点穷举 + ΔAIC>10 双门槛 | families.Piecewise | 分段幂律不虚报 |
| 交替最小二乘(分块精确 LS) | families.RationalLog | 多特征有理拟合 |
| 蒙特卡洛误差传播 exp(N(0,u)) | uncertainty.py | 测量噪声 → 指数区间 |
| PSI(人口稳定指数) | monitoring.py | 特征分布漂移 |
| Wilson 分数区间 | 全实验层 + 构建器 | 比例的 95% CI |
| Mann-Whitney rank AUC | srsd_analyze.py | 报警→外推误差关联 |
| McNemar 配对检验 | HX 论文 | 配对功效差异 |
| brentq 隐式求根 | validate_friction | PK 真值数据生成 |
| Fortran 记数解析 + 双通道断言 | parse_nist.py | 认证数据 provenance |
| 四求解路径互检 + 三级缩放翻转 | run_nist.py | 数值伪影 vs 检测现象分离 |

## 10.4 全仓工程纪律(逐条可指认到代码)

1. **不猜**:单位必须 pint 可解析、角色必须白名单、不确定度必须显式声明(否则 propagate_uncertainty 拒绝)、场景未见拒绝预测(schema 期/数据期/预测期三道闸)。
2. **显式降级,绝不硬输出**:失败也是结构化结果(success=False + reason + 分类 fail 字段);CLI 退出码承载;Level<2 直接降级。
3. **不悄悄改数据**:anomaly 默认 report_only;剔除超 1/3 拒绝;dropped_features/excluded/flagged_sources 全留痕;monitor 只告警不接管。
4. **统计同源**:唯一 _ols 实现;TestStdErrorsParity 锁 parity;实验脚本自带迷你 OLS 以保持被试只读。
5. **预注册**:H0/H1 分类写死于代码;判据先于实验;预期翻转本身留痕(validate_friction C2、泵判据修正)。
6. **公平比较**:同一 p 尺度、同一零分布、匹配经验 FPR 下比功效、逐格不总平均、配对分析。
7. **失败可见**:fail 字段逐方法逐行;N/A 单元注册;论文把失败案例写进正文(F10/L1)。
8. **可复现**:seed 规则入配置(网格位置入种子);核心文件 sha256 入 env.json;断点续跑;同 seed 逐位一致有测试锁定。
9. **产品/实验隔离**:列名隔离(p_fc vs p_qda/p_md*/p_hx)+ docstring 声明"analogue 不得表述为产品成绩";实验目录外零写入。
10. **审计审计者**:statistic_mismatch 哨兵;raw/scaled 决策翻转检测;双通道数据裁决;C5b 事故存档。

## 10.5 一句话总评

dimreflex 是一个把"**诚实地知道自己不知道什么**"做成架构的量纲感知工程模型发现框架:π 定理给出可解释的搜索空间,函数族+CV 给出有证据的候选,form_check 及其后继(HX)给出"拟合得好"与"形式正确"的区分,守卫与监视层保证这条链路从数据进入持续到部署运行——而整个仓库的测试、验证脚本、基准实验与论文,本质上都在为同一件事积累审计记录:**每一条主张都能指认到产生它的代码、数据与种子**。
