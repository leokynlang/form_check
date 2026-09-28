# 第三部分 守卫与运维层:异常、不确定度、多场景、监视

> 本部分覆盖 `anomaly.py`(463 行)、`uncertainty.py`(167 行)、`scenarios.py`(383 行)、`monitoring.py`(477 行)。共同点:它们都不发明新拟合,而是在拟合的**周围**回答四个工程问题——数据里的坏点怎么处理?仪器误差如何动摇结论?多个工况是同一机理吗?上线之后怎么持续监视?

---

## 3.1 `anomaly.py` — 后端异常守卫层

**文件定位**:套在"log 空间线性拟合"外面的守卫。只关心函数形式,不关心量纲(量纲合法性已由上游 π 空间保证)。三层检测 + 五种调整策略,全程可审计。

### 3.1.1 `check_data_health(data, target, features)`(anomaly.py:48-97)

- **【计算】** 数据层体检,五类问题各返回 `{'type', 'columns', 'severity'}`:missing_column / non_numeric / non_finite / non_positive / zero_variance,前三类与后两类均 fatal(非正值注明"log 空间拟合要求正值域";零方差注明"常数列不携带信息,且使杠杆值不可辨识")。
- **【原因】** fatal 即拒绝整个拟合——"先修数据再看异常"。零方差列 special-cased 是因为常数列会让 hat 杠杆值失去意义(任何点都同等极端)。
- **【导向】** `fit_with_guard` 第一步;health 清单进入返回 dict 的 `health` 键;有 fatal 时 success=False 直接返回。

### 3.1.2 `detect_residual_outliers(residuals, threshold=3.5)`(anomaly.py:104-121)

```python
med = float(np.median(r))
scale = 1.4826 * float(np.median(np.abs(r - med)))
z = (r - med) / scale
idx = [int(i) for i in np.where(np.abs(z) > threshold)[0]]
```

- **【计算】** MAD 稳健 z 分数:`z_i = (r_i − median)/(1.4826·MAD)`,|z|>3.5 判离群;MAD 退化(<1e-12,残差近乎常数)时显式 note、不判任何离群。
- **【原因】** 3σ 法的均值与标准差本身被离群点带偏(掩蔽效应);MAD 的 breakdown point 是 50%,离群点几乎不影响中位数与 MAD——"对最小二乘被离群点带偏这件事本身不敏感"。1.4826 是把 MAD 换算成正态一致尺度常数。
- **【导向】** `fit_with_guard` 的 `anomalies['residual_outliers']`;exclude 策略的剔除清单来源;audit 中每个被标记点的 `residual_z` 证据。

### 3.1.3 `detect_leverage(X, threshold=None)`(anomaly.py:128-145)

- **【计算】** hat 矩阵对角 `h_ii = x_i'(XᵀX)⁺x_i`(用 `np.einsum('ij,jk,ik->i')` 一次算完),默认阈值 `2p/n`(经验规则);n≤p+1 时杠杆不可辨识,显式 note。
- **【原因】** 高杠杆=该点在特征空间位置极端,会把拟合面拉向自己——与残差离群正交:杠杆点可以残差很小却支配整条拟合。语义区分(docstring/README):残差离群=疑似测量错误(可剔);高杠杆=条件极端的设计特征(默认只报告)。
- **【导向】** `anomalies['leverage']`;`exclude_leverage=False`(默认)时高杠杆点只进 detected/报告、不进剔除清单;`fit_with_guard` 返回的 audit 中 `flagged_sources` 标注 'residual'/'leverage'/'residual+leverage'。

### 3.1.4 `detect_collinearity(X, cond_threshold=30.0)`(anomaly.py:148-165)

- **【计算】** 特征列(不含常数列)逐列**标准化**(减均值除标准差)后求条件数 `np.linalg.cond(Z)`;>30 返回 warning 级问题(附条件数值)。
- **【原因】** 标准化消除单位/量级影响,只留下共线性本身;条件数>30 时指数估计不稳(Re/Pr 共享 μ 时常见)。CONTEXT.md 记录了相关的统计教训:"写死常数阈值(如 anomaly 条件数>30)同属脆弱类别"——这是已知保守近似,靠调用方在特定设计上复核。
- **【导向】** `anomalies['collinearity']`;`fit_with_guard(adjustment='auto')` 中"仅共线→自动选 ridge"的触发条件。

### 3.1.5 `robust_irls(X, y, k=1.345, max_iter=50, tol=1e-10)`(anomaly.py:172-200)

```python
w = np.where(np.abs(u) <= k, 1.0, k / np.maximum(np.abs(u), 1e-300))
sw = np.sqrt(w)[:, None]
new_coef, *_ = np.linalg.lstsq(X * sw, y * np.sqrt(w), rcond=None)
```

- **【计算】** Huber IRLS:初始 OLS → 每轮用残差的 MAD 尺度标准化 → |u|≤k 权重 1、其余 k/|u| → 加权 lstsq → 收敛判据 `max|Δcoef| < tol·(1+max|coef|)`。k=1.345 是 95% 渐近效率的标准取值。
- **【原因】** 降权而非剔除:"数据一个不丢"。与 exclude 的区别是语义性的:离群点仍参与拟合,只是话语权下降——适合"不确定是不是真错误"的灰区。
- **【导向】** `fit_with_guard(adjustment='robust')`;返回的 `weights` 让 after 候选携带 `n_downweighted/min_weight` 参数;audit['irls'] 记录迭代数/收敛/降权数。

### 3.1.6 `_ridge_coefficients` 与 `ridge_fit`(anomaly.py:207-264)

```python
alphas = np.logspace(-6, 3, 25)
for alpha in alphas:
    for valid in folds:            # 确定性 k-fold(默认 5)
        coef = _ridge_coefficients(X[train], y[train], alpha)
        squared_error += Σ(y_valid − X_valid @ coef)²
alpha = alphas[argmin(errors)]
```

- **【计算】** `_ridge_coefficients`:特征标准化后加岭罚 `β = (ZᵀZ + α·P)⁻¹Zᵀy`(P 对角阵、截距分量置 0 不罚),再把系数换回原尺度(斜率除尺度、截距补偿均值)。`ridge_fit`:25 点对数网格 α∈[1e-6,1e3],确定性 5-fold CV 累计平方误差选 α,全量重拟合;返回 `coef/alpha/cv_rmse_log/alpha_grid`。
- **【原因】** 共线性让 OLS 指数在正负间大幅摆动;岭罚以少量偏差换方差稳定。α 由数据(CV)选而非拍脑袋——与 form_check 的 bootstrap 阈值同一哲学:阈值是设计条件的函数,不写死。
- **【导向】** `fit_with_guard(adjustment='ridge')`;after 候选带 `ridge_alpha/cv_rmse_log`;audit['ridge'] 注明 'deterministic_kfold'。

### 3.1.7 `_logspace_candidate(...)`(anomaly.py:271-291)

- **【计算】** 把任意 log 空间系数包装成同构 `CandidateModel`(family 可为 power_law/power_law_robust/power_law_ridge/power_law_excluded);带权重时记录 `n_downweighted/min_weight`。
- **【原因】** before/after 必须同构才能逐字段比较。
- **【导向】** `fit_with_guard` 的 `before`/`after`;comparison 的逐特征指数位移计算。

### 3.1.8 `fit_with_guard(data, target, features, adjustment='report_only', ...)`(anomaly.py:294-463)

```python
if applied_adjustment == 'exclude' and not flagged:
    ... "检测到 N 个高杠杆点(条件极端, 非数据错误), 默认只报告不剔除"
if len(flagged) > max(1, n // 3):
    result['audit']['refused'] = "...拒绝自动剔除——这多半是模型形式错了, 不是数据错了。"
```

- **【计算】** 总装:①数据层 fatal → 拒绝;②log 变换构造设计矩阵 `[1, log X]`;③基准 OLS(共享引擎)作 `before`;④三层检测;⑤策略路由:auto 按检测到的异常类型选 robust(有残差离群)>ridge(仅共线)>report_only(无事不调);report_only 提前返回;⑥各策略执行调整并生成 `after`;⑦before/after 对比:`r2_before/r2_after/exponent_shift/max_exponent_shift`。两条保护性提前返回:检测到的只有高杠杆且策略是 exclude → 只报告不动手;被标记点超过样本 1/3 → **拒绝自动剔除**并在 audit['refused'] 说明原因与替代方案。
- **【原因】** 默认 report_only 对应"不悄悄改数据"哲学;1/3 上限的理由写在拒绝文本里:"这多半是模型形式错了,不是数据错了"——大规模剔除是形式失配的信号,该回到函数族选择层面(这正是 form_check 的职责,两层守卫在语义上衔接)。auto 的优先序:残差离群(点级问题)优先于共线性(设计问题)。
- **【导向】** ①独立 API(README 示例);②`discover_auto` 推荐幂律后的 `backend_adjustment`;③`EngineeringMonitor` 未直接使用它,但同一套检测思想在监视层以残差 z 分数复现;④validate_airfoil.py 的三路径之一。

---

## 3.2 `uncertainty.py` — 测量不确定度蒙特卡洛传播

**文件定位**:docstring 开宗明义:"`Variable.uncertainty` 自 v0.1 起声明,但从未接入计算。std_errors 只回答'回归残差噪声有多大'(条件于测量无误差),回答不了'仪器精度会如何动摇指数'。两者是独立的噪声来源,必须分开报告。"

### 3.2.1 `propagate_uncertainty(data, system, target, n_sims=200, random_state=42)`(uncertainty.py:38-167)

```python
for name, u in declared.items():
    perturbed[name] = cols[name] * np.exp(rng.normal(0.0, u, len(data)))
pi = compute_groups(perturbed, groups)
fit = fit_power_law(pi, y_col, x_cols)
exp_samples.append(list(fit['exponents'].values()))
...
out['exponents_mc'][name] = {'mean', 'std', 'q025', 'q975', 'summary'}
```

- **【计算】** 五步:①收集所有声明了 uncertainty 的变量(含 target);无任何声明 → 显式拒绝,reason 给出正确声明示例("框架不在无测量噪声信息时虚构传播结果");②数据门卫:`check_data` Level<2 → 拒绝(在 compute_groups 之前,绝不抛异常);③π 群定义只算一次(群由量纲决定,不随数值扰动变化);④MC 主循环 n_sims 次:每个不确定变量乘 `exp(N(0,u))`(log 空间加性高斯,与 log-log 拟合假设自洽且中位数无偏)→ 重投影 π → 重拟合幂律 → 收集指数与 log(prefactor) 样本;某次扰动把群打进负值/溢出 → 该次作废计数(`n_failed`);有效样本 < max(30, n_sims/10) → 失败并提示"不确定度声明可能过大";废弃率 >10% → 在 note 里警告"区间可能低估";⑤汇总:每指数 mean/std(ddof=1)/2.5%-97.5% 分位数/人读摘要;`regression_std_errors` 从无扰动基准拟合透传作对照。
- **【原因】** 方法论核心是**两种噪声的语义分离**(返回 dict 的 note 键内置):回归 std_errors 是拟合残差折算的系数标准误,n 越大越小;MC 传播的是测量噪声,反映仪器极限,不随样本量消失——"不应相加或直接比较"。乘性扰动用 exp(N(0,u)) 是因为声明的是**相对**误差。EIV(误差-in-变量)衰减是诚实的附带效应(docstring):扰动特征变量会让 MC 均值略低于无噪拟合值,这不是 bug 而是测量噪声对指数的总体影响。
- **【导向】** 独立 API(README 示例);换热器实跑结论记录在 CONTEXT.md(h±3% 下 Re 指数 0.805±0.0064,回归 se 0.0042——两个噪声源量级可互比但独立);测试 6 个(test_uncertainty,含"量纲独立变量零空间维 0"的玩具系统教训)。同 seed 逐位可复现(default_rng(42))。

---

## 3.3 `scenarios.py` — 多场景结构比较

**文件定位**:docstring:"同一套设备数据经常混合不同机器、材料、站点或运行模式。fit_scenarios 不预设'工况一定相同'或'每个工况必须完全独立',而是同时比较三种结构。"

### 3.3.1 `_validate(data, target, features, scenario)`(scenarios.py:27-57)

- **【计算】** 缺列/特征非空/场景列 NaN 检查;**标签碰撞检查**:不同场景标签字符串化后重名 → 拒绝(防 `machine_1` 与字符串 `"machine_1"` 混流);非数值/NaN/inf/非正值 → 拒绝;每场景样本数 ≥ len(features)+3,不足的场景列在 reason 里。
- **【原因】** shared_exponents 设计矩阵每场景需要一列指标变量,独立策略每场景需要独立拟合——样本下限保证折内可辨识。标签字符串化是 `predict_scenarios` 的键匹配语义,提前统一。
- **【导向】** `fit_scenarios` 的结构化失败。

### 3.3.2 `_fit_strategy(data, target, features, scenario, strategy)`(scenarios.py:71-179)

```python
elif strategy == "shared_exponents":
    indicators = np.column_stack([(scenario_values == key).astype(float) for key in label_keys])
    design = np.column_stack([indicators, log_x])   # 场景哑变量 + 共享 log 特征
```

- **【计算】** 三种结构一次设计矩阵实现:
  - `pooled`:`[1, log X]`,n_params = p+1;
  - `shared_exponents`:`[场景哑变量, log X]`——每场景独立截距(=独立倍率 exp(截距)),指数全场景共享;n_params = 场景数+p;
  - `independent`:逐场景子矩阵各自 `[1, log X]`,rank 累加、条件数取最大;n_params = 场景数×(p+1)。
  - 收尾统一:秩不足(n_params > rank)→ 结构化失败"无法辨识该多场景结构";R²/AIC/AICc/BIC/train_rmse_log/design_rank/condition_number。
- **【原因】** 三种结构是三种机理假设的显式编码:同一机理不同标定(shared_exponents,工程最常见)vs 完全同质(pooled)vs 完全异质(independent)。用哑变量截距实现"共享指数+独立倍率"使三者都能走同一 OLS 引擎与同一 AICc 公式,可比性免费获得。
- **【导向】** `_assess_strategy` 的模型;`fit_scenarios` 的 final_model(params 键按策略不同:prefactor+exponents / scenario_prefactors+exponents / scenario_models)。

### 3.3.3 `predict_scenarios(model, data)`(scenarios.py:182-227)

- **【计算】** 按策略预测;训练未见场景 → `ValueError("出现训练时未见场景 ...")`(显式拒绝,绝不静默套用旧模型);特征含非有限值/非正值同样拒绝。
- **【原因】** README:"训练时未见的场景不会被静默套用到旧模型上"——这是框架"不猜"哲学在预测端的体现。
- **【导向】** `_assess_strategy` 的折内预测;用户侧预测 API(README 示例);monitoring 对场景模型的预测。

### 3.3.4 `_stratified_folds(data, scenario, n_splits, random_state)`(scenarios.py:230-241)

- **【计算】** 每个场景组内部独立置换后 array_split,再把各场景的第 k 片拼接成第 k 折。
- **【原因】** 保证**每折都包含所有已知场景**——否则 independent 策略在缺场景的折上根本无法拟合,评估会系统性偏向 pooled。
- **【导向】** `_assess_strategy` 的折划分。

### 3.3.5 `_assess_strategy` 与 `fit_scenarios`(scenarios.py:244-384)

- **【计算】** `_assess_strategy`:拟合全量模型 → 场景分层 k-fold(逐折重拟合+预测,单折失败记入 failures 不中断)→ out-of-fold log 空间 `cv_r_squared/cv_rmse_log/cv_r2_standard_error/fold_r_squared`。`fit_scenarios`:三种策略全评 → best 取 CV 最高 → 一标准误差容差(max(0.005, best 的折间标准误))内取竞争集 → 竞争集内按 (n_params, AICc, BIC, 策略先验序 pooled<shared<independent) 选推荐;排序与 rank/recommended 标注;附 `scenario_counts` 与 `selection.rule` 说明文本。
- **【原因】** 与 autodiscovery 完全同构的选择逻辑(一标准误差规则偏好简单结构):"进入一标准误差范围后优先选择参数更少的结构,再比较 AICc/BIC,从而避免把随机噪声误判为场景差异"(README)。
- **【导向】** `fit_scenarios` 是独立 API + `EngineeringMonitor.fit_baseline`(有场景列时)的发现步骤;`predict_scenarios` 的模型来源;monitoring 的 `_scenario_probabilities` 从其模型键读取。

---

## 3.4 `monitoring.py` — 工程异常持续监视

**文件定位**:docstring:"把四个互补信号组合,而不是把每个异常行都当成坏测量":schema 硬失败 / 运行包络越界 / 特征与场景漂移 / 残差漂移;外加连续异常批次计数。

### 3.4.1 `DEFAULT_THRESHOLDS` 与 `_robust_scale`(monitoring.py:25-47)

- **【计算】** 13 个默认阈值(PSI 0.20/0.35、越界率 0.10/0.25、残差 z 3.5、离群率 0.05/0.15、中心偏移 z 2.5/4.0、RMSE 比 1.5/2.5、场景构成 TV 0.20/0.35、持续批次 3)。`_robust_scale = max(1.4826·MAD, std, 1e-8)`。
- **【原因】** 阈值全部可被 `thresholds={...}` 覆盖且未知键报错(防拼错静默生效);尺度取 MAD 与标准差的较大者——双保险,残差分布异常宽时退化为普通尺度,极窄时下限 1e-8 防除零。
- **【导向】** `EngineeringMonitor` 全部判级的分母;`snapshot()` 透出 thresholds 供审计。

### 3.4.2 `population_stability_index(reference, current)` 与 `_psi_reference`(monitoring.py:50-69)

```python
quantiles = np.quantile(values, [0.2, 0.4, 0.6, 0.8])
edges = np.concatenate(([-np.inf], internal, [np.inf]))
...
return float(np.sum((actual - expected) * np.log(actual / expected)))
```

- **【计算】** 基线分位数(20/40/60/80%)做 5 桶边界,基线概率存档;PSI = Σ (actual−expected)·ln(actual/expected),两侧概率截断在 1e-6 防除零/ln(0)。
- **【原因】** PSI 是风控行业标准的分布漂移量:对桶数不敏感、对称、对大块偏移敏感。用基线自身分位数分桶保证基线均匀 5 桶——漂移量对桶边界的依赖降到最低。
- **【导向】** `inspect` 的逐特征 `metrics['features'][f]['psi']`;`feature_distribution_drift` 告警;测试(test_monitoring)构造已知漂移验证 PSI 数值行为。

### 3.4.3 `EngineeringMonitor.__init__`(monitoring.py:95-172)

- **【计算】** 构造期校验与基线标定:①模型类型判别(多场景模型 dict 带 strategy 键 vs CandidateModel);失败模型、模型/监视器的场景列-特征-目标不一致全部拒绝;②阈值合并;③基线预测必须有限且为正;④残差标定:`residual_center = median(log y − log ŷ)`、`residual_scale = _robust_scale(residual)`、`baseline_rmse_log`;⑤逐特征运行包络:log 空间 1%/99% 分位数外扩 0.5×IQR 作上下界 + 存 PSI 参考分桶;⑥场景构成基线概率。
- **【原因】** 监视器的一切判级都相对基线:残差中心/尺度、包络、PSI 参考都是"正常状态"的操作定义。IQR 外扩而不是 min/max 作包络,避免基线极端值把包络撑得永远不报警。
- **【导向】** 全部内部状态(`_residual_center/_residual_scale/_baseline_rmse_log/_feature_reference/_scenario_probabilities`)只在本类使用;`fit_baseline` 工厂把发现+标定合成一步。

### 3.4.4 `EngineeringMonitor.fit_baseline(...)`(monitoring.py:174-212)

- **【计算】** 无场景列 → `discover_auto`;有 → `fit_scenarios`;失败透传结构化原因;监视器构造异常折叠为 `{'success': False, 'reason': ...}` 不抛给调用方。
- **【原因】** "一键发现+建监视器"是 README 主打的入口;监视器构造失败也结构化,与全框架"绝不抛异常"的纪律一致。
- **【导向】** 用户侧 `built['monitor'].inspect(batch)` 工作流;discovery 结果随返回值透出供审计。

### 3.4.5 `EngineeringMonitor.inspect(batch, batch_id, require_target=True)`(monitoring.py:287-455)

```python
residual = np.log(batch[target]) - np.log(prediction)
z = (residual - self._residual_center) / self._residual_scale
...
if alerts: self.consecutive_abnormal_batches += 1
else: self.consecutive_abnormal_batches = 0
if alerts and self.consecutive_abnormal_batches >= persistent:
    alerts.append(self._alert("persistent_abnormal_batches", "critical", ...))
```

- **【计算】** 单批次五段检查(全程持 RLock):
  1. `_validate_frame`:缺列/空批/非数值/NaN/inf/非正值/场景 NaN → **失败报告**(critical,连续计数 +1);
  2. 未知场景 → 失败报告;已知场景 → 构成总变差 `TV = 0.5·Σ|p_now − p_base|` → `scenario_mix_drift` 告警(可选关闭);
  3. 逐特征:log 值超出包络比例 → `operating_envelope_excursion`;PSI → `feature_distribution_drift`;均按 warning/critical 双阈值;
  4. `require_target=True` 时:模型预测 → 残差 → 稳健 z 分数,行级异常位置 `row_anomaly_positions`;三个批次级指标:离群率→`residual_outliers`、残差中心偏移 z→`residual_center_drift`、RMSE 比→`model_error_degradation`;
  5. 告警合并:有任何告警则连续计数 +1,无告警清零;连续 ≥3 批追加 `persistent_abnormal_batches`(critical);status = critical/warning/normal;`_record` 写入有界历史(history_limit=100)。
- **【原因】** 四信号分离的动机(docstring):"坏数据(缺列/NaN)、超出适用范围(包络)、工况人口改变(漂移)、相关式失效(残差)"是四种不同的工程事件,处理方式不同,混在一起会产生"每个异常行都是坏测量"的误判。连续性升级的原因:"弱但持久的信号应当升级,而不应对单个噪声批次过度反应"。无目标预警(`require_target=False`)对应现场目标回传滞后的真实工况——先监视输入侧。
- **【导向】** 返回报告 dict(`status/alerts/metrics/row_anomaly_positions/consecutive_abnormal_batches`)是用户决策依据;监视器**只告警和留痕**,不删样本、不重训、不接管控制(README:"这些动作应由工程审批流程决定")——这是"仅服务工程认知,不接管生产端控制"哲学的最终落点。`snapshot()`/`reset_history()` 提供运行期状态管理与"清历史不清标定"的区分。

---

## 3.5 本部分小结

| 模块 | 守护的对象 | 核心统计量 | 输出的决策语义 |
|---|---|---|---|
| anomaly.py | 单次拟合 | MAD z、hat h_ii、标准化条件数、Huber 权、CV 岭 α | 检测→(可选)调整→审计 |
| uncertainty.py | 指数的可信度 | MC 均值/std/95% CI vs 回归 se | 两种噪声分开报告 |
| scenarios.py | 机理同质性 | 场景分层 CV R²、AICc/BIC | 三结构择一,宁简勿繁 |
| monitoring.py | 上线后的模型 | PSI、包络越界率、残差 z/RMSE 比、连续批次数 | 只告警留痕,不接管 |
