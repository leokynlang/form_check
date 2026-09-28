# 第二部分 发现层:从 π 空间到候选模型

> 本部分覆盖 `discover.py`(518 行)、`autodiscovery.py`(322 行)、`families.py`(650 行)、`templates.py`(141 行)。这四个文件回答"在无量纲空间里如何找到候选函数关系、如何验证它、如何诚实地说出它的局限"。

---

## 2.1 `discover.py` — 发现层总装

**文件定位**:端到端流水线。`_prepare_pi` 是所有 π 入口的公共前半段;`form_check` 是 v0.5.0 新增的"形式合适性检查"(源于负对照暴露的缺口);三个入口 `discover(method=...)`、`discover_multi_pi`、`discover_auto_pi` 共享前半段。

### 2.1.1 `compute_groups(data, groups)`(discover.py:24-41)

```python
with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
    for i, group in enumerate(groups):
        col = np.ones(len(data))
        for var, exp in group.items():
            col = col * (data[var].values ** exp)
        out[f"PI_{i}_{name}"] = col
```

- **【计算】** 逐群计算无量纲数值 `∏ data[var]^exponent`,列名 `PI_{i}_{var^exp_...}`。
- **【原因】** π 群定义(指数字典)到可回归数值表的唯一步骤。零点的负指数会产生 inf——代码注释明确:这是**合法情形**(如攻角 0°),由 `_prepare_pi` 的逐特征降级分类处理,此处不告警、不拦截(`errstate` 显式抑制 numpy 噪声)。
- **【导向】** 产出 `pi_data` DataFrame;被 `_prepare_pi`(discover.py:135)与 `uncertainty.propagate_uncertainty`(每轮 MC 重投影,uncertainty.py:125)调用。列名同时充当 `y_col`/`x_cols` 的标识符,贯穿所有后续报告。

### 2.1.2 `fit_power_law(data, target, features)`(discover.py:44-101)

- **【计算】** log-log 线性拟合 `log y = log a + Σ bᵢ·log xᵢ`。两道前置防线:非有限值 → 显式失败;非正值 → 显式失败(discover.py:57-67)。设计矩阵 `[1, log X]` 送入共享 `_ols` 引擎;指数取 `exp(coef[0])` 为倍率、coef[1:] 为各特征指数;输出公式字符串、`std_errors`、`design_rank`、`identifiable`、`condition_number`、`n_samples`。
- **【原因】** 这是框架的基线假设:工程传热/流动世界大量满足幂律(Dittus-Boelter 类)。失败也是结构化结果(`success=False + reason`),绝不抛异常——调用方拿到的是"为什么不能拟合",不是 traceback。
- **【导向】** 返回 dict 被三方消费:① `discover(method='power_law')` 的 `result['fit']`;② `uncertainty.propagate_uncertainty` 的基准拟合与每轮 MC 拟合(uncertainty.py:108,133);③ `discover_auto` 经 `PowerLawFamily` 间接复用同一数值路径。`std_errors` 是泵相似律验证(validate_pump.py)判"指数落入 95% t 区间"的直接输入。

### 2.1.3 `_prepare_pi(data, system, target)`(discover.py:104-207)

**【计算】五步流水线**,每步都"留痕":

1. **数据就绪检查**(L117-126):`system.check_data`;Level<2 → note 非空,调用方显式降级。
2. **π 群生成**(L128-133):`dimensionless_groups(target=target)`;空群 → note。
3. **π 数值表与 Y/X 定位**(L135-150):`compute_groups`;恰好一个群含 target 作 `y_col`,其余作 `x_cols`;只有一群 → "无法进行回归"。
4. **逐特征降级**(L152-181):逐特征群检查——源数据含 NaN/inf 的群**不剔除**(留给拟合层整体失败);π 值含非正值/inf 的群剔除进 `dropped_features`,警告写入 `report['warnings']`;全部被剔 → note。
5. **共享变量检测**(L183-198):目标群变量集(去掉 target 本身)与每个特征群变量集求交,交集写入 `shared_vars` 并追加伪相关警告。
6. **形式检查钩子**(L200-205):`out['form_check'] = form_check(pi_data, y_col, x_cols)`,`_apply_form_check` 把报警写进 `report['warnings']`。

- **【原因】**(代码注释原文)"零点往往是合法物理量(如攻角 alpha=0 deg),让一个群毁掉全部函数族粒度过粗;改为剔除该群、留痕、用剩余群继续拟合。刻意的不对称:源数据含 NaN/inf 时不剔除……零是物理现实,NaN 是数据质量事故,后者框架保持最大保守。" 共享变量检测针对 π 群伪相关陷阱(如 Nu 与 Re 共享 D,R² 被数学性抬高)。形式检查钩子放在公共前半段是因为"首版只接 discover() 就漏了另外两个入口"——五入口共享防止绕行。
- **【导向】** 返回 dict(`report/groups/pi_data/y_col/x_cols/note/dropped_features/shared_vars/form_check`)是三个入口的共同输入;`dropped_features`/`shared_vars`/`form_check` 全部透传到最终结果字典供用户审计。

### 2.1.4 `form_check(data, target, features, n_boot=200, random_state=42)`(discover.py:253-340)

**这是本仓库方法论上最重要的单个函数**,docstring 完整给出了动机-方法-阈值三段论证:

```python
if features and len(features) == 1:
    quad = (lx[:, 0] - lx[:, 0].mean()) ** 2
else:
    quad = ((lx - lx.mean(axis=0)) ** 2).mean(axis=1)
...
c2_obs = _curv(ly, design)                      # 偏残差曲率(二次项系数)
for i in range(n_boot):
    boot_y = pred + rng.permutation(resid)      # 残差重排伪数据集
    null_c2[i] = _curv(boot_y, design)
t_val = (c2_obs - null_c2.mean()) / null_c2.std()
if abs(t_val) > 3.0: alarm
```

- **【计算】** log-log 空间中幂律=直线。步骤:①log 变换后构造设计矩阵 `[1, log x, 二次列]`,二次列以均值为中心的 `(log x − x̄)²`(多特征取各特征二次列的行均值);②`c2_obs` = 完整设计的 lstsq 最后一个系数,即"线性拟合之后残余的弯曲量"(偏残差曲率);③残差 bootstrap 零分布:在线性拟合预测值 `pred` 上加**打乱的残差**重建 n_boot=200 个伪数据集——继承本数据的 x 网格/样本量/噪声水平,唯一丢掉真实曲率;④`t = (c2_obs − mean_null)/std_null`,`|t|>3` 报警。
- **【原因】** 动机(docstring [动机]):框架原有检查全部关于"数据够不够"(L0-L3),没有"幂律形式合不合适"的。负对照(validate_friction.py,Prandtl-Kármán 对数律)实锤:幂律路径对非幂律真值静默输出 R²=0.99/Level 3,外推误差 -20% 而标准误极小——"统计精确 ≠ 形式正确"。阈值用残差 bootstrap 的原因([阈值为什么用残差 bootstrap]):曲率统计量的噪声尺度是数据设计的函数(σ ∝ x 量程跨度²,实测 5.9 倍),写死常数阈值必然某些设计误报、另一些漏报。代码中还有一段实现教训注释:二次列**必须平方**——首版漏写 `**2`,设计矩阵秩亏、lstsq 最小范数解给出任意系数,"零分布精确骑在观测值上,t 恒为 0——检查永远沉默",调试指纹是 `null.mean ≈ c2_obs`。
- **【导向】** 返回 `{checked, alarm, curvature, t_stat, n_samples, n_boot, reason}`:①挂到 `_prepare_pi` 结果的 `form_check` 键,透传到三个入口的最终结果;②报警时 `_apply_form_check` 写入 `data_report['warnings']`(措辞含"形式失配");③`discover(method='auto')` 中,推荐幂律族但报警时把矛盾**追加进** `recommended['selection_reason']`(discover.py:419-424)——用户不该自己去挖 form_check 才知道推荐有争议;④CLI 非量纲模式也调用它(cli.py:314-318);⑤它本身是 form_check_vs_reset 基准实验的被试对象(见第六部分),该实验把它与 Ramsey RESET 逐场景对比。

### 2.1.5 `_cv_r_squared` 与 `_promote_level_3`(discover.py:210-250)

```python
_LEVEL3_CV_R2 = 0.9
```

- **【计算】** `_cv_r_squared`:对指定函数族调 autodiscovery 的 `assess_family` 取 out-of-fold CV R²,失败/异常返回 None。`_promote_level_3`:CV R²>0.9 → `report['level']=3` 并写"Level 4 仍需跨工况外部验证"的警告;CV 不可用绝不升级;样本内 R²≥0.9 但 CV R²≤0.9 → 把"看似达标实则未升级(可能过拟合或形式失配)"写进 warnings。
- **【原因】** Level 3 的统一判据历史上有三条入口各说各话(样本内 R²),负对照证明样本内 R² 会被形式失配骗到 0.99。统一为确定性 k-fold CV R²>0.9 后,三个入口口径一致;矛盾警告专门捕捉"统计精确但形式错误"的情形。
- **【导向】** `report['level']` 升级到 3 会出现在所有 discover 入口的 `data_report` 里;CLI 打印"数据等级: Level N";测试(test_form_check)断言报警文本。矛盾警告是论文主张 1("fitted-but-wrong 现象学")的框架内对应物。

### 2.1.6 `discover(data, system, target, method='power_law')`(discover.py:351-434)

- **【计算】** 三条路径共用 `_prepare_pi` 前半段:
  - `power_law`:`fit_power_law` → 成功则算 CV 并尝试升 Level 3;
  - `multi`:`discover_multi` 全函数族 → 最优成功候选放 `result['fit']`;
  - `auto`:`_discover_auto` → 推荐幂律+form_check 报警时在 `selection_reason` 后追加矛盾说明。
  未知 method → 结构化失败。
- **【原因】** 单一入口承载三种认识论姿态:只给基线(幂律)、给全景(候选列表)、给自主判断(验证+推荐)。
- **【导向】** 返回 dict 是 README 快速入门、CLI 量纲模式、验证脚本(validate_pump/friction/airfoil)的总出口。

### 2.1.7 `discover_auto_pi` 与 `discover_multi_pi`(discover.py:437-518)

- **【计算】** 与 `discover` 的 auto/multi 路径相同的逻辑,但作为独立函数暴露完整参数(`families` 白名单、`n_splits`、`random_state`、`adjust_backend`)。
- **【原因】** README:"函数族漏斗与量纲感知引擎的合流点——π 空间中所有运算符天然量纲合法,量纲一致性免费获得"。
- **【导向】`discover_multi_pi` 返回 `candidates/best`;`discover_auto_pi` 返回 `autonomous/fit`;两者都含 `dropped_features/shared_vars/form_check` 全套审计键,供 `validation/` 脚本与测试直接消费。

---

## 2.2 `autodiscovery.py` — 验证与自主选择层

**文件定位**:docstring:"现有函数族层回答'什么能拟合这些数据',本模块添加决策层:确定性 k-fold 验证、复杂度诊断、一标准误差规则偏好简单模型、可选的后端调整。"

### 2.2.1 `predict_candidate(candidate, data, features)`(autodiscovery.py:17-82)

- **【计算】** 分家族把训练好的候选模型作用到新数据行:
  - `power_law*`:y = prefactor·∏ x^exponent;
  - `piecewise_power_law`:按 breakpoints 切分,每段支持新 schema(exponents 字典,多特征)与旧 schema(标量 exponent);
  - `rational_log`:L=log10(x)−log_center,log10 预测 = num(L)/den(L) + Σ b_j·(log10 x_j − center),再 10^ 回原空间;
  - 其他:调用 params 里的 `predict` 闭包;无接口 → `ValueError`(无法 CV 的家族显式声明)。
- **【原因】** k-fold 交叉验证需要在验证折上预测——每个家族必须有跨样本预测契约。统一入口让 CV 循环不关心家族内部结构。
- **【导向】** 被 `assess_family` 的每折预测(autodiscovery.py:166)与完整模型预测(L200)调用;也被 `EngineeringMonitor._predict` 用于批次残差计算(monitoring.py:248)。

### 2.2.2 `_validate_input`(autodiscovery.py:85-103)

- **【计算】** 检查缺列/特征非空/样本数≥4/非数值/NaN/inf/非正值,返回可读 reason。
- **【原因】** auto 层在 log 空间工作,正值域与有限值是硬前提;显式拒绝而非让下游崩溃。
- **【导向】** `discover_auto` 的第一道结构化失败。

### 2.2.3 `_folds(n_samples, n_splits, random_state)`(autodiscovery.py:106-110)

- **【计算】** 确定性 k-fold:`rng = default_rng(random_state)` 全量置换后 `array_split`。
- **【原因】** "确定性"是设计承诺:同 seed 同数据 → 同折 → 同 CV 分数(可复现性)。折数夹在 [2, n]。
- **【导向】** `assess_family` 与 `ridge_fit`(anomaly.py)使用同一风格的确定性折叠。

### 2.2.4 `_information_criteria(y, prediction, n_params)`(autodiscovery.py:113-129)

```python
sse = max(float(np.sum(residual ** 2)), np.finfo(float).tiny)
aic = n * np.log(sse / n) + 2 * k
correction = (2*k*(k+1)/(n-k-1) if n > k+1 else float("inf"))
aicc = aic + correction
bic = n * np.log(sse / n) + k * np.log(n)
```

- **【计算】** 在 **log 空间残差**上计算 AIC/AICc/BIC:残差 = log y − log ŷ。小样本校正项在 n≤k+1 时为 inf(AICc 不可信即拒绝比较)。
- **【原因】** CV R² 衡量预测力,IC 衡量复杂度代价——两者结合才能执行"一标准误差内偏好简单模型"。log 空间口径与拟合口径一致(幂律 = log 空间线性,残差正态假设在 log 空间更自洽)。
- **【导向】** `assessment` 的 `aicc/bic/train_rmse_log` 键;`discover_auto` 的推荐排序键(autodiscovery.py:265-274);`fit_scenarios` 的独立实现用同一公式(scenarios.py:158-176)。

### 2.2.5 `assess_family(data, target, features, family, n_splits=5, random_state=42)`(autodiscovery.py:132-218)

```python
for fold_index, valid in enumerate(folds):
    fold_model = FAMILIES[family]().fit(data.iloc[train], ...)
    predicted = predict_candidate(fold_model, data.iloc[valid], features)
    ...
    fold_r2.append(1 - ss_res / ss_tot)
...
cv_se = float(np.std(fold_r2, ddof=1) / np.sqrt(len(fold_r2)))
```

- **【计算】** 完整的单族评估:①全量拟合;②逐折重拟合+out-of-fold 预测(每折失败记录 `fold_failures` 而不中断);③汇总 out-of-fold log 空间 `cv_r_squared`、`cv_rmse_log`、折间标准误 `cv_r2_standard_error`、逐折 `fold_r_squared`;④完整模型预测算 train RMSE + AIC/AICc/BIC。任何一折失败或预测含非有限值/非正值 → 整族判"交叉验证未完成"(不给出半吊子分数)。
- **【原因】** 这是"自主发现不再用训练集 R² 决定最佳模型"的核心机制。折级失败留痕保证"验证未完成"与"验证通过但分数低"在语义上严格分开。
- **【导向】** ①被 `discover_auto` 对每个注册族调用;②被 `discover._cv_r_squared` 复用(三入口统一的 Level 3 判据来源)。

### 2.2.6 `discover_auto(data, target, features, families=None, n_splits=5, random_state=42, adjust_backend=True)`(autodiscovery.py:221-323)

```python
competitive = [item for item in successful
               if item["validation"]["cv_r_squared"] >= best_cv - tolerance]
recommended = min(competitive, key=lambda item: (
    round(item["aicc"], 10), round(item["bic"], 10),
    item["model"]["n_params"], family_priority.get(item["family"], 100),
    -item["validation"]["cv_r_squared"]))
```

- **【计算】** 全族评估 → 取 CV R² 最高者为 best → tolerance = max(0.005, best 的折间标准误) → 竞争集 = 落在 best 一标准误差内的族 → 在竞争集内按 (AICc, BIC, 参数量, 家族先验序 power_law<piecewise<rational, −CV R²) 字典序取最小者。排序后给每个 assessment 标 `rank`/`recommended`。推荐族是 power_law 且 `adjust_backend=True` → 调 `fit_with_guard(adjustment='auto')`,有 `after` 时 `final_model` 换成调整后候选;非幂律推荐 → backend 记 `not_applicable` 及原因。
- **【原因】** 一标准误差规则是经典的保守模型选择(Breiman/Acerentos):验证性能统计不可区分时,更简单的模型优先——这直接实现了"输出候选与证据,不造单一赢家"的哲学,但又要给用户一个可执行的推荐。family_priority 让工程可解释的幂律在同等验证力下自然胜出。后端调整只对幂律做,因为守卫层调整(降权/岭)定义在 log 线性设计上。
- **【导向】** 返回 dict(`success/assessments/recommended/final_model/backend_adjustment/selection`)是:①`discover(method='auto')` 的 `autonomous` 键;②`discover_auto_pi` 的主体;③`EngineeringMonitor.fit_baseline` 的发现步骤(monitoring.py:184-190)。`backend_adjustment` 键提供自动调整决策审计(applied_adjustment + decision 文本)。

---

## 2.3 `families.py` — 函数族体系

**文件定位**:统一候选 schema(`CandidateModel`)+ 模板方法基类(`FunctionFamily`)+ 四个内置族 + 注册表 + `discover_multi`。

### 2.3.1 `CandidateModel`(dataclass,families.py:35-60)

- **【计算】** 字段 `{family, formula, params, success, r_squared, n_params, n_samples, reason}`;`to_dict`/`from_dict`(宽松反序列化:忽略未知字段向前兼容)。
- **【原因】** 所有函数族产出同构候选,上层(排序、验证、报告)无需关心具体家族——这是注册表模式能工作的前提。失败候选的 `reason` 字段让"显式降级"有统一载体。
- **【导向】** 贯穿全框架:`discover_multi` 的列表元素、autodiscovery 的 assessment['model']、anomaly 的 before/after、monitoring 的模型持有形式。`from_dict` 被 `EngineeringMonitor._as_candidate` 用于持久化模型重建。

### 2.3.2 `FunctionFamily` 基类(families.py:67-99)

- **【计算】** 模板方法 `fit()`:调用子类 `_fit()`,捕获**一切**异常折叠为 `success=False` 的候选(带异常类型与消息);无 reason 的失败补默认原因。
- **【原因】** docstring:"这把'数据不足时显式降级,绝不硬输出'的框架哲学从调用方责任变成了基类保证:任何函数族实现都不可能向调用方抛数学异常。" 上层循环因此可以对所有族无差别遍历。
- **【导向】** 四个内置族与任何用户 `@register_family` 自定义族共享此契约;测试(test_families)会注入故意抛异常的族验证兜底。

### 2.3.3 `PowerLawFamily`(families.py:112-170)

- **【计算】** 与 `fit_power_law` 相同的 log-log OLS(共享引擎),额外把 `std_errors/design_rank/identifiable/condition_number` 放进 `params`。
- **【原因】** "现有 fit_power_law 的插件化迁移"——让基线幂律以统一候选形态参与多族竞争。
- **【导向】** `FAMILIES['power_law']`;autodiscovery 推荐链的默认优先族;守卫层调整的唯一作用对象。

### 2.3.4 `PiecewisePowerLawFamily`(families.py:183-372)

```python
m = max(int(self.min_segment_frac * n), p + 2)   # 每段最小样本
...
for c in cand:                                    # 变点网格 ≤200 点
    left = self._fit_segment(log_X, log_y, 0, c)
    right = self._fit_segment(log_X, log_y, c, n)
    aic = n * np.log(sse / n) + 2 * two_seg_k
...
d_aic = base_aic - best['aic']
d_r2 = (1 - best['sse'] / ss_tot) - base['r2']
if d_aic <= 10 or d_r2 < self.min_r2_gain: → 如实报告单段
bp = float(np.sqrt(x1[bp_cut - 1] * x1[bp_cut]))  # 间隙几何中点
```

- **【计算】** 穷举单变点搜索:①按第一特征排序,log 变换;②基准 = 单段完整多特征幂律,`base_aic = n·ln(SSE/n)+2(p+1)`;③变点网格(每段至少 m=max(0.2n, p+2) 个样本,大样本封顶 200 个候选切点)上两段各自 OLS,累计 SSE 最小者胜;④双门槛:`ΔAIC>10`(主门槛)且 `ΔR²≥0.01`(防退化下限),否则如实报单段并把增益写进 `params['note']`;⑤变点取"左段最后一点与右段第一点"间隙的**几何中点**,间隙本身存入 `params['breakpoint_interval']`。单特征时保留旧 schema 标量 `exponent` 键向后兼容。
- **【原因】** 工况转变(层流/湍流 Re≈2300)导致分段关联式。注释解释门槛取舍:"宽动态范围数据上 log 空间 R² 天然偏高,绝对 ΔR² 不能作为主门槛"——所以主门槛是带参数惩罚的 ΔAIC>10。诚实声明的取舍:变点仅沿第一特征(工况变量)搜索,其余特征以完整幂律进入每段——"这是设计取舍而非遗漏"(docstring)。报告区间而非假装精确:"数据只能把过渡约束在左右两点之间的间隙里"。
- **【导向】** `FAMILIES['piecewise_power_law']`;负对照 validate_friction 中 PK 对数律滑窗指数跨度证据的正面回应;predict_candidate 的分段预测分支。

### 2.3.5 `RationalLogFamily`(families.py:379-567)

```python
V = np.vander(L, N=self.num_degree, increasing=True)
def fit_rational(response):
    A = np.hstack([V, -W * response[:, None]])
    sol, _, _, _ = np.linalg.lstsq(A, response, rcond=None)
    den = 1.0 + W @ q
    if ... float(np.min(np.abs(den))) <= 0.05: return None   # 分母过小 → 发散
```

- **【计算】** log10 空间有理函数 `log10 y = num(L)/den(L) + Σ b_j·Z_j`:
  - 一次性线性化:`(1+W·q)·r̂ = V·c → [V, −r·W]·[c;q] = r`(把非线性有理拟合变成单次 lstsq);
  - 分母常数项固定为 1(可辨识性:有理函数分子分母同乘常数不变);
  - 多特征:分块交替最小二乘——(a)固定协变量系数 b,对 `log10y − Z·b` 拟合有理部分;(b)固定有理部分,对**中心化残差**做 Z 上 OLS 更新 b(残差均值交给 c0 吸收);每块是条件精确 LS,SSE 单调不增;收敛判据 `|ΔSSE| ≤ tol·(1+SSE)`,上限 50 轮;
  - 分母在数据上 `min|den| ≤ 0.05` → 显式失败;`alternating_ls_converged` 状态写入 params。
- **【原因】** Gnielinski 类关联式不是纯幂律,幂律族拟合它必然形式失配(负对照 C3:PK 台架上 rational_log CV 0.99989 击败 power_law 0.99036——负对照意外成为本族首个实战正对照)。docstring 诚实声明两条已知限制:交替 LS 收敛到局部平稳点而非全局最优;线性化在噪声大时有偏(errors-in-variables)。
- **【导向】** `FAMILIES['rational_log']`;predict_candidate 的 rational 分支;NIST 认证实验中 Filip(度-11 多项式真值)路径的天然候选。

### 2.3.6 `SymbolicRegressionFamily`(families.py:574-600)

- **【计算】** 只做一件事:尝试 `import pysr`;未安装 → 显式失败候选(说明安装方式);已安装 → 显式失败候选("占位尚未实现搜索配置")。
- **【原因】** docstring:"绝不假装成功。这符合框架'数据不足/能力不足时显式降级'的哲学。" 负对照 C4 要求 symbolic 路径保持显式失败。
- **【导向】** `FAMILIES['symbolic']`;所有 discover_multi 调用默认带出一个带原因的失败候选,提示用户扩展路径存在但未启用。

### 2.3.7 `FAMILIES` 注册表 + `register_family` + `discover_multi`(families.py:607-651)

- **【计算】** 注册表以类属性 `name` 为键(装饰器用法);四个内置族模块加载时注册。`discover_multi`:对每个(未注册的族也生成带 reason 的失败候选)族调用 `fit()`,成功且 R² 有限的按 R² 降序排前,失败候选排后供诊断。
- **【原因】** 注册表模式(项目评价第四节):第三方可以 `@dr.register_family` 无侵入扩展;排序只看 R² 是本层的决定——"复杂度惩罚在 autodiscovery 层实现,本层不改"(CONTEXT.md)。
- **【导向】** `discover(method='multi')`/`discover_multi_pi`/CLI 量纲模式/CLI 演示模式全部经此;`list(FAMILIES)` 是 discover_auto 的默认族清单。

---

## 2.4 `templates.py` — 工程知识显式注入

**文件定位**:docstring 点明动机:"Buckingham π 定理告诉我们零空间的维数,但不告诉我们应该选哪组基。数学上等价的基,工程意义可能完全不同。Re、Pr、Nu 是几代工程师从无穷多种合法组合中筛选出的'好基',这种工程知识必须显式地告诉框架。"

### 2.4.1 `Group` 类(templates.py:25-52)

- **【计算】** `evaluate(data)`:∏ data[var]^exp 逐行计算群数值,缺列 KeyError;`check_dimensionless(system)`:Σ exp·dimension_vector = 0(atol=1e-9)验证群确实无量纲,引用未声明变量抛 KeyError。
- **【原因】** 用户声明的无量纲群是**工程约束**——框架的职责是验证合法性并投影数据,而不是替用户发明基。
- **【导向】** `fit_with_templates` 的输入;`validate_pump.py`(Q/H/P 群)、`validate_friction.py`(f·Re 群)、`validate_airfoil.py`(BPM 先验群)都以 Group 声明工程假设。

### 2.4.2 `fit_with_templates(data, system, target_group, feature_groups)`(templates.py:55-142)

```python
overlap = {g.name: sorted(target_vars & set(g.exponents)) for g in feature_groups}
if overlap:
    warnings.append(f"目标群 '{target_group.name}' 与特征群共享变量 {overlap}:...R² 不宜直接比较")
```

- **【计算】** 三步:①所有群 `check_dimensionless`(非法指数 → 结构化失败);②共享变量检测——目标群与特征群变量集求交,重叠写入 `warnings`(教科书情形 Nu∩Re={D} 也会提示,措辞是提示性的);③`Group.evaluate` 计算数值 → 有限值/正值防线 → 共享 `_ols` 引擎做 log-log 拟合 → 返回与 `fit_power_law` 同构的结果(含 std_errors/design_rank/condition_number/warnings)。
- **【原因】** "工程上推荐使用模板路径"(README):自动 π 路径的基选择有任意性,模板路径把领域知识变成显式输入。共享变量提示的措辞刻意是"不宜直接比较"而非禁止——Re 与 Nu 共享 D 是教科书常态。统计口径复用共享引擎,注释明确"避免与 fit_power_law/PowerLawFamily 的统计定义漂移"。
- **【导向】** 三个 validate_*.py 正/负对照脚本的核心拟合调用;`TestStdErrorsParity` 测试锁定它与 fit_power_law 的标准误逐位一致;`exponents`/`std_errors` 键直接被泵脚本判"恢复指数落入 95% t 区间"。

---

## 2.5 本部分小结:发现层的证据链设计

```
_prepare_pi(公共前半段)
 ├─ check_data            → Level 门控(数据够不够)
 ├─ dimensionless_groups  → 目标隔离/化简(可解释性)
 ├─ compute_groups        → π 数值表
 ├─ 逐特征降级            → dropped_features(零是物理,NaN 是事故)
 ├─ 共享变量检测          → shared_vars(π 群伪相关警示)
 └─ form_check            → 形式失配报警(统计精确 ≠ 形式正确)
        ↓ 三入口分流
 power_law / multi(R² 排序) / auto(k-fold CV + 一标准误差 + AICc/BIC)
        ↓
 _promote_level_3(CV R²>0.9 才升级,矛盾写 warnings)
```

每个分支的存在都有对应的"事故记忆":逐特征降级来自 airfoil 的 α=0;共享变量检测来自泵/教材情形;form_check 来自 PK 负对照;Level 3 统一判据来自同一次负对照。代码不是抽象设计,是验证史的直接沉淀。
