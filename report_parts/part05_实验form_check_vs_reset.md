# 第五部分 实验一:form_check vs Ramsey RESET 基准(`experiments/form_check_vs_reset/`)

> 本部分覆盖 15 个 Python 文件(约 3,300 行):`common.py`(302)、`scenarios.py`(319)、`run_experiment.py`(187)、`analyze.py`(301)、`run_srsd.py`(184)、`srsd_data.py`(217)、`srsd_analyze.py`(255)、`run_families.py`(248)、`families_analyze.py`(192)、`battery.py`(132)、`run_battery.py`(170)、`battery_analyze.py`(162)、`battery_probe.py`(212)、`battery_probe2.py`(147)、`make_report.py`(212)。
>
> **实验总目标**(论文主张 1 的证据链):量化"现有指标(R²/CV/AICc/std_error)集体静默放行 fitted-but-wrong 模型"的现象,把 dimreflex 的 `form_check` 与统计学标准答案 **Ramsey RESET(1969)** 的三个变体在同一公平尺度上逐场景对比,最终按数据改写论文 claim("盲区制图 + 组合门卫",不再宣称单一检测器优越)。实验纪律在 `common.py` docstring 中不可违背:方法 A 调用 form_check **本体**不复制不暗改;任何失败写行内 failure 字段绝不静默丢弃;本目录外零写入。

---

## 5.1 `common.py` — 公共实验台(公平性设计的核心)

### 5.1.1 常量与环境快照

- `ALPHA = 0.05`、`B_PERM = 999`:所有方法统一的显著性水平与排列次数(p 值分辨率 1/1000)。**原因**:修掉旧 pilot"直接比 |t|>3 vs p<0.05"的不公平——两种量纲不同的判据不可比,必须全部折算到同一 p 值尺度。
- `sha256_file(path)`(common.py:41-46):1MB 分块流式 SHA-256。**原因**:核心代码审计。
- `capture_env()`(common.py:49-82):快照 python/numpy/scipy/pandas/dimreflex 版本、平台、`discover.py`/`_ols.py`/`validate_friction.py`/`__init__.py` 四个核心文件的 sha256、git HEAD 与工作区脏条目数、实验脚本自身与 scenarios.py 的 sha256。**原因**:结论必须能定位到"哪一版核心代码上得出"——git 不可用时也记入元数据不中断实验。**导向**:每个 run 目录的 `env.json`,进入论文报告 §8 复现信息。

### 5.1.2 `_ols(X, y)`(common.py:86-93)

- **【计算】** 实验自带的迷你 OLS:lstsq 系数、残差、RSS、hat 对角线(`einsum`)、`(XᵀX)⁺`。
- **【原因】** 实验脚本刻意**不 import dimreflex._ols**——只读调用被试对象 form_check 本体,其余统计自己实现,保证被试代码零改动(纪律"只读 dimreflex")。
- **【导向】** `form_check_with_p`、`reset_test`、`extrapolation_error`、`run_srsd.eval_errors` 全部复用。

### 5.1.3 `form_check_with_p(data, target, features, seed, n_boot=200, b_perm=999)` — 方法 A(common.py:107-172)

```python
fc = dimreflex_form_check(data, target=target, features=..., random_state=seed)
...
pinv_row = np.linalg.pinv(design)[-1, :]          # c2 = pinv_row @ y
idx = np.argsort(rng.random((b_perm, n)), axis=1)
null_c2 = c2_base + resid[idx] @ pinv_row
p_emp = (1 + #{|c2_perm − mean| ≥ |c2_obs − mean|}) / (b_perm + 1)
```

- **【计算】** ①调用 `dimreflex.form_check` 本体(n_boot=200,原生 |t|>3 阈值只留痕);②包装层用同一设计重算 c2_obs,与本体返回的 curvature 比对(rtol=1e-6),不符 → `fail='statistic_mismatch'`(实现漂移哨兵:意味着 discover.py 漂移、排列零分布失效);③排列零分布向量化——**c2 是 y 的线性泛函**,故 `c2(pred+perm) = pinv末行 @ (pred+perm)`,999 次排列一次矩阵乘完成(与逐次 lstsq 差 ~1e-13 相对);④零残差退化(零噪声+精确幂律)→ `fail='degenerate_zero_residual'`(检验对象不存在);⑤经验 p 值 + 拒绝判定 `p < 0.05`。
- **【原因】** 公平性的关键装置:四方法在**同一 α 尺度**上比较。一致性自检保证"被试没有在实验期间被改动"是可验证的命题而非假设。
- **【导向】** 行内字段 `p_fc/rej_fc/fc_t_native/fc_t_emp/curvature/fail_fc/time_fc` → raw.parquet → 一切下游分析。

### 5.1.4 `reset_test(Xlog, ly, powers=(2,3), variant, vcov)` — 方法 B/C/D(common.py:176-252)

```python
if variant == "fitted":
    added = np.column_stack([yhat ** k for k in powers])      # Ramsey(1969) 教科书形式
else:
    added = np.column_stack([Xlog[:, j] ** k ...])            # exog 变体
...
F = ((rss_r − rss_u)/q) / (rss_u/(n − k))                      # 经典 F
# HC3: V = (X'X)^-1 [Σ (e_i/(1−h_i))^2 x_i x_i'] (X'X)^-1
W = b_u[-q:] @ solve(Vq, b_u[-q:]);  pval = stats.f.sf(W/q, q, n−k)
```

- **【计算】** Ramsey RESET 全家:受限模型 `[1, log X]` → 无约束模型再加 q 个增广列(fitted:yhat²、yhat³;exog:每个特征列的 2、3 次幂);三种统计量:
  - **fitted_F**(方法 B):经典 F 统计量,`F(q, n−k)` 参考分布;
  - **fitted_HC3**(方法 C):HC3 夹心协方差的 Wald 统计量 W/W 参考 `F(q, n−k)`(与 R 的 `car::linearHypothesis` 稳健检验同口径);
  - **exog_HC3**(方法 D):exog 增广 + HC3。
  - 七种失败显式登记:`rank_deficient / insufficient_df / degenerate_zero_residual / degenerate_perfect_fit / numerical_nonfinite / hc3_leverage_degenerate(1−h ≤ 1e-12) / hc3_singular_covariance`,外加 `cond`(无约束设计条件数)留痕。
- **【原因】** RESET 是形式失配检验的统计学标准答案,三个变体覆盖三个正交维度:fitted vs exog(增广列来源)、经典 vs HC3(协方差估计)。HC3 变体的存在是为了回答"异方差下经典 F 超尺寸,稳健协方差能否修复"(分析阶段的 Q3)。
- **【导向】** `_RESET_TAGS` 三个配置实例化;`run_all_methods` 逐方法写入 `p_/rej_/stat_/fail_/cond_/time_` 字段。

### 5.1.5 `run_all_methods(df, target, features, seed)`(common.py:260-277)

- **【计算】** 对一个数据集跑 A+B+C+D,返回宽行 dict,每方法独立计时、独立失败字段。
- **【原因】** 单一入口保证每个数据集上四方法见到**完全相同**的输入(同数据同 seed),差异只可能来自方法本身。
- **【导向】** run_experiment/run_srsd/run_families/run_battery 共同的最内层调用。

### 5.1.6 `extrapolation_error(df, target, features, truth_logmean, far=3.0, m_pts=12)`(common.py:281-302)

- **【计算】** 训练域内 OLS → 沿第一特征在 log 域边界外推 ×far(低端/高端各布点,其余特征固定在中位数)→ 真值由场景注册的闭式真函数给出 → 相对误差 `|exp(pred−true)−1|`(log 差稳定形式防溢出)→ 返回(中位误差, 最远点误差)。
- **【原因】** 回答第四个实证问题:"报警能否预测外推危险?"——统计精确但形式错误的模型的代价集中显现在域外,把检验统计量与外推误差关联起来是"报警有没有用"的落地判据。
- **【导向】** `run_experiment.one_cell` 的 extrap.parquet;analyze.py 图 5 与外推 AUC;srsd_analyze 的 Q4。

---

## 5.2 `scenarios.py` — 19 个预注册场景库

**设计原则**(docstring):每个场景 = `seed -> DataFrame(target, x1, [x2])` 函数 + 元数据;真函数经 `y = exp(f_true(lx))` 进入,噪声以相对尺度加在 log 空间(乘性进入 y,与 validate_friction 语义一致)。

### 5.2.1 公共构件

- `_finalize(ly, rng, rel)`(scenarios.py:18-23):`ly += rel·N(0,1)` → `y = exp(ly_c)`,同时返回干净真值 `exp(ly)`。**原因**:噪声在 log 空间加性 = 原空间乘性,保证"相对噪声"语义;返回 clean 值供外推真值对照。
- `_domain(n, seed, ...)`(scenarios.py:25-31):`rng = default_rng(seed·7919 + 104729)`,x1∈log-uniform[0.1,10]、x2∈log-uniform[0.05,2]。**原因**:量程跨度 ~4.6 与 pilot/validate_friction 同量级(曲率零分布 σ ∝ 量程跨度² 的教训——域必须可比);同 seed 全场景同网格、噪声流独立(控制变量)。
- `_reg` 装饰器 + `SCENARIOS` 注册表 + `SCENARIO_META`(scenarios.py:65-319):`(类别 H0/H1, 族, 真值函数, 描述)` 四元组注册。**原因**:**预注册**——H0/H1 分类在代码写入时固定,不依据检验结果调整(run_experiment 的 config.json 里原文声明)。真值函数单独实现(不以生成器代码复用)是刻意的:外推真值用固定域中点而非样本均值(注释:`_t_thresh` 用 `_M` 非 `ctx["m1"]`),避免评估随样本抖动。

### 5.2.2 H0 场景(严格单项式真值,拒绝=误报)

| 场景 | 真值 | 特殊设计 |
|---|---|---|
| `H0_BASELINE` | y=1.5·x1^0.7 高斯噪声 | 基线校准格 |
| `H0_MULTIVAR` | y=1.2·x1^0.5·x2^0.9 | 多特征基线 |
| `H0_CORR` | 同上,**特征相关 ρ=0.9** | log 空间线性混合 `lx = 2z − lz` 构造相关设计,独立噪声流(seed+31) |

- **【原因】** 三个 H0 覆盖校准的三个压力面:单变量、多变量、病态设计。`H0_CORR` 专门压 HC3-Wald 的共线脆弱性(后来实锤为 42-56% 超尺寸的来源之一)。

### 5.2.3 结构 H1 场景(均值结构失配)

| 场景 | 真值结构 | 设计要点(注释原义) |
|---|---|---|
| `QUAD_LIN` | log 空间二次曲率 c=+0.10 | 最基础的曲率信号 |
| `CUBIC_LIN` | log 空间三次(标准化)0.05 | 三次项——后来发现 form_check 的"随 n 反降"盲区 |
| `INTERACT` | lx1·lx2 交叉项 0.30sd | "两种 RESET 都不建模的 power 缺口" |
| `THRESHOLD` | 指数 0.4→1.0 阈值(域中点) | 分段结构,四方法皆弱(不在平滑方向张成空间) |
| `EXP_MODEL` | y=exp(0.5·x1) | log 空间强凸 |
| `LOG_MODEL` | y=0.8+1.1·ln x1(域平移 +2.4 保正值) | 凹形;真值函数有 `_t_logm_fixed` 平移修正 |
| `RATIONAL` | y=x1/(0.2+0.1·x1) 饱和 | log 空间凹;fitted-RESET 的结构性盲区 |
| `OMITTED_M` | x2 与 x1 相关且**二次**进入,拟合只见 x1 | 注释的方法学要点:"x2 独立线性进入不构成形式失配(E[ly|x1] 仍线性)"——所以必须让遗漏变量二次进入 |
| `DILUTION` | x1 曲率 0.10,x2 无 | form_check 平均平方列的信号稀释实验 |
| `OPPOSING` | +0.12 与 −0.12 反向曲率 | **对抗性场景**:"平均曲率列的总体系数精确为 0,检验统计量期望携带零信号而真实失配存在"——form_check 的理论盲区 |
| `COLLINEAR_CURV` | ρ≈0.95 + 曲率 0.10 | H1 双重压力:信号 + 病态(测失败登记与 p 值稳定性) |

### 5.2.4 噪声层 H1 场景(结构仍是严格幂律——拒绝=对噪声假设的误报)

| 场景 | 噪声设计 | 设计要点 |
|---|---|---|
| `H1_HEAVYTAIL` | Student-t(3)(标准化方差 1) | 重尾下排列零分布稳健性 |
| `H1_VAR_GROWING` | SD 随 log 预测值线性增大(0.3→1.7×rel) | 异方差 |
| `H1_AUTOCORR` | AR(1) ρ=0.6,burn-in 200 丢弃、单位方差标准化 | "打破排列零分布的可交换性假设 → 测各法 FPR 失真" |
| `H1_OUTLIER` | 2% 单点 8σ 单侧冲击 | 点污染 |
| `H1_OUTLIER_GRP` | 连续 6% 区段 +5σ | 成组污染 |

- **【原因】** 噪声层场景与结构 H1 在分析里**严格分开**(analyze.py 的 NOISE_H1 集合):"它们的拒绝对'形式'而言是误报,不能与结构 H1 功效混在一张总平均里"——防止"总平均宣布胜者"的错误推理。
- **【导向】** 全部场景经 `SCENARIO_META` 进入 run_experiment 的网格与 config.json 的预注册声明;`SYNTHETIC_TRUTH`(battery.py)又引用它们做真值族注册。

---

## 5.3 `run_experiment.py` — 合成台架网格 runner

- **固定网格**:`n ∈ {50,200,1000} × rel ∈ {0, 0.5%, 1%, 3%, 5%} × 19 场景 × seeds`。
- **种子规则**:`base_seed = 100_000 + s·97 + n_idx·13 + rel_idx·5_000_000`(run_experiment.py:39-40)。**原因**:网格位置入种子,任何格可独立复现、跨格不重叠。
- **`ctx_for`**(run_experiment.py:43-55):从训练数据估计真值函数所需的上下文(样本均值/中位数/交叉项标准差)。
- **`one_cell`**(run_experiment.py:58-86):逐种子生成场景 → `run_all_methods`(四方法)→ `rel ≥ 0.5%` 时加算外推误差(零噪声场景残差退化无外推意义,`EXTRAP_REL_MIN=0.005`)→ 行字典带 scenario/klass/family 元数据。
- **断点续跑与落盘**(run_experiment.py:89-177):启动读旧 raw.parquet,该格 `(scenario,n,rel)` 种子数已齐则跳过;每 4 格 flush 到 parquet+CSV;config.json 两次写入(started/finished+wall_seconds),内嵌预注册声明原文。**原因**:500 种子全量 14 万+行需要跨会话分片跑;预注册声明随配置落盘使"先注册后实验"可审计。
- **【导向】** `runs/<tag>/raw.parquet, extrap.parquet, raw.csv, config.json, env.json` → analyze.py。

---

## 5.4 `analyze.py` — 汇总分析(校准/功效/外推/胜负)

### 5.4.1 基础统计件

- `wilson(k, n, z=1.96)`(analyze.py:37-44):Wilson 分数区间 95% CI。**原因**:拒绝率是二项比例,小样本(n=50 种子)时正态近似不可靠,Wilson 区间在边界附近不失真。
- `rate_ci(mask)`(analyze.py:47-51):`k/n, lo, hi, n` 四元组。

### 5.4.2 分片合并(main 入口,analyze.py:54-79)

- `--shards full_r0,full_r1,...`:合并各分片 raw.parquet/extrap.parquet,输出目录用首个 tag。**原因**:500 种子跑不完一个会话,分片是运行单位,分析必须能合并;`os.makedirs(run_dir, exist_ok=True)` 使合并模式输出目录可不存在。

### 5.4.3 逐格统计与经验 FPR 校准(analyze.py:81-115)

```python
for n in sorted(raw.n.unique()):
    sub = h0[h0.n == n]
    tau[m] = float(np.quantile(p, 0.05)) if len(p) >= 50 else np.nan
```

- **【计算】** ①逐 `(scenario,klass,n,rel)` 格:各方法拒绝率+Wilson CI、mean p、平均耗时、失败率;②**经验 FPR=5% 校准阈值**:每 n 用全部 H0 格 p 值的 5% 分位数定 τ;③`power_at_tau`:H1 格上 `p ≤ τ` 的比率。
- **【原因】** docstring 约束 10:"power@empFPR≤5% 回答'相同经验 FPR 下谁更有功'"——名义 α=0.05 下经典 F/HC3 超尺寸(误报率高),直接比拒绝率会让超尺寸者虚赢;校准是公平比较的必需设施。
- **【导向】** summary.json 的 `calibrated_tau_at_fpr5`、`power_at_empirical_fpr5_rel1pct`;论文"fc 校准后功效小中样本 2× 领先(0.41/0.51 vs 0.19/0.28)"的直接来源。

### 5.4.4 五张图(analyze.py:117-227)

1. **校准图** `{tag}_calibration.png`:H0 池 p 值经验 CDF per n,对角线=完美校准。**原因**:p 值一致性是检验有效性的第一道检查。
2. **功效热图** `{tag}_power_heat.png`:11 结构 H1 × 4 方法,n=200/rel=1%,数值标注。**导向**:论文主图。
3. **n 扩缩** `{tag}_n_scaling.png`:5 关键场景(QUAD/CUBIC/INTERACT/DILUTION/OPPOSING)功效 vs n(log 轴),0.8 参考线。**原因**:捕捉"随 n 反降"(CUBIC_LIN 的 form_check 盲区)这类非单调现象。
4. **运行时间** `{tag}_timing.png`:每方法每 n 的平均秒(log 轴)。**原因**:工程可用性(成本维度)。
5. **外推图** `{tag}_extrap.png`:报警 vs 远端外推相对误差散点 + p 值对"误差>10%"的 rank-AUC(Mann-Whitney)。**原因**:Q4——报警预测外推危险的能力是"检测器有用"的独立证据。

### 5.4.5 胜负矩阵与 summary(analyze.py:229-291)

- **【计算】** 逐 `(scenario,n,rel)` 格比较四方法拒绝率,`rates[a] > rates[b]+0.02` 记 a 胜(2pp 容差防噪声平局);summary.json 汇总 H0 误报/结构 H1 功效/噪声层拒绝/校准功效/外推 AUC/胜负矩阵/失败率。cell_stats.csv 与 win_matrix.csv 落盘。
- **【导向】** make_report.py 的输入;runs/full/report.md;PAPER_PLANS 的 gate 判定。

---

## 5.5 `srsd_data.py` + `run_srsd.py` + `srsd_analyze.py` — SRSD-Feynman 真实公式基准(P1)

### 5.5.1 `srsd_data.py` — 官方数据接入(零改动加载)

- `_ensure_upstream()`(srsd_data.py:23-29):浅克隆 `omron-sinicx/srsd-benchmark` 到 %TEMP%。**原因**:用官方任务池而非自制公式,增强外部效度;浅克隆省时。
- `_install_stubs()`(srsd_data.py:32-58):注入 torchdistill/graphviz 桩模块(仅 `make_parent_dirs`/yaml 读取/占位 Digraph)。**原因**:120 个公式类 import 时依赖这两个包的边角工具,公式代码本身一行不改——"零改动加载"与实验纪律一致。
- `_is_monomial(sympy_eq)`(srsd_data.py:74-89):**严格单项式判定**——逐符号计算 `d/ds[ln f]·s`,能全为常数 → 单项式。**原因**:H0 的数学定义("幂律 = log 空间线性")必须机器可判,不能人眼;sympy 病态式异常统一按非单项式处理(保守)。
- `classify(eq_cls)`(srsd_data.py:114-142):**预注册分类**——64 个官方采样点检验有限/正值 → `EXCLUDE(nonfinite/nonpositive)`;变量数>3 → EXCLUDE(harness 限制);`_is_monomial` → H0,否则 H1。注释:"判据只依赖公式结构与官方采样范围,与任何检验结果无关"。**原因**:这是"先分类后实验"的纪律载体——分类函数写死于实验前。
- `_draw(so, rng, n)` 与 `make_dataset_fn(eq_cls)`(srsd_data.py:92-111,165-217):按官方采样类语义(DefaultSampling=log10 域均匀/SimpleSampling=线性/IntegerSampling=整数,正负号语义保留)抽点;`y<=0 或非有限 → 抛 ValueError`(上游数据性质,登记失败不静默);hetero=True 换异方差噪声(SD 0.3→1.7×rel,与合成台架 H1_VAR_GROWING 同构);range_scale 对官方域乘性缩放(敏感性实验用,注释:"2 倍域可能使原本正值函数出负值——那是真实性质,登记 not-applicable,不做平移")。
- **【导向】** task_registry.csv(逐任务分类+理由,交付物);run_srsd/run_families/run_battery 的数据源。

### 5.5.2 `run_srsd.py` — 32 任务网格 runner

- **【计算】** 任务池 = registry 中 klass∈{H0,H1} 的公式(13 H0 单项式 + 19 H1 非单项式);网格同合成台架;每种子每任务:训练集(n 样本)→ 四方法;**test 集**(2000 样本,官方域,同 rel 噪声,独立流 seed+777777)→ log 空间测试 RMSE;**外推集**(log 域 ×3)→ 中位/远端相对误差(真值用闭式 `eq_func` 逐点算)。种子规则 `200_000 + s·97 + ...`;断点续跑同前。
- **【原因】** docstring 四重点问题:合成台架的结论要在**真实物理公式**上复现,否则"合成外推到真实公式族会得出错误的排名"(families_analyze 后来的实锤)。
- **【导向】** `runs/srsd_*/raw_srsd.parquet` → srsd_analyze.py。

### 5.5.3 `srsd_analyze.py` — 合并分析(四重点问题量化)

- `rank_auc(p, y)`(srsd_analyze.py:44-56):Mann-Whitney rank AUC,平均秩处理并列。
- **【计算】** 数据完整性(每格种子数);H0 FPR(Wilson CI,仅 rel≥0.005);经验 FPR=5% 校准阈值与校准功效;逐任务功效(n=200/rel=1%,`h1_power_per_task.csv`);外推 AUC(全局+**按 rel 分层**——注释:"池化会混入'高噪声→高误差但低功效'的 Simpson 伪影");报警 precision/recall;Q3 专用 `hc3_vs_classical_fpr`(fitted_f vs fitted_hc3 vs exog_hc3 按 n);胜负矩阵;Q2 数据驱动:`q2_gap` 找 exog-HC3 比 form_check 优势最大的前 3 多变量任务(正则 `i\.26|i\.27|...|ii\.|iii\.` 筛交互类)。
- **【原因】** docstring 把四个研究问题(Q1 低自由度功效/Q2 平均稀释盲区/Q3 HC3 异方差修复/Q4 报警-外推关联)逐一映射到具体统计量——分析不是图表堆砌,是问题的量化落点。
- **【导向】** `summary_srsd.json + report_srsd.md(SRSD_SECTION) + srsd_calibration.png + srsd_per_task_power.png`;500 种子正式跑结论(RESET-HC3 超尺寸 42-56%;fc 校准后小中样本 2×;外推告警 fc 高召回 0.63 / RESET-F 高精确 0.56)。

---

## 5.6 `run_families.py` + `families_analyze.py` — 函数族分层覆盖(P2)

### 5.6.1 `run_families.py`

- `FROZEN`(run_families.py:31-33):`powers=(2,3), B=999, α=0.05, rel∈{0.5,1,3,5}%, n_train=200, n_seeds=30, range_scales=[1.0,0.5,2.0]`。**原因**:"P1 冻结配置,P2 不得重调"——跨实验可比性。
- `_monomial_exponents(f)`(run_families.py:37-56):廉价单项式判据(顶层 Mul 分解,每因子须为 Symbol 或 Symbol**常数),无需 sympy simplify。
- `classify_family(sympy_eq)`(run_families.py:59-102):按表达式结构分 8 族(exact_power_law/trigonometric/logarithmic/exponential/rational/multivar_additive/additive/fractional_power/other),`preorder_traversal` 扫描节点类型,判据透明写死先于实验。
- `build_family_table` + `select_representatives`(run_families.py:105-136):全部 120 方程建档;分层选择每族至多 6 个(族内排序 eligible 优先→变量数少→名称,"避免 EXCLUDE 成员挤占同族可检方程名额"),总数超限按族轮转裁剪。**原因**:30 个代表覆盖 8 族,兼顾覆盖度与算力。
- `eval_cell` + main(run_families.py:139-244):网格 = H1 代表方程 × rel × {homo,hetero} × 3 种域缩放;H0 族只跑官方域("它们的角色是各族 FPR 校准,不测 scale 交互");种子 base=300_000+格序号×1e6;data_fail 登记照旧。
- **【导向】** `family_table_all120.csv / selection_table.csv / raw_families.parquet` → families_analyze;selection_table 还被 run_battery 复用(冻结复用)。

### 5.6.2 `families_analyze.py`

- **【计算】** 逐方程检测表(per_equation.csv);族汇总(official_homo/official_hetero/scale_half/scale_double 四情形);**N/A 注册**(数据非正/非有限整体失败的格,`not_applicable_cells.csv`——"不平移数据"的审计);采样范围敏感性(>15pp 登记为"随采样范围明显变化");经验 FPR=5% 校准功效(分 homo/hetero,用同格 H0 池校准);报告正文直接嵌入四个重点问题的**人工解读**(哪些族被 fc 漏检:exponential 0.11;被 fitted-F 漏检:rational 0.13-0.19——"yhat²/yhat³ 作为有理函数的局部泰勒展开恰好吸收饱和曲线的主弯曲";量纲诚实留白;采样稳健性)。
- **【原因/导向】** P2 的存在理由写进"冲突/新现象"一节:合成台架排名(fitted-F 最强)在真实公式上反转——"合成外推到真实公式族会得出错误的排名,这是 P2 存在的理由";同时发现 form_check 的 OPPOSING 盲区在真实方程池未见("平均平方列的盲区需要两个变量曲率恰好抵消,真实物理公式中未见此结构")。

---

## 5.7 `battery.py` + `run_battery.py` + `battery_analyze.py` + 两个探针 — 形式合格性电池(P3)

### 5.7.1 机制与动机(battery.py docstring 原文要点)

- **【机制】** 每个可线性化函数族 F 一组坐标 `(u(x), v(y))`:`power=(log,log)`、`exponential=(恒等,log)`、`logarithmic=(log,恒等)`、`rational=(1/x,1/y)`;F 的检验 H0 = "条件均值在 F 坐标里线性";统计量 = [1,u,u²] 设计下二次列系数;零分布 = **wild bootstrap**(Rademacher 翻转残差),经 pinv 末行向量化。
- **【为什么 wild 不是排列】**(探针 1→探针 2 的因果):battery_probe.py 用排列零分布实锤 `v(y)=y 与 v(y)=1/y` 坐标下乘性噪声残差异方差、可交换性破裂 → **rational 假升级 30-34%**;battery_probe2.py 换 wild bootstrap 修复到 4.5-9%。battery.py docstring:"坐标决定零分布合法性,不是实现细节"。两个探针的 docstring 都带**预注册预测**(P1-P5 / W1-W5,运行前写死用于对答案)。
- **【语义】** 输出 `viable={F: p_F ≥ α}` 候选列表而非单一赢家(项目宪法);三种错误分开计数不合并:假排除(真族∉viable/空集)、假包含(错族∈viable)、正确升级(truth_family=none 时空集)。多重性声明:每族各自 α=0.05 不做族间校正——"非拒绝语义下多重性不膨胀第一类错误"。
- `spec_test_wild`(battery.py:63-112):u 逐列变换→标准化→quad=(u²).mean(axis=1)→lstsq 统计量→Rademacher 零分布→经验 p;失败五类显式登记(nonfinite/nonpositive/transform_failed/transform_nonfinite/degenerate_*)。
- `battery_run`(battery.py:115-132):四族各跑(种子 `seed·1_000_003 + 族序号·7919 + 17` 格间不重叠),`viable` 字符串 + `n_viable`。
- `SYNTHETIC_TRUTH`(battery.py:44-53):合成场景真值族注册(先于正式跑写死)——8 个噪声层+3 个 H0 场景=power;EXP/LOG/RAT 各归其族;9 个结构场景=None("线性化加性形式装不下→升级是正确结局,同时量化假包含")。

### 5.7.2 `run_battery.py`

- **【计算】** 配置冻结(n=200, 500 种子, rel∈{0.5,1,3,5}%, B=999 wild);合成侧 19 场景×rel(无 hetero 维——"噪声型已编为独立场景");SRSD 侧**复用 P2 的 selection_table.csv**("选择先于检验结果,冻结复用");种子 `900_000_000 + cell_idx·10_000_000 + seed·97`;按 cell_id 断点续跑;SRSD 侧 truth="none"(H1 方程无电池族——那是 P2 语义不是电池语义)。
- **【导向】** `runs/battery_full/raw_battery.parquet`。

### 5.7.3 `battery_analyze.py`

- **【计算】** 三层汇总:层 1 合成校准(结构 H0+噪声层上 power 族尺寸,池化 vs 0.05 目标);层 2 合成判别(EXP/LOG/RAT 的 containment/假升级/假包含 top3 按 rel);层 2b truth=none 场景(正确升级率+假包含组合分布);层 3 SRSD(每方程升级率/中位|viable|/top viable)。生成 `summary_battery.json + report_battery.md`。
- **【原因】** "语义对齐 battery.py docstring:三种错误分开计数,绝不合并成一个数字"——viable 集合语义下,TPR/FPR 这类单标签指标不适用,必须按错误类型分层。

---

## 5.8 `make_report.py` — 最终报告装配

- **【计算】** 读 summary.json/cell_stats.csv/config.json/env.json/raw.parquet,生成 8 节 report.md:0 公平性设计(约束 4/5/6 的落实——逐条写明)、1 H0 校准表、2 结构 H1 功效表(**逐条解读约束:"不得只报总平均"**——每场景一行文字点评)、3 噪声层 H1(**语义警示**:"本节场景的均值结构是严格单项式——任何拒绝都是对噪声假设的误报,不是发现形式失配;反之 H0 格的不拒绝也不构成'幂律得证'")、4 相同经验 FPR=5% 下的功效、5 胜负矩阵、6 运行时间、7 失败登记(逐类型说明:degenerate_zero_residual=零噪声/精确幂律"检验对象不存在,显式登记"、statistic_mismatch="实现漂移哨兵")、8 复现信息(git HEAD/核心文件 sha256/版本/种子规则);可选 SRSD_SECTION 拼接。
- **【导向】** `runs/full/report.md`(500 种子正式跑的最终报告,论文 §实验 的素材)。

---

## 5.9 实验一小结:因果链与证据流向

```
负对照缺口(C2:幂律静默输出 R²=0.99)
   → form_check 进入产品(discover.py)
   → 需要回答"它与统计学标准答案 RESET 比怎么样"
   → common.py 公平装置(同 p 尺度 + 本体调用 + 漂移哨兵)
   → P1 合成台架(19 预注册场景,网格 3n×5rel×500 seeds)
      + SRSD 真实公式池(32 任务,官方采样)
      + P2 函数族分层(120 方程 8 族)
      + P3 电池(4 可线性化族,wild bootstrap)
   → analyze/srsd_analyze/families_analyze/battery_analyze
      (校准 → 校准功效 → 盲区制图 → 外推关联 → 失败审计)
   → 结论改写 claim:删"单一检测器优越性",
     立"盲区制图 + 组合门卫(OR)+ 校准对照"
```

每个分析脚本都在 docstring/报告里写明**推理约束**(不得只报总平均/噪声层拒绝是误报不是检出/不得依据结果调整分类),这是"预注册 + 审计留痕"方法论的代码化。
