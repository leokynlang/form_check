# 第一部分 核心引擎:量纲分析与最小二乘

> 本部分覆盖 `core.py`(438 行)、`_ols.py`(85 行)。这两个文件是全框架的数学地基:前者把"带单位的物理变量"变成"无量纲群",后者是全框架**唯一**的 log 空间最小二乘实现。报告中每个计算按 **【计算】做什么 → 【原因】为什么存在 → 【导向】结果流向哪里** 三段组织。

---

## 1.1 `core.py` — 量纲引擎

**文件定位**:定义 `Variable`(带单位/量纲/角色/误差的变量"身份证")与 `System`(变量集合 → 量纲矩阵 → π 群)。模块 docstring 声明三条设计哲学:变量必须显式声明单位和角色(框架不猜);量纲检查在数据进入框架的最早阶段进行;数据不足时显式报告而不是默默处理。

### 1.1.1 模块级常量

```python
_ureg = pint.UnitRegistry()
SI_BASE_DIMS = ['[mass]', '[length]', '[time]', '[temperature]',
                '[current]', '[substance]', '[luminosity]']
```

- **【计算】** 创建全局共享的 pint 单位注册表;定义 7 个 SI 基本量纲的顺序。
- **【原因】** 所有变量共享一个 registry,避免多 registry 之间的单位转换混乱;每个变量的量纲被表达为这 7 维上的整数指数向量(单位 ≠ 量纲:kg/m³ 与 g/cm³ 单位不同、量纲相同)。
- **【导向】** `SI_BASE_DIMS` 的顺序决定 `Variable.dimension_vector()` 返回向量的分量含义,进而决定量纲矩阵每一行的物理含义;`_ureg` 被 `Variable.__post_init__` 用于单位解析。

### 1.1.2 `_rational_null_space(matrix, max_denominator=10000)`(core.py:33-93)

```python
rational = [[Fraction(float(value)).limit_denominator(max_denominator)
             for value in row] for row in matrix]
reconstructed = np.array([[float(value) for value in row] for row in rational], dtype=float)
if not np.allclose(reconstructed, matrix, rtol=0.0, atol=1e-12):
    return None
```

- **【计算】** 把量纲矩阵 D 的每个浮点元素用 `Fraction.limit_denominator` 有理化,重建浮点矩阵并与原矩阵逐元素比对(atol=1e-12)确认有理化无损;然后做**精确有理高斯-约当消元**(core.py:59-83):逐列选主元、行交换、主元归一、其余行消去,记录主元列;最后对每个自由列构造零空间基向量——自由列分量取 1,主元列分量取 `-rational[row][free_col]`(core.py:85-93)。
- **【原因】** pint 的量纲指数通常是整数或小分母有理数。数值 SVD(`scipy.linalg.null_space`)能找到正确的子空间,但其正交归一基是任意旋转的结果——指数形如 0.7071/-0.4082,完全不可读。精确有理消元给出形如 `Re = ρ¹v¹D¹μ⁻¹` 的可读指数;无法忠实有理化时(自定义量纲)返回 `None` 走数值回退,而不是硬猜。
- **【导向】** 返回值被 `System.dimensionless_groups()`(core.py:313)消费:非 `None` 时走精确路径生成 π 群,并在返回前做防御性验证 `np.allclose(D @ vector, 0)`(core.py:326-327);`None` 时走 SVD 数值回退。

### 1.1.3 `_isolate_target_exact(basis, target_index)`(core.py:96-118)

- **【计算】** 在精确有理零空间基上做目标隔离:找出所有 target 分量非零的基向量,取 target 分量绝对值最大者为主元,缩放使其 target 指数恰为 +1;其余向量减去自身 target 分量 × 主元向量,使它们不再含 target(core.py:101-118)。
- **【原因】** 没有隔离时,零空间基中多个向量都含 target,输出会变成 `h^11·…` 这类不可解释组合,甚至"用 target 预测 target"。隔离后:恰好一个群含 target 且指数归一(+1,如 Nu = h·D/k),其余群不含 target——保证回归中 Y 与 X 的分离是物理可解释的。
- **【导向】** 被 `dimensionless_groups(target=...)` 调用(core.py:318);其输出决定哪个 π 群成为回归的因变量 `y_col`、哪些成为特征 `x_cols`(`_prepare_pi` 的 target_groups 逻辑,discover.py:139-146)。

### 1.1.4 `_simplify_target_exact(basis, target_index, coefficient_limit=6)`(core.py:121-173)

```python
def score(vector):
    nonzero = [abs(float(value)) for value in vector if value != 0]
    return (len(nonzero), sum(nonzero), max(nonzero, default=0.0))
```

- **【计算】** 利用"目标隔离不唯一"的自由度做坐标下降搜索:任何不含 target 的群可以整数倍加减到目标群上而不破坏 target 指数与无量纲性。以(非零分量个数、绝对值之和、最大绝对值)三元组为字典序评分,在 coefficient_limit=±6 的整数系数范围内反复尝试"目标群 ± c·特征群",直到无法再降低评分(core.py:141-169)。
- **【原因】** 隔离后的目标群合法但不唯一,可能仍携带冗余变量。工程上 Nu=hD/k(3 个变量)远比 h⁰·⁵·D²·k⁻¹·ρ·…(7 个变量)可读。用精确有理运算搜索"变量更少、指数更小"的等价表达,把教科书式的简洁群还给用户,全程无舍入误差。
- **【导向】** 在 `dimensionless_groups` 中紧跟 `_isolate_target_exact` 之后调用(core.py:320-321);最终决定用户看到的 π 群公式形态,并影响 `compute_groups` 的列命名(`PI_i_var^exp_...`)。

### 1.1.5 `_canonical_fraction_vector(vector, target_index=None)`(core.py:176-201)

- **【计算】** 把精确有理向量规范成"标准指数向量":若 target 分量非零,全向量除以 target 分量(归一到 +1,core.py:181-184);否则取所有分母的最小公倍数通分成整数,再除以最大公约数、并翻转符号使首个非零分量为正(core.py:186-201)。
- **【原因】** 同一个零空间向量的 (2,4,-6) 与 (1,2,-3) 是同一个 π 群,必须唯一化;target 归一与"首非零分量为正"让群定义确定、可复现、可比较。
- **【导向】** `dimensionless_groups` 中对每个基向量调用(core.py:322-325),产出的浮点指数字典直接构成返回的 groups 列表。

### 1.1.6 `Variable`(dataclass,core.py:204-244)

- **【计算】** `__post_init__` 立即用 `_ureg.parse_expression(self.unit)` 解析单位字符串,失败折叠为 `ValueError`(core.py:225-230);校验 role ∈ {state, material_property, geometry, target, derived}(core.py:232-233)。`dimension_vector()` 从 pint 的 dimensionality 字典取出 7 维指数向量(core.py:235-241)。
- **【原因】** "没有单位的数字、没有角色的列、没有误差范围的测量,都不被接受"——这是框架"不猜"哲学的第一道闸。单位解析失败比延迟失败好(错误在声明处暴露,而非在发现流程深处);role 集合校验保证下游对"谁是 target"的判断有据可依。
- **【导向】** `dimension_vector()` 的输出按列堆叠成 `System.dimension_matrix()` 的量纲矩阵 D;`uncertainty` 字段(相对误差)由 `uncertainty.py::propagate_uncertainty` 消费(见第三部分 3.2)。

### 1.1.7 `System.add / add_variable`(core.py:263-272)

- **【计算】** 向变量字典添加变量,重名抛 `ValueError`;`add_variable` 是免构造 `Variable` 的便利方法。
- **【原因】** 变量集合是命名空间,重名意味着量纲矩阵列的歧义,必须在声明期拒绝。
- **【导向】** `self.variables` 是后续所有计算(量纲矩阵、check_data、π 群)的唯一变量来源;`add_variable(..., uncertainty=u)` 的 kwargs 透传给 `Variable`。

### 1.1.8 `System.dimension_matrix()`(core.py:274-285)

- **【计算】** 把每个变量的 7 维量纲向量按声明顺序列堆叠为矩阵 D(行=基本量纲,列=变量),同时返回列名顺序。
- **【原因】** Buckingham π 定理的代数对象:无量纲组合 ⇔ D 的零空间中的向量。D 是 π 群计算的唯一定义性输入。
- **【导向】** 传入 `dimensionless_groups` → `_rational_null_space`。

### 1.1.9 `System.dimensionless_groups(target=None)`(core.py:287-364)

```python
exact_basis = _rational_null_space(D)
if exact_basis is not None:
    ...
    if all(np.allclose(D @ vector, 0.0, rtol=0.0, atol=1e-9) for vector in vectors):
        return [{var_names[j]: float(vector[j]) ...} for vector in vectors]
```

- **【计算】** 主路径:精确有理零空间 → 目标隔离 → 目标群化简 → 规范化,并在返回前验证每个群满足 D@x=0(atol=1e-9)。回退路径(core.py:336-364):SVD 数值零空间 → 数值版 `_isolate_target`(target 指数归一,否则按最大绝对值缩放并定号)→ 同样的 D@x=0 验证,失败抛 `RuntimeError`(这是全文件唯一抛运行时异常的位置——生成的"π 群"连无量纲性都不满足属于内部不变量破坏,不可静默)。target 不在系统中抛 `ValueError`(core.py:307-308)。
- **【原因】** 这是 π 定理的直接实现:零空间的每一维给出一个独立的无量纲群;目标隔离与化简保证输出的物理可解释性(见 1.1.3/1.1.4);返回前的验证防止任何路径(尤其数值回退)输出非法群。
- **【导向】** 群定义字典列表流向四处:① `compute_groups`(discover.py:24)数值化;② `_prepare_pi` 的 y_col/x_cols 定位与共享变量检测;③ `uncertainty.propagate_uncertainty` 的固定群重投影;④ 直接呈现给用户(README 快速入门的输出)。

### 1.1.10 `System._isolate_target(ns, t)`(core.py:366-386)

- **【计算】** 数值回退路径的目标隔离:选 target 行绝对值最大的列作主元,除以该分量使 target 指数=1;其余列减去 target 行分量×主元列。若所有列的 target 分量 <1e-12(target 本身无量纲)则原样返回。
- **【原因】** 与 1.1.3 相同的数学目标,作用在浮点矩阵上;是精确路径不可用时的镜像实现。
- **【导向】** 仅被 `dimensionless_groups` 的数值回退分支调用。

### 1.1.11 `System.check_data(data)`(core.py:388-435)

```python
if report['n_samples'] < 10 * n_vars:
    report['warnings'].append(f"样本数 {report['n_samples']} 相对于变量数 {n_vars} 偏少,...")
else:
    report['level'] = 2
```

- **【计算】** Data Readiness Assessment:报告缺失变量(declared − present)、多余变量(present − declared)、样本数;可信度分级——缺变量 Level 0,齐备 Level 1,样本数 ≥ 10×变量数 则 Level 2(core.py:416-434)。
- **【原因】** 这是可信度分级(L0-L5)的第一关:单位/量纲检查通过、变量齐全,才有资格升到更高等级。10×变量数是"支持稳健拟合"的粗启发式——太少的样本连 π 群拟合的自由度都不够,框架选择提前警告而不是让回归悄悄失败。
- **【导向】** 返回的 report(dict:`level/missing_vars/extra_vars/n_samples/warnings`)被 `_prepare_pi` 用作 `data_report`(Level<2 直接显式降级,discover.py:122-126);`warnings` 列表后续被逐特征降级、共享变量、form_check 三处追加写入,最终透传给所有 discover 入口的调用方与 CLI 打印。

---

## 1.2 `_ols.py` — 共享 OLS 引擎

**文件定位**:模块 docstring 讲清了它存在的历史原因——fit_power_law、fit_with_templates、PowerLawFamily、anomaly 的基准拟合、scenarios 的内层拟合曾经各自维护一份几乎相同的 `lstsq + R² + 标准误` 代码,每复制一份,统计补丁(VIF、岭、自由度修正)就要同步改多处,并用对等性测试人工锁定漂移。收拢到一个模块后,各调用方只保留自己的 API 形态。

### 1.2.1 `ols(design, response)`(_ols.py:26-75)

```python
coef, _, rank, _ = np.linalg.lstsq(design, response, rcond=None)
...
if n > p:
    sigma2 = ss_res / (n - p)
    cov = sigma2 * np.linalg.pinv(design.T @ design)
    std_errors = np.sqrt(np.diag(cov))
else:
    sigma2 = float("nan")
    std_errors = np.full(p, np.nan)
```

- **【计算】** 一次普通最小二乘返回**全部**诊断:
  - 系数 `coef`(SVD 基 lstsq)、拟合值 `prediction_log`;
  - `ss_res`/`ss_tot`/样本内 `r_squared = 1 − ss_res/ss_tot`(ss_tot=0 时记 0);
  - `rank`(设计矩阵秩)、`condition_number = np.linalg.cond(design)`;
  - 系数标准误:经典公式 `σ² = ss_res/(n−p)`,`Cov = σ²·(XᵀX)⁺`,`se = √diag(Cov)`;n≤p 自由度为负时显式 NaN,不假装有标准误(_ols.py:53-60)。
- **【原因】** 统一口径:全框架所有 log 空间拟合的 R²、标准误、条件数必须同源同值,否则"两个模块报出的标准误差 1e-9"这类漂移会污染下游所有推断。`TestStdErrorsParity` 测试(tests/test_discover.py)锁定迁移前后逐位一致。rank<p 时仍按伪逆给标准误(与迁移前行为一致),可辨识性判定交给调用方(identifiable 键)。
- **【导向】** 被五方消费:`fit_power_law`(discover.py:74)、`fit_with_templates`(templates.py:112)、`PowerLawFamily._fit`(families.py:140)、`fit_with_guard` 的基准 OLS 与 exclude 重拟合(anomaly.py:343,437)、`_fit_linear`/`_fit_strategy`(scenarios.py:61-68,84-144)、`PiecewisePowerLawFamily` 的 `_logspace_ols`(families.py:177-180)。`std_errors` 键最终出现在泵相似律验证、不确定度 MC 的回归对照里。

### 1.2.2 `ols_r2(design, response, coef)`(_ols.py:78-85)

- **【计算】** 对**任意**给定系数(OLS/Huber/岭均可)计算样本内 R²,不重新求解。
- **【原因】** 守卫层的 before/after 对比(anomaly.py:271-291 的 `_logspace_candidate`)需要用"调整后的系数"在原设计矩阵上算可比 R²——直接复用 ols() 会把系数覆盖掉。
- **【导向】** 仅被 `anomaly._logspace_candidate` 消费,写入 before/after 候选的 `r_squared` 字段,最终进入 `fit_with_guard` 的 comparison 报告。

---

## 1.3 本部分小结:一条数学主轴

两个文件合起来是框架的"公理层":

```
Variable(单位声明)
   → System.dimension_matrix()        # D 矩阵
   → _rational_null_space(D)          # 精确零空间(可读指数的来源)
   → _isolate_target_exact / _simplify_target_exact   # 物理可解释性
   → compute_groups()                 # π 数值表
   → _ols(log π)                      # 幂律指数 = 零空间上的回归系数
```

每一环的"原因"都对应框架哲学的一条:不猜(单位必须声明)、可解释(目标隔离+化简)、显式降级(有理化失败→数值回退→验证失败→RuntimeError,绝不带病输出)、统计同源(唯一 OLS 实现)。
