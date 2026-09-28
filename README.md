# dimreflex

**量纲感知的工程模型发现框架**

Dimensional-Aware Engineering Model Discovery Framework

---

## 简介

dimreflex 是一个基于 **Buckingham π 定理**的工程模型自动发现框架。它从带单位的实验/仿真数据中自动推导无量纲群，并在无量纲空间中发现候选函数关系。

**核心理念：**

- 数据进入框架前必须携带单位、量纲、变量角色等元信息
- 函数发现在无量纲空间进行，而非原始变量空间
- 数据不足时框架显式降级，绝不硬输出
- 仅服务工程认知，不接管生产端控制

## 当前状态

> ⚠️ **研究原型 (v0.5.0)**：核心链路（量纲矩阵 → π 群 → 自主模型选择 → 多场景拟合 → 批次监视）已可用并通过合成数据验证；正对照（泵相似律）与负对照（Prandtl–Kármán 对数律，`validation/validate_friction.py`）均已完成。验证与运行脚本统一在 `validation/` 目录。

- 示例数据 `heat_exchanger_data.csv` 是**合成数据**（由 Dittus-Boelter 公式 + 2% 噪声生成），用于验证框架链路，不代表真实实验。
- 工程上推荐使用**模板路径**（`fit_with_templates`，用户显式声明 Re/Pr/Nu 等无量纲群）；`discover()` 自动路径已做目标隔离（target 指数归一化、其余群不含 target），输出物理可解释，但仍建议结合领域知识核对。
- 符号回归等更高级的发现方法尚在计划中。

## 安装

```bash
# 开发安装(推荐,项目已提供 pyproject.toml)
pip install -e .

# 或按兼容范围安装依赖
pip install -r requirements.txt

# 可复现部署:锁定版本(2026-09-04 实测,Python 3.13)
pip install -r requirements-lock.txt
```

## 快速入门

```python
import dimreflex as dr
import pandas as pd

# 1. 定义变量系统
sys = dr.System("pipe_flow")
sys.add_variable("rho", "kg/m**3", "material_property")
sys.add_variable("v", "m/s", "state")
sys.add_variable("D", "m", "geometry")
sys.add_variable("mu", "Pa*s", "material_property")
sys.add_variable("h", "W/(m**2*K)", "target")

# 2. 加载数据
data = pd.read_csv("data.csv")

# 3. 自动发现（自动对 target 做目标隔离）
result = dr.discover(data, sys, target="h")
print(result["fit"]["formula"])
```

**数据格式**：CSV 每一列对应系统中的一个变量，列名与变量名一致，数值为 SI 单位。例如换热器示例（`heat_exchanger_data.csv`）需要列：`rho, v, D, mu, cp, k, h`。

## 架构

```
Variable    — 带单位、量纲、角色的变量对象
System      — 变量集合，构建量纲矩阵，计算无量纲群（支持 target 隔离）
discover    — 端到端发现流程（无量纲化 → 目标隔离 → 幂律/多函数族拟合）
Group       — 用户声明的无量纲群模板（如 Re, Nu, Pr）
fit_scenarios — 比较共享、部分共享和独立的多工况相关式
EngineeringMonitor — 运行包络、分布漂移、残差漂移与连续告警
```

### 函数族体系（families.py）

所有函数族继承同一基底类 `FunctionFamily`，产出同构的 `CandidateModel`
（统一 schema：`{family, formula, params, success, r_squared, n_params,
n_samples, reason}`）：

| 函数族 | 覆盖形式 | 对应改进建议 |
|---|---|---|
| `PowerLawFamily` | y = a·∏x^b（Dittus-Boelter 类） | 现有能力插件化 |
| `PiecewisePowerLawFamily` | 分段幂律 + 变点检测（层流/湍流） | 2.1 第一步 |
| `RationalLogFamily` | log10(y) = 有理式(L)，Gnielinski 类 | 2.1 第二步 |
| `SymbolicRegressionFamily` | PySR 符号回归（占位，显式失败） | 2.1 第三步 |

基底类契约：`fit()` 是模板方法，统一兜底——任何家族的数学异常都被
折叠为 `success=False` 的失败候选，**绝不向调用方抛异常**。

```python
import dimreflex as dr

candidates = dr.discover_multi(data, target="Nu", features=["Re", "Pr"])
for c in candidates:                      # 成功者按 R² 降序在前
    print(c.family, c.success, c.formula or c.reason)

# 注册自定义函数族
@dr.register_family
class MyFamily(dr.FunctionFamily):
    name = "my_family"
    def _fit(self, data, target, features): ...
```

**与量纲引擎合流（推荐入口）**：`discover_multi_pi` 先做量纲分析
（`System` → π 群 → 目标隔离），再把全部函数族跑在 π 空间上——π 空间里
所有运算符天然量纲合法，量纲一致性免费获得：

```python
sys = dr.System("pipe_flow")
sys.add_variable("rho", "kg/m**3", "material_property")
# ... 声明全部变量后:
result = dr.discover_multi_pi(data, sys, target="h")
for c in result["candidates"]:            # 同构候选,含 π 空间公式
    print(c.family, c.success, c.formula or c.reason)

# 或经 discover() 一键走 multi 路径(等价,最优候选在 result["fit"])
result2 = dr.discover(data, sys, target="h", method="multi")
```

多特征支持与诚实声明：

- `PiecewisePowerLawFamily`：变点**仅沿第一个特征**（工况变量，如 Re）
  搜索；每段是完整的多特征幂律（段间指数允许不同）。其余特征不参与
  分段边界判定——这是设计取舍而非遗漏。变点报告 `breakpoint_interval`
  （数据可识别的过渡区间）而非假装精确。
- `RationalLogFamily`：主特征（`features[0]`）上拟合有理函数
  `num(L)/den(L)`，其余特征以幂律协变量项 `b_j·(log10(x_j)−c_j)` 进入
  （Gnielinski 类的自然推广）。多特征拟合用分块交替最小二乘，每块都是
  条件精确最小二乘、SSE 单调不增、全程确定性；已知限制：收敛到局部
  平稳点而非全局最优（`params['alternating_ls_converged']` 报告收敛状态）。
- 单特征输入下两个家族都退化为原有数值路径，`params` 中保留
  `exponent` 标量键与 `predict_log10` 闭包以向后兼容。

降级与警告（真实数据验证驱动的两条规则）：

- **逐特征降级**：某 π 特征群含非正值（零是合法物理值，如攻角 0°）时，
  框架显式剔除该群、用剩余群继续拟合，并在结果 `dropped_features`
  与 `data_report['warnings']` 中留痕。刻意的不对称：源数据含
  NaN/inf 时仍整体失败——零是物理现实，NaN 是数据质量事故。
- **共享变量提示**：目标群与特征群共享变量（如自动群 p/(f·μ) 与含 f
  的特征群，或教科书情形 Nu 与 Re 共享 D）会在数学上抬高 R²
  （π 群伪相关）。结果 `shared_vars`（模板路径为 `warnings`）显式
  标注此类情形——这类 R² 不宜与无共享变量的拟合直接比较。

## 异常守卫（anomaly.py, v0.2.0）

套在任意 log 空间拟合外面的守卫层——先检测、后调整、全程留痕：

- **三层检测**：数据层（非有限值/非正值/零方差，致命即拒绝）；残差层（MAD 稳健 z 分数判离群，对最小二乘被带偏免疫）；几何层（hat 杠杆值 > 2p/n、标准化条件数 > 30 告警共线性）
- **五种调整策略**（`adjustment`）：`report_only`（只报告）、`robust`（Huber IRLS 降权）、`ridge`（CV 选择岭强度）、`exclude`（带审计剔除）、`auto`（按异常类型自动决策）
- **语义区分**：残差离群 = 疑似测量错误（可剔）；高杠杆 = 条件极端的设计特征（默认只报告，`exclude_leverage=True` 才连剔）

```python
from dimreflex import fit_with_guard

g = fit_with_guard(pi_data, y_col, x_cols, adjustment='exclude')
g['anomalies']   # 三层检测清单
g['before']      # 被带偏的基准
g['after']       # 调整后的候选
g['audit']       # 剔除审计: 索引/来源/得分
```

## 自主发现与自动后端调整（v0.3.0）

`discover_auto` 不再用训练集 R² 直接决定“最佳模型”。它会对注册函数族执行
确定性的 k-fold 验证，同时报告 AICc/BIC；验证性能落在一标准误差范围内的
模型优先选择复杂度更低者。每个候选都保留验证分数、折间结果、复杂度指标
和推荐标记。

```python
result = dr.discover_auto(
    data,
    target="Nu",
    features=["Re", "Pr"],
)
result["recommended"]       # 经验证的推荐候选及选择依据
result["assessments"]       # 所有函数族的 CV、AICc、BIC 与失败原因
result["backend_adjustment"]  # 自动调整决策与审计
```

量纲感知场景推荐使用 `discover_auto_pi`，它先生成并验证 π 群，再进入上述
自主选择链路：

```python
result = dr.discover_auto_pi(data, sys, target="h")
# 等价快捷入口:
result = dr.discover(data, sys, target="h", method="auto")
```

后端守卫新增 `adjustment="auto"`：

- 检测到残差离群：自动使用 Huber IRLS 降权；
- 仅检测到共线性：通过 k-fold 自动选择岭强度并稳定指数；
- 数据健康且无异常：保持原模型，不为“调整”而调整；
- 高杠杆点默认只报告，不自动删除；
- 所有决定写入 `applied_adjustment` 和 `audit.decision`。

```python
guarded = dr.fit_with_guard(pi_data, y_col, x_cols, adjustment="auto")
```

## 多场景拟合（scenarios.py, v0.4.0）

同一套设备数据经常混合不同机器、材料、站点或运行模式。`fit_scenarios`
不会预设“工况一定相同”或“每个工况必须完全独立”，而是同时比较三种结构：

- `pooled`：所有场景共享倍率和指数；
- `shared_exponents`：共享机理指数，但每个场景有独立倍率；
- `independent`：每个场景都有完整的独立幂律模型。

选择过程使用场景分层 k-fold 验证，保证每折都包含所有已知场景；进入一标准误差
范围后优先选择参数更少的结构，再比较 AICc/BIC，从而避免把随机噪声误判为场景差异。

```python
result = dr.fit_scenarios(
    data,
    target="Nu",
    features=["Re", "Pr"],
    scenario="machine_id",
)
model = result["final_model"]
print(model["strategy"])           # pooled/shared_exponents/independent
print(result["assessments"])      # 每种结构的 CV、AICc、BIC

prediction = dr.predict_scenarios(model, new_data)
```

训练时未见的场景不会被静默套用到旧模型上，而会显式拒绝预测；每个场景样本不足、
非正值、非数值或 NaN/inf 也会返回结构化失败原因。

## 工程异常持续监视（monitoring.py, v0.4.0）

`EngineeringMonitor` 把“坏数据”“超出适用范围”“工况人口改变”和“相关式失效”
分开监视，并保留最近批次历史：

- 数据有效性：缺列、空批次、非数值、NaN/inf、非正值、未知场景；
- 运行包络：逐特征超界比例，识别外推风险；
- 输入漂移：逐特征 PSI，以及可选的场景构成总变差；
- 模型漂移：稳健残差 z 分数、残差中心偏移、log-RMSE 恶化；
- 持续性：异常连续达到默认 3 批后升级为 `critical`，正常批次自动清零；
- 无目标预警：现场目标尚未回传时，可先仅监视输入和运行包络。

```python
built = dr.EngineeringMonitor.fit_baseline(
    baseline,
    target="Nu",
    features=["Re", "Pr"],
    scenario_column="machine_id",  # 单场景数据可省略
)
if not built["success"]:
    raise RuntimeError(built["reason"])

monitor = built["monitor"]
report = monitor.inspect(batch, batch_id="2026-09-01/shift-2")
print(report["status"])                  # normal / warning / critical
print(report["alerts"])                  # 代码、严重度、解释和证据
print(report["row_anomaly_positions"])   # 批次内异常行位置

# 目标值尚未回传时，先做输入侧监视
early = monitor.inspect(feature_batch, require_target=False)

# 查看/清除运行历史；校准基线和模型不会被清除
state = monitor.snapshot()
monitor.reset_history()
```

默认阈值可通过 `thresholds={...}` 在构建时覆盖。监视器只告警和留痕，不会自行
删除样本、重训模型或接管生产控制；这些动作应由工程审批流程决定。

## 测量不确定度传播（uncertainty.py, v0.5.0）

`Variable.uncertainty`（相对误差）不再只是声明。`propagate_uncertainty`
按声明对数据加乘性扰动（log 空间高斯），重新投影 π 群并重拟合，
把仪器精度传播成指数的置信区间：

```python
sys.add_variable("h", "W/(m**2*K)", "target", uncertainty=0.03)  # ±3%
res = dr.propagate_uncertainty(data, sys, target="h", n_sims=300)
res["exponents_mc"]    # 每个指数的 mean/std/95% CI 与人读摘要
res["regression_std_errors"]  # 对照:回归残差噪声,与本结果语义独立
```

两种噪声分开报告、不可相加：`std_errors` 是回归残差噪声（条件于测量
无误差），MC 传播的是测量噪声（条件于拟合形式正确）——后者反映仪器
极限，不随样本量消失。没有任何变量声明 uncertainty 时显式拒绝，
框架不在无测量噪声信息时虚构传播结果。

## 运行示例

```bash
# Python 分析脚本（数据检查 + 自动发现 + 模板拟合）
python validation/data_analysis.py

# Node.js 验证脚本（独立复刻，交叉验证结果）
node validation/analyze_ht.js
```

## 命令行 / 独立 exe

量纲感知模式：在 CSV 同目录放 `<数据文件名>.schema.json`（或用 `--schema`
显式指定）。schema 声明每列 unit/role/target，CLI 才会走 `System → π 群 →
discover` 的完整链路；`heat_exchanger_data.schema.json` 是根目录示例。

```bash
# 有 schema 时：自动进入量纲感知模式
dimreflex_app data.csv

# 显式指定 schema
dimreflex_app --schema units.json data.csv
```

无 schema 时是**诚实兜底模式**：仅做原始变量空间的 log-log 幂律拟合，并
在输出顶部明确声明"未校验单位/量纲、不代表物理量纲结论"。

```bash
# 无 schema：非量纲模式（仅数值曲线拟合）
dimreflex_app data.csv h rho v D mu cp k

# 省事模式：最后一列自动作目标，其余作特征（拖拽 CSV 到 exe 上等同此模式）
dimreflex_app data.csv

# 双击 exe：回车 = 内置双工况演示；或输入/拖入 CSV 路径
```

打包（Windows）：

```bash
pip install pyinstaller
py -3.13 -m PyInstaller --onefile --name dimreflex_app launcher.py
# 产物: dist/dimreflex_app.exe（自带 Python 运行时，无需安装即可运行）
```

## 测试

```bash
pip install pytest
pytest tests/
```

## 依赖

- numpy, pandas, scipy, pint
- （仅测试）pytest

## 代码质量

- 静态检查：`ruff check .`（规则集见 `pyproject.toml` `[tool.ruff]`，2026-09 清零）
- 测试：`pytest tests/`（144 项，含与验证脚本同源的正/负对照）

## 协议

MIT
