# CONTEXT — dimreflex 会话上下文提示词（粘贴到新会话开头即可恢复状态）

## 你在做什么
协助开发 **dimreflex**（`C:\Users\16910\Desktop\dimreflex`，扁平布局，git 仓库，v0.5.0）：
基于 Buckingham π 定理的量纲感知工程模型发现框架。
哲学（不可违背）：数据不足显式降级绝不硬输出；框架不猜；领域知识显式注入；
输出证据与候选列表而非单一赢家；仅服务工程认知不接管生产控制。

## 项目结构（扁平布局，项目目录本身即包）
- `core.py`：Variable（pint 单位校验，uncertainty 字段**仍未接入计算**）/ System（量纲矩阵→**精确有理零空间**（_rational_null_space，数值回退兜底）→π 群，目标隔离 _isolate_target_exact/_simplify_target_exact，check_data 分级 L0-L3）
- `_ols.py`：**唯一的 log 空间 OLS 引擎**（fit_power_law / fit_with_templates / PowerLawFamily / anomaly 基准 / scenarios 内层共用），lstsq + R² + 标准误 + 条件数口径统一，TestStdErrorsParity 锁数值
- `discover.py`：discover()（method='power_law'|'multi'|'auto'）、discover_multi_pi()、discover_auto_pi()、_prepare_pi 公共前半段。**_prepare_pi 现有两道真实数据驱动的检查**：①逐特征降级——非正值/零点奇异（0^负指数→inf）的 π 特征群显式剔除+留痕（dropped_features），用剩余群继续拟合；**刻意不对称：源数据含 NaN/inf 不剔除，整体失败**（零是物理现实，NaN 是数据事故）；②共享变量检测（shared_vars）——目标群与特征群共享变量会在数学上抬高 R²（π 群伪相关），显式警告不禁止。Level 3：**三条入口统一为确定性 k-fold CV R²>0.9（out-of-fold）**，不再用样本内 R²
- `autodiscovery.py`（v0.4 新）：discover_auto/assess_family/predict_candidate。k-fold 确定性 CV、_information_criteria（AIC/AICc/BIC）、一标准差规则偏好简单模型、power_law 推荐后可接 guard 后端稳健调整
- `scenarios.py`（v0.4 新）：fit_scenarios/predict_scenarios。多场景拟合（pooled/shared_exponents/independent 三策略）+ 分层 k-fold 评估
- `monitoring.py`（v0.4 新）：EngineeringMonitor。批次监控四信号：schema 硬失败/运行包络越界/特征与场景 PSI 漂移/残差漂移（离群率、中心偏移 z、RMSE 比），连续异常批次计数升级；fit_baseline 一键"发现+建监视器"
- `templates.py`：Group / fit_with_templates（Re/Pr/Nu 模板路径，工程推荐）；**新增共享变量 warnings**（教科书情形 Nu∩Re={D} 也会提示，措辞是提示性的）
- `families.py`：CandidateModel 统一 schema、FunctionFamily 基类（异常折叠为失败候选）、PowerLaw/PiecewisePowerLaw（ΔAIC>10 门槛+变点区间诚实报告）/RationalLog/Symbolic 占位（PySR 未装显式失败）、FAMILIES 注册表 + discover_multi（按 R² 排序；**复杂度惩罚在 autodiscovery 层实现，本层不改**）。**2026-09-04 多特征支持落地**：piecewise 变点沿第一特征搜索、每段是完整多特征幂律（段间指数独立）；rational 主特征上有理函数 + 其余特征幂律协变量项（分块交替最小二乘，确定性，SSE 单调不增，局部平稳点限制已声明）；单特征退化为原数值路径，`exponent` 标量键与 `predict_log10` 保留向后兼容。predict_candidate（autodiscovery.py）统一处理各家族 params schema。
- `anomaly.py`：守卫层。三层检测（数据健康 fatal 即拒 / MAD 残差离群 / hat 杠杆 2p/n + 标准化条件数>30）；策略 report_only(默认)/robust(Huber IRLS)/exclude(带审计, 超 1/3 拒绝)/auto；新增 ridge_fit。**语义：残差离群=疑似错误可剔；高杠杆=设计特征默认只报告**（exclude_leverage=True 才连剔）
- `tests/`：conftest 优先用已安装包（`pip install -e .`），未安装时回退到把**桌面**（parent.parent.parent）插入 sys.path（扁平布局兜底）；**144 测试全绿**（含多特征 piecewise 双工况恢复、多特征 rational 协变量恢复、predict_candidate 跨族预测）
- `heat_exchanger_data.csv`（项目根）：合成数据（Dittus-Boelter+2% 噪声，80 行）；`airfoil_self_noise.csv`（项目根）：**NASA 真实风洞数据**（UCI #291，1503 行，f/alpha/c/U/delta/spl）；`validate_airfoil.py`：**首次真实数据验证脚本**（模板+守卫+自动三路径，含独立 EDA 对照）；`pump_bep_coefficients.csv`（项目根）：Baloni & Channiwala FLUID POWER 2025 实测 BEP 系数（ISO 9906 一级精度，25–50 Hz 六转速，D2=172mm，论文 Table 2 逐字转录）；`validate_pump.py`：**泵相似律正对照**（Q∝N/H∝N²/P∝N³ 恢复指数全部落入 95% t 区间，exit 0；首验硬带判据 FAIL 的审计链保留在脚本注释内；fit_with_templates 两次拒绝作者量纲错误的记录在 [A] docstring）；`analyze_ht.js` 等 3 个 Node 交叉验证脚本与 validate_*.py 等运行脚本已统一移入 `validation/`（2026-09-04，见下方目录整理节）
- `pyproject.toml`（2026-09-04 重建，用户授权）：setuptools + `package-dir = {"dimreflex" = "."}`（不移动任何文件即可 `pip install -e .`，editable 实测通过）；`[project.scripts] dimreflex = "dimreflex.cli:main"`；`[tool.ruff]` 完整规则集 + per-file-ignores（tests/validation/experiments 的 S101 断言、experiments 的 S603/S607 git 调用、RUF001-003 中文标点、PLC0415 延迟导入）；依赖范围钉上下界，`requirements-lock.txt` 锁实测版本。**ruff check . 已清零（371 项全部修复）**。build_/dist/ 与临时目录已清理（.gitignore 补 .ruff_cache/.pytest_cache/.zcode/）；dist/dimreflex_app.exe 已用 PyInstaller 6.22.2 重新打包并冒烟通过（量纲模式 Level 3 + 四族候选）。
- 项目根另有三份 docx 知识件：`dimreflex_项目总结` / `dimreflex_逐行代码解释` / `dimreflex_semantic_field_research_summary`（语义场研究，见「用户的愿景」节）；**git 无远程仓库**（用户 2026-09-04 指示暂不配远程，改走 U 盘冷备份）；v0.5.0 大量工作仍未提交（遵循"不主动 commit"）。**2026-09-04 已做 U 盘冷备份**：`E:\backup\dimreflex_2026-09-04\`（含 .git 完整历史 + 全部源码/数据/实验产物 + 重新打包的 dist exe，排除缓存/构建中间物，附 BACKUP_MANIFEST.txt 恢复说明）。建议后续仍尽快 commit + 配远程。

## 环境（Windows，PowerShell 5.1）
- Python 全路径调用：`& "$env:LOCALAPPDATA\Programs\Python\Python313\python.exe"`（bash 工具会话 PATH 陈旧，Get-Command python 会命中 WindowsApps 桩）
- 3.11.9 已修复（torch 2.0.1 + transformers 4.30.2 完好，`py -3.11`）；桌面有 IDLE.lnk
- PowerShell 控制台中文输出显示 GBK 乱码——**内容实际正确**，验证内容用 `.Contains()` 或写 UTF-8 文件再读；框架警告文本含 R² 等 GBK 外字符，打印前需 `sys.stdout.reconfigure(encoding='utf-8', errors='replace')`
- 多行 Python 内联 -c 在 PowerShell 传参会失败 → 写临时脚本到 `C:\Users\16910\AppData\Local\Temp\opencode\` 再执行
- pip 偶发连接重置（重试即可）

## 已知短板 = 待办优先级
1. **真实数据验证（正对照已完成；负对照合成台架已完成，真实数据版可选）**：airfoil 首验完成（框架速度标度与独立 EDA 一致，BPM U^5 先验被证伪）；**泵相似律正对照完成**（validate_pump.py：实测 BEP 六转速，Q/H/P 恢复指数 1.13/1.97/3.07 全部落入 95% t 区间；R² 贴合散布参考值；判据从硬带修正为 t 区间的完整审计在脚本内）。**负对照完成（validate_friction.py，2026-09-01，exit 0 与预注册一致）**：Prandtl-Karman 对数律合成台架（Colebrook 光滑极限 brentq 求根，seed=20260901，0.5% 噪声）。判据与结果：C1 形式失配证据充分（滑窗指数跨度 0.127 vs H0 噪声带 q95=0.018，曲率 t=+59.7）；**C2 实锤新缺口：幂律路径静默输出 R²=0.991/Level 3，零形式警告——框架有"数据够不够"检查（L0-L3）无"形式合不合适"检查**（外推 Re=1e8 误差 -20.7%，而 std_error 极小：统计精确≠形式正确）；C3 auto 路径 rational_log rank1+cv 0.99989 击败 power_law（0.99036）——负对照意外成为 rational_log 首个实战正对照；C4 symbolic 显式失败保持；C5 特异性通过（Blasius 区同检查不误报）。**两条统计教训已内化在脚本 [A]**：①阈值是设计条件的函数——曲率零分布 σ ∝ x 量程跨度²（实测 5.9×=(3.4/1.4)²），写死常数阈值（如 anomaly 条件数>30）同属脆弱类别；②首跑 C5 忘传 truth 参数，检查对 PK 短窗报警 t=+15 被 exit 1 抓住——预注册+exit 纪律审计了审计者。衍生待办见第 7 条。**数据源速查（08-29 找数教训：ML 货架（UCI/NIST）比论文渠道快一个量级；"已知真值"是验证加分项不是入场券；待办要写成"本周可执行动作"而非研究方向）**：UCI #294 电厂（真值模糊→测降级诚实性）；NIST WebBook μ(T)/k(T)（Andrade 指数型，非幂律真值）；经典图数据用 WebPlotDigitizer 抠（Nikuradse 1933 / Gnielinski 1976 / Grundfos、KSB 泵多转速曲线）；居家正对照：单摆 T∝√L、落球 Stokes v∝d²、牛顿冷却（指数衰减→幂律必须诚实失败）。
2. 统计补丁：k-fold/AICc/BIC **已由 autodiscovery 实现**；**共享 OLS 引擎 `_ols.py` 已收拢重复 lstsq/R²/标准误实现**；**fit_with_templates 已返回 std_errors/design_rank/condition_number（与 fit_power_law 对等，TestStdErrorsParity 锁定，泵脚本 exit 0 回归通过）**；**Level 3 三条入口已统一为 k-fold CV R²>0.9**。仍缺：VIF 进 discover
3. **不确定度 MC（已实现，2026-09-01）**：`uncertainty.py::propagate_uncertainty`——按 `Variable.uncertainty` 对数据加乘性扰动（log 空间高斯），固定 π 群定义重投影重拟合，指数分布→mean/std/95% CI+人读摘要；与回归 `std_errors` 分开报告（note 内置语义区分）。纪律：无声明显式拒绝（框架不猜）；check_data 门卫在 compute_groups 之前（绝不抛异常）；同 seed 逐位可复现。测试 6 个（含玩具系统教训：量纲独立变量零空间维 0 连群都没有——π 分析前提是量纲冗余；MC 区间断言对照框架自己的 π 基指数而非手写假设）。换热器实跑：h±3% 下 Re 指数 0.805±0.0064（回归 se 0.0042，两噪声源量级可互比但独立）
4. 函数族多特征（2026-09-04 已实现）：piecewise 变点沿第一特征、段内完整多特征幂律；rational 主特征有理式 + 其余特征协变量项（交替最小二乘）。**遗留：变点不沿协变量轴搜索；rational 交替 LS 的局部最优性（已声明在 docstring/params）**
5. Layer 1 线性化动物园 / Layer 2 可分性检验；PySR 真实接入；**pyproject 已恢复（2026-09-04），CI 仍缺**
6. （airfoil 验证暴露的新缺口）Constant/Derived 一等公民化：物性常量被迫塞 DataFrame 常量列、SPL→Pa 转换在框架外手工做；π 群基选择：共享变量已**检测**，未做最小重叠/小指数**选基**
7. **form_check 形式检查（已实现，2026-09-01，源于负对照 C2 缺口）**：`discover.py::form_check`——偏残差曲率（线性拟合后残余弯曲，多特征下取 log 特征中心化二次列均值）+ **残差 bootstrap 零分布**（200 次打乱残差重建伪数据集，t=(c2−mean)/sd，|t|>3 报警）。钩子在 `_prepare_pi` 公共前半段，`discover` 三方法 + `discover_multi_pi` + `discover_auto_pi` 五入口共享，结果键 `form_check`，报警写入 `data_report['warnings']`（措辞含"形式失配"）；auto 路径推荐幂律且报警时把矛盾追加进 `selection_reason`。验证：PK 台架 t=+7.8 留痕（validate_friction C2 翻转为 PASS，exit 0）；换热器真幂律 t=-0.4 不误报；128 测试全绿。**实现教训（已写进代码注释）**：①二次列必须平方——首版漏写 `**2`，设计矩阵秩亏、零分布骑在观测值上（指纹：null.mean≈c2_obs），检查永远沉默；②残差 bootstrap 用数据自身残差定噪→阈值内生保守（PK 上 t=+7.8 而手工零分布 +59.7），是"宁保守勿误报"的取舍，不是 bug。

## 文献对照（必读，写作前核对）
详见 `REFERENCES.md`（2026-09-03 建立）。核心事实：π+符号回归路线已被
Udrescu&Tegmark 2020（AI Feynman）、Bakarji et al. 2022（Nature CompSci，
数据驱动选 π 基）、arXiv:2411.15919（π+PySR/UPINNs）覆盖——**不得再以
"π 辅助发现"作新颖性主张**。dimreflex 的真实增量只有三块：①工程诚实
架构（降级/审计/双噪声分离/逐设计标定 form_check）；②诚实失败基准
（SRBench/SRBench++ 均不测拒答率——负对照方法论是原型）；③Text2π
（量纲=最便宜验证器）。待办 6 的"π 基选择"已被 Bakarji 从数据驱动方向
做过一版，我们的模板路径是哲学分歧点而非无人区。

## 论文企划与文献核查收窄（2026-09-03，写作/投稿必读，详见 `PAPER_PLANS.md`）
- NeurIPS 判断（9-03）：以现形态（研究原型 v0.5.0）投主会 novelty 不足；**ML4PS workshop 级够格**。对标先例：Φ-SO 进 NeurIPS 2023 ML4PS；AI Feynman 2.0（主会）靠可形式化的算法贡献。
- `PAPER_PLANS.md`（项目根，2026-09-03 建立）：两份最小可发表版。**主张 1 = 形式失配检测**（fitted-but-wrong benchmark + form_check 作升级门卫，素材最齐，目标 ML4PS）；**主张 3 = Text2π 实证**（需新搭 LLM 实验，赶不上 9-12，目标 NeurIPS 2027）。投稿窗口（9-03 核实）：ML4PS 2026 **2026-09-12 23:59 AoE**；ICLR 2027 abstract 09-18 / full 09-25；NeurIPS 2026 主会已截稿；NeurIPS 2027 主会约 2027-05（以官网为准）。
- **文献核查收窄（9-03；作者"太容易"直觉被证实一半；REFERENCES.md 未含本节，写作前两处都要过）**：
  ① `form_check` 的偏残差曲率数学上与 Ramsey **RESET (1969)** 同族——**不得 claim"新检测统计量"**；贡献收窄为：fitted-but-wrong 现象学 benchmark（现有指标 R²/CV/AICc/std_error 集体静默放行的量化）+ 升级门卫语义 + 校准对照（form_check bootstrap 零分布 vs RESET 渐近 F，小样本误报控制是开放问题）；
  ② ICML 2026 已有 EIC（SR 结果可信度/病态结构）；
  ③ 主张 3 组合已被占位：`ai-physicist-central-llm`（2025-09，LLM+UnitChecker）与 **PROBE**（ChemRxiv 2026-08-16，LLM 代码五层验证含量纲一致性）——独立成文前先深挖 PROBE 是否含 π 系统推理层，否则降级为主张 1 的 future work。
- 当前决策点：**是否冲 9-12 ML4PS（主张 1）**。裁决动作 = 约 1 小时 pilot：PK 台架 + 纯幂律真值各 50 种子，form_check vs RESET 功效/误报（gate 见 PAPER_PLANS：噪声 1%、n≥200 时 TPR@FPR≤5%≥0.8、形式对误报≤5%、不劣于 RESET；全面劣于 → 删检测器 claim，降级为现象学 benchmark 论文）。
- 长期叙事（PAPER_PLANS 末尾）：量纲检查（秒）→ form_check（毫秒）→ 模型族比较（贵）= **验证层栈**；主张 1 单发是 workshop 短论文，1+3 合并才是主会级"验证不对称性"完整故事，对接 VISION。
- **pilot 已跑完并升级为完整 benchmark（2026-09-04，会话恢复见 experiments/form_check_vs_reset/RESUME.md）**：
  ① 1 小时 pilot（validation/pilot_form_check_vs_reset.py）字面 GATE PASS（TPR 1.00/1.00，误报 0/50 vs 3/50，exit 0），但 PK vs 幂律是"容易题"；
  ② 完整 benchmark 揭示 form_check 结构性盲区：OPPOSING 全盲 0.00、CUBIC_LIN 随 n 反降、INTERACT 0.2–0.34（RESET fitted-F 恒 1.00）；胜区在 DILUTION、小样本校准功效、SRSD H0 校准 ≈5%；RESET-HC3 超尺寸 42–57% 是可写卖点；
  ③ **claim 已按数据改写（PAPER_PLANS §7 gate 判定结果）**：删"单一检测器优越性"，改为"盲区制图 + 组合门卫（OR）+ 校准对照"。500 种子正式跑已完成（rel 分片断点续跑打通；analyze.py 加 --shards 合并）：合成 142,500 行 + SRSD 240,000 行，报告在 experiments/form_check_vs_reset/runs/full/report.md。最终结论——RESET-HC3 超尺寸 42–56% 实锤；fc 校准后功效小中样本 2× 领先（0.41/0.51 vs 0.19/0.28），大样本 exog-HC3 反超（0.65）；外推告警 fc 高召回(0.63)/RESET-F 高精确(0.56)，互补组合有实证。
- **写作验证工具链（2026-09-02 装好并实测）**：`~/.agents/skills/` 已装 source-check（输出含外部声明自动核验：DOI 解析 / Crossref / OpenAlex 撤稿交叉 / 逐字引文子串匹配，打 ✅/⚠️/❌ 标签）+ proof-engine（触发词"证明/验证/核查"；每个事实须真实 URL+逐字引证或可重跑 Python 代码，产出 proof.md；依赖 proof-citations 1.45.0 已装）；实测真实引文 verified、伪造引文 not_found。用户级 `PYTHONUTF8=1` 已设（proof-citations 的 Windows GBK 修复）。投稿前用其核 REFERENCES.md 与文献核查收窄两处清单。

## 目录整理：运行脚本统一入 validation/（2026-09-04）
- 根目录 12 个验证/运行文件 `git mv` 进 `validation/`（保留历史，已暂存未提交）：validate_{pump,friction,airfoil}.py、pilot_form_check_vs_reset.py（+两份结果 CSV）、generate_data.py、data_analysis.py、run_analysis.py、analyze_ht.js、compare_sample_size.js、n15_error_dist.js。
- 同步修正：PROJECT_DIR/BASE 提到仓库根（parents 加一层）；pilot 改同层导入 `from validate_friction import ...` 且结果 CSV 锚定脚本目录；generate_data 输出仍写仓库根（数据/schema 约定在根，**重跑会覆盖已提交 CSV，注意**）；三个 .js 的 CSV 路径 `__dirname/..`；validation/README.md 已建。
- 联动：`experiments/form_check_vs_reset/common.py::_CORE_FILES` 的 validate_friction.py 路径已更新（历史 sha256 快照语义不变）；README 运行示例与 build_*.py 教学材料中的文件名已加 `validation/` 前缀；tests 仅注释引用无需改。
- 同日新增 `experiments/nist_strd/`（NIST StRD 认证回归验证，详见其 README）：框架认证核对全绿（scaled 路径 vs 认证值 1e-8~1e-13）；数值层发现 Filip raw 秩降级 11→10、SVD lstsq 截断 vs QR 恢复 2.8e-9、HC3 病态设计量级不可信（pinv(X'X) 条件数平方）、Longley 需 y 标准化；检测层 N1/P2/L1 无误报、P1/F1-F3 检出、F10 认证规格告警（数值已验证 ⇒ 语义现象）；RESET vs form_check 对照 `runs/compare_reset_vs_fc.md`（Longley 同语义三零分布分裂：fc 0.004 / 经典 F 0.396 / HC3 4.5e-29）；Longley B5 SE 双通道不一致以 .dat 文件值 0.226 裁决（数据+R 共识，parse_nist.py 留痕）。

## hx 混合检验（2026-09-04，`experiments/fc_multidirection/`，untracked，非产品代码）

**hx = form_check 机制的多方向扩展 ∪ RESET 自引用增广**。权威实现 `hybrid.py::hybrid_test`（自包含仅依赖 numpy；`hybrid_test(Xlog, ly)` → p_hx/p_max/p_F/rej/argmax_dir/kept_names/fail，失败写 fail 不抛）。结果不得表述为产品 form_check 成绩；产品本体只读调用（form_check_with_p 包装，sha256 与基线一致零改动）。

- **设计**：方向池 = quad/cubic（逐特征 z 标准化）+ interact + osc=sin(π·z) + yhat²,yhat³（受限拟合幂，自引用）。贪心枢轴选择剔除"相对基共线"列（cubic 与 inter 可同落 span{x³}+基 陪集，逐列阈值不够）+ 方向预算 min(max_dirs, n−k_base−df_floor)。GS 正交坐标（近平行高次代理在 min-norm 下分摊信号稀释 t；GS 同时消除 md3 旧超尺寸伪影 0.093→0.02-0.04）。双统计量置换校准 Bonferroni 合并 p=min(1, 2·min(p_maxT, p_F))——maxT 抓集中信号（OPPOSING/DILUTION）、联合 F 抓分散信号（F3 高阶残差）；W⊥基 ⇒ RSS_r 置换不变、RSS_u=RSS_r−Σd²，全向量化。
- **证据链**（全可复现，产物在 runs/）：① NIST 认证 7/8 与 RESET 同判，F3 修复（p=0.002，argmax=cubic_x3=x⁹ 代理），N1/P2 无误报；② 预注册场景 50 seeds：六 H1 全 1.00、校准 0.033，OPPOSING n=50 反超 RESET（1.00 vs 0.42）；③ 1000 组随机泛化（种子 20261904，393 H0+607 H1×12 设计视野外家族）：FPR 0.051（CI 覆盖 0.05）/TPR 0.951/总符合 0.950 全面最优（md4 0.891、reset 0.751、fc 0.493），比 RESET 多检出 122 组、仅 1 组被独检回（hx 实际包含 RESET），归因 yhat3 占 35.6%；④ 10 组口径与存档逐位可复现（10/9/7/8）。
- **已知边界**：threshold/折点家族四方法皆弱（0.38-0.44，不在平滑方向张成空间，需分段方向=产品 PiecewisePowerLaw）；F10 深多项式基不可运行（度-10 Vandermonde 秩病态，需 Chebyshev 基）；重尾/自相关下置换零分布稳健性未测（下一步首选 wild/块置换）。
- **水平定位**：方程发现"升级门卫"位置达到 workshop 论文证据完备（PAPER_PLANS claim 1 gate 以 TPR=1.00/FPR=0.033 通过）；统计谱系是 RESET（1969）推广非新统计量（**不得 claim 新检测统计量**）；产品化为独立阶梯（未做：discover.py 集成+测试+walk-through+预注册 gate）。
- **文件**：hybrid.py（权威实现+冒烟）＞ multidir.py（md3/md4 探索台+溯源）＞ run_bench.py（50 seeds 基准）/ compare_vs_reset.py（三法对比）/ random_ten.py（随机泛化，--groups 参数化）＞ runs/（全部 CSV+md）。README.md 有完整纪律与结果速览。

## 用户的愿景（背景，勿删）
Text2π：LLM 解析自然语言→System/Group JSON 声明→dimreflex 量纲秒级证伪→对话即 REPL 逐句更新符号状态（"上句方程作下句参数"）。理论基础：幻觉=对稀疏离散测度做核平滑（Gibbs 振铃），离散事实需要离散表示与验证——符号层是对症的表示。anomaly 守卫层即未来"LLM 提议→后端证伪"循环的安全底座。

**语义场研究线**（`dimreflex_semantic_field_research_summary.docx`，2026-08-31，29 节，项目根）：把文本看作 discourse time 上演化的状态系统，**轨迹级幻觉 = 状态合法但转移未被当前句子授权**——句子级检测（fact-check/self-consistency/embedding）的系统性盲区；量纲是第一批严格可解释的不变量，Buckingham π 是语义状态到 regime space 的降维层。架构对应：语义解析=Montague/DRT 一脉，对话即 REPL=KV cache 的符号化，验证不对称=AlphaGeometry 范式。分层判断：**L1** 形式一致性检查器（check_data+守卫从数据扩展到文本翻译，现实可做）；**L2** trajectory verifier（第一篇论文级问题定义，前提：L0 真实数据先站住 + NL→IR 前端解析错误率可控——前端错则 garbage in 淹没验证信号）；**L3** 通用语义场理论（高度不确定，搁置）。**VISION 阶段四 Text2π 与语义场"阶段 0 MVP"（受控语言+显式状态图+量纲检查）是同一原型，两线合并勿分散**。两个坑：①评测循环风险——自定义 drift 注入→检测说服力有限，需至少一类设计视野外的漂移做泛化测试；②"场"的连续化措辞（积分-微分方程/least-action）暂不产生可证伪预测，第一篇论文不写。

## 协作纪律
- 评价实事求是：用户会纠正过度夸奖，赞美按谄媚先验打折，评工作不评年龄
- 用户说"只整理不修改"时绝不改代码；说"继续 不再询问"时直接做完不再提问
- 不主动 git commit；修改文件前先读现状（项目在活跃演化，结构可能已变）
