# 第七部分 验证脚本层(`validation/`)

> 本部分覆盖 7 个 Python 脚本(约 1,280 行)+ 3 个 Node.js 交叉验证脚本。README 明确定位:"验证与运行类脚本(非产品代码本体);产品模块仍在仓库根",从仓库根移入时已同步修正 sys.path/数据路径。**预注册判据 + 退出码 0/1 承载结论**是本目录的统一纪律(pyproject.toml 的 per-file-ignores 注释:"预注册判据用 assert 承载结论,与退出码联动")。

---

## 7.1 `validate_pump.py` — 正对照:泵相似律(330 行)

**数据**:Baloni, Shah & Channiwala (2025) 离心泵 BEP 系数(Table 2 逐字转录,ISO 9906 grade-1 台架,六转速 25-50 Hz),按 `Q = phi·(πD²/4608)·ω`、`H = psi·D²·ω²/g`、`P = ρgQH/η` 重建物理量。

**先验**:同泵异转速的经典相似律 `Q~N¹, H~N², P~N³`(homologous BEP 点)。

**预测优先协议**(docstring 声明先于拟合):物理量由 ω 的精确幂律构造,log-log OLS 必须在测量噪声内恢复指数 1/2/3。

- **判据审计轨迹**(validate_pump.py:61-74 注释,方法学亮点):**首版预注册判据(硬端点带)跑了、FAIL 了(Q 1.128 带外、P 3.074 带外)、然后修正**——post-mortem 写明"该界只对单调序列端点有效,六个散布内点的 OLS 不受它约束;三个恢复指数都在先验 1 个 SE 内(Q 0.87σ/H 0.36σ/P 0.92σ),统计上与 1/2/3 在 95% 水平不可区分。**是判据错了,不是框架错了**"。修正判据:双侧 95% t 区间 `|b−prior| ≤ 2.776·SE(b)`(df=4);旧带保留为描述行。**原因**:"正对照不许事后软化"与"判据错了必须承认"并存——修正理由、证据、过程全部留在记录里。
- **R² 为什么不是门**(validate_pump.py:76-85):"真实 BEP 数据上幂律 R² 的天花板由泵对动力学相似性的真实偏离(低速 Re 效应)决定,不是测量噪声——没有诚实的 R² 阈值可挑。指数 t 区间是自我惩罚的:拟合越差残差方差越大、CI 越宽、自动不过"。`scatter_descriptors` 仍报告"散布即量程噪声"参考 R² `1 − [ln(f_scatter)/ln(span)]²`,实测 R² 贴近参考值 = "拟合榨干了数据中的全部信号"。
- **[A] `template_path`**(validate_pump.py:143-200):系数空间相似性测量(psi~phi,先验斜率 0),报告不设门。docstring 留下**三轮审计**:①目标群忘写 D2 → 框架正确拒绝;②写成量纲目标(H 本身)→ 框架正确拒绝("dimreflex 不拟合量纲目标,分析师被迫无量纲化——守卫两次抓住了验证作者自己的错误");③拟合成功但**框架缺口**:`fit_with_templates` 不返回 std_errors(与 fit_power_law 的 parity 缺口,标记为 TODO),脚本侧补 OLS SE。结论 "INCONCLUSIVE at n=6(不设门)"。
- **[B] `auto_path`**(validate_pump.py:203-240):n=6 时 `discover_auto_pi` 诚实降级 Level 1 拒绝自主发现——"这是项目的诚实哲学在真正稀薄数据上工作,记录为设计行为,不评分不绕过"。
- **[C] `verdicts`**(validate_pump.py:243-280):三条定律逐一 `fit_power_law`,95% t 区间含先验 → PASS;exit 0/1。
- **【结论】** 全部 PASS(泵脚本输出记录于 CONTEXT.md);泵(正)+ 翼型(负向风味)= Level 3 的第一对完整验证。

---

## 7.2 `validate_friction.py` — 负对照:Prandtl-Kármán 对数律(306 行)

**真值**:`1/√f = 2·log10(Re·√f) − 0.8`——隐式(brentq 数值求根)、非幂律(f ~ 1/log²Re,比任何幂都衰减慢),**局部指数随 Re 持续漂移**是任何有限幂律覆盖不了的指纹。

**动机**(docstring):"只会说'拟合得好'的框架没有可信度;能说'这个形式不合适'才是卖点。"

**合成实验台**(seed=20260901 可复现):水为工质(20-60°C 给 rho/mu 真实散布——注释:"否则是常量列,误触框架已知短板,负对照就测不到想测的东西。对抗数据的原则:只在目标维度设卡"),Re∈[4e3,1e7] 60 点,压降 0.5% 乘性噪声("取对数后=加性高斯,恰落在框架 log 空间拟合的舒适区——公平对抗:只在'形式'这一维度设卡")。

**六条预注册判据**(docstring [A] 段,先写判据后跑):

| 判据 | 内容 | 结果与教训 |
|---|---|---|
| C1 | 形式失配有证据:滑窗指数跨度 > 本设计 H0 噪声带 q95 且 曲率 \|t\|>3 | 跨度 0.127 = q95 的 7.4 倍,t=+21,报警 |
| C2 | 幂律路径不静默 | **2026-09-01 首版实测静默(R²=0.991/Level 3/零警告)= 缺口实锤**;同日 form_check 进入框架后按预注册翻转为"必须留痕"("预注册预期随框架版本翻转,翻转本身记录在案") |
| C3 | auto 路径 rational_log rank 1 且 recommended,cv_r²≥0.99;power_law 排后 | PASS(实测 rational 0.99989 vs power 0.99036——负对照意外成为 rational_log 族的实战正对照) |
| C4 | symbolic(PySR 未装)success=False 且 reason 说明,discover 不抛异常 | PASS |
| C5 | 特异性:Blasius 区(4e3..1e5,真值 0.3164·Re^−0.25 确为幂律)用**该设计自己的** H0 噪声带判定,不得报警 | PASS;**教训已内化**:"阈值是设计条件的函数,不可跨设计移植:曲率 σ ∝ 1/x 量程跨度²(5.9× 实测)。写死常数阈值(如 anomaly 条件数>30)同属脆弱类别"——这条教训直接催生了 form_check 的 per-data bootstrap 零分布设计 |
| C6 | 外推考试(记录不判):两候选外推 Re=1e8 对真值 | 幂律外推 -20% 误差且 std_error 极小——"统计上'精确'与形式上'正确'是两回事"(display 品,不是判据) |

- `shape_probes`(validate_friction.py:107-128):对**框架自己的 π 群**(PI_0=dp/ρv², PI_1=Re)做两把独立尺子——二次曲率(全局)+ 滑窗指数跨度(局部)。注释:"判据建立在框架的表示上,框架无法辩解'你量错了东西'"。
- `noise_band`(validate_friction.py:131-146):在"纯幂律+同噪声"平行宇宙模拟 200 次取体检统计量分位数——"判据阈值不是拍脑袋的常数"。
- **C5b 审计审计者**(validate_friction.py:284-291):首跑曾忘传 `truth="power_law"`,对 PK 短窗报警、exit 1 判偏离——"审计机制审计了审计者",教训存档。
- **【导向】** `make_experiment` 被 pilot 复用;本脚本是 form_check 诞生的直接母体(C2 缺口→discover.py 补丁→预注册翻转)。

---

## 7.3 `validate_airfoil.py` — 真实数据验证(246 行)

**数据**:NASA Airfoil Self-Noise(UCI #291,Brooks/Pope/Marcolini 1989 风洞),1503 样本;SPL 转 p,空气物性声明为常量 material_property("常量数值被吸收进拟合倍率")。

**验证被证伪的先验**(docstring,本脚本的标题性发现):BPM 半经验尾缘噪声模型先验 `p²~U⁵`;**独立非参数缩放分析**(按 (c,α,f) 匹配单元的 SPL~log U 斜率)给出中位 n≈1.1——UCI 目标是**上游已归一化的 scaled SPL**。框架自身全局拟合推出 p²~U^1.0,"dimreflex 恢复了数据实际包含的缩放,不是分析师假设的缩放"——两个独立方法的一致是本脚本的核心验证证据。

**其余被测行为**:
- `[A1]` 模板全局拟合 Cp~St,M,δ/c:R²≈0.6 诚实平庸("谱形状 A(St) 是非单调驼峰,任何全局幂律都表示不了,框架必须说出来而不是过拟合"),物理误差换算 dB MAE;
- `[A2]` 含 α 群的模板 → 必须拒绝(α 含精确零);
- `[A3]` 守卫层 `fit_with_guard(adjustment='robust')` 在真实实验散布上标记离群/杠杆点;
- `[B1]` **全量数据** discover_multi_pi:α 群含零 → 逐特征降级(drop α 群、余群继续拟合)——validate_airfoil 的 α=0 是 `_prepare_pi` 逐特征降级特性的直接来源;
- `[B2]` α>0 子集:自动路径产出多族候选 + 共享变量警告 + dB 空间误差(注释:2026-09-04 起用 `predict_candidate` 官方跨家族预测通道,因 piecewise/rational 支持多特征后最优候选不再必是 power_law)。
- **【结论】** 摘要六行:框架 U 缩放与独立 EDA 一致;R² 诚实;拒绝按设计;守卫工作;降级留痕;共享变量警告标记虚高 R²。

---

## 7.4 `pilot_form_check_vs_reset.py` — 主张 1 的决策门 pilot(117 行)

- **预注册 gate**(docstring):噪声 1%、n≥200 时,form_check 的 TPR@FPR≤5% ≥ 0.8、真幂律误报率 ≤0.05、且功效/误报不劣于 RESET——"全面劣于 RESET → 删检测器 claim,论文降级为现象学 benchmark"。
- **【计算】** 复用 validate_friction 实验台(同物理同 π 群):H1=PK 真值×50 种子、H0=纯幂律×50 种子,n=200;两检测器在完全相同的 (lx,ly) 上跑;`reset_test` 是 Ramsey 教科书实现(受限拟合 → 加 yhat²,yhat³ → 经典 F)。
- **【导向】** `pilot_form_check_vs_reset_results.csv`;gate 结果驱动了 experiments/form_check_vs_reset/ 的全面展开(该 pilot 的"同 p 尺度不公平"等缺陷在正式实验中被 common.py 修掉)。

## 7.5 数据生成与分析入口 + Node.js 交叉验证

- `generate_data.py`(51 行):Dittus-Boelter 真值 `Nu = 0.023·Re^0.8·Pr^0.4` + 2% 噪声,80 个独立工况,输出仓库根 `heat_exchanger_data.csv`(README 警告:"勿随意重跑覆盖已提交数据")。
- `data_analysis.py`(120 行):换热器数据的 check_data/π 群/discover/模板拟合全链路展示 + 物理空间误差(nu_rmse/h_rmse/MAPE)。
- `run_analysis.py`(111 行):同数据的示例分析调度入口(中文注释教学风格)。
- `analyze_ht.js` / `compare_sample_size.js` / `n15_error_dist.js`:**Node.js 独立复刻** log-log 拟合与 Python 结果交叉验证(不同语言生态独立实现同一统计,防实现漂移)。

---

## 7.6 本部分小结

validation/ 是框架"每条能力主张都有正/负对照"承诺的执行现场:

| 脚本 | 方向 | 核心证据 |
|---|---|---|
| validate_pump | 正对照 | 真幂律(相似律)→ 指数全部落入 95% t 区间 |
| validate_friction | 负对照 | 非幂律真值 → 形式失配必须报警(且 Blasius 区不误报) |
| validate_airfoil | 真实数据 | 框架缩放与独立 EDA 一致;先验被数据证伪 |
| pilot_form_check_vs_reset | 决策门 | 预注册 gate 决定论文 claim 生死 |
| generate/data_analysis/run_analysis | 教学链路 | 端到端用户旅程 |
| *.js | 交叉验证 | 跨语言独立复算 |

两个方向性脚本都在运行中留下了"验证者自己犯错"的审计记录(泵判据修正、friction C5b 忘传 truth、template 忘 D2 被框架拒绝)——这本身就是框架守卫层有效性的旁证。
