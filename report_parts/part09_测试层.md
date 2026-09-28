# 第九部分 测试层(`tests/`,10 个文件,约 1,935 行,116 个测试)

> 单元测试统一放 tests/(pyproject 的 `testpaths = ["tests"]`,避免把实验归档中的历史样例误当正式测试)。`conftest.py` 提供 import 兜底:优先已安装的 editable 包,未安装时把项目父目录插入 sys.path(扁平布局要求)。测试的总体特征:**每个测试类的 docstring 都在讲"为什么"**(对应哪次事故/哪个待办/哪条哲学),断言对象刻意选在"框架自己的表示"上。

---

## 9.1 `test_core.py` — 量纲引擎(120 行,10 测试)

- `TestVariable`:速度/密度的 7 维量纲向量逐位断言;非法单位/非法角色抛 ValueError。
- `TestSystem`:7 变量系统秩 4 → π 群数 3(Buckingham: n−r);每个群 `D@vec=0`;目标隔离后恰好一群含 h 且指数=1、群定义精确等于 `{"k":-1,"D":1,"h":1}`。
- **`test_canonicalization_never_breaks_dimensionless_constraint`**(test_core.py:64-88,回归测试的典范):一个特制的 4×7 量纲矩阵,注释直说"旧的逐分量取整会让这个矩阵产生约 0.049 的量纲残差"——测试锁住精确有理零空间重构的修复,断言群指数全为整数且 D@vec=0(atol=1e-9)。

## 9.2 `test_discover.py` — 发现层集成(333 行,14 测试)

- `test_template_recovers_dittus_boelter`:合成数据(真值 0.023·Re^0.8·Pr^0.4)经模板路径恢复——prefactor rel=10%、指数 abs=0.05 容差。
- **`TestNonFiniteGuard`**(4 测试):NaN/inf 在 fit_power_law/templates/discover 端到端全部显式失败不崩。注释记录历史:"fit_power_law 曾因 NaN 直接抛 LinAlgError 崩溃"。
- `TestPiSpaceMultiFamily`(4 测试):method='multi' 与 discover_multi_pi 在 π 空间跑全族、目标隔离保持、Level 3 升级、缺列降级。
- `TestPerFeatureDegradation`(3 测试):**零点只剔除受影响群并留痕、NaN 维持整体失败**——"不对称设计"逐条断言(dropped_features 为空 + fit 失败)。注释:"空气动力学真实数据验证驱动的架构改进"(即 validate_airfoil 的 α=0)。
- `TestLevel3UnifiedGating`(3 测试):CV≥0.9 升级且警告含 "k-fold";样本内 0.99 但 CV 0.85 → 不升级 + "未升级"警告;CV 不可用 → 不升级且无警告。
- **`TestStdErrorsParity`**(2 测试):锁住 templates 与 power_law 的 std_errors/exponents 逐位一致(rtol=1e-9)。注释直白:"泵相似律正对照时模板路径缺 std_errors,脚本被迫手算兜底(待办 2)。本类锁住补丁"。
- `TestSharedVariableDetection`(3 测试):tiny 系统(特征群 b/a 与目标群 u·c/a 共享 a)与教科书 Nu~Re,Pr 情形都应提示;不相交时无警告。

## 9.3 `test_form_check.py` — 形式检查(119 行,7 测试)

- 三件事钉死(docstring):真幂律不误报(C5 的单元化)、对数律真值报警(C1 的单元化)、全链路留痕。`make_log_law_data` 用 brentq 现场解 PK 隐式方程(与 validate_friction 同源)。
- `test_structured_rejection_nonpositive`/`test_too_few_samples_skipped`:非正值与 n=8 → checked=False + 原因分类;`test_deterministic`:同 seed 逐位一致。
- `TestDiscoverIntegration`:`test_power_law_path_flags_log_law` 断言**拟合本身依旧漂亮(R²>0.98)但 form_check 报警且 warnings 含"形式失配"**——这正是"统计精确≠形式正确"的单元化;`test_power_law_data_no_false_alarm` 保证不破坏干净输出。

## 9.4 `test_anomaly.py` — 守卫层(~200 行,19 测试)

- `TestDataHealth`(5):五类数据层问题各就各位;`TestDetection`(5):干净数据零标记+audit note "无需调整";注入 20× 离群点被 MAD z 抓住;重复特征触发共线性;n≤p+1 杠杆不可辨识。
- `TestAdjustments`(8):**robust 比 OLS 更接近真值**(注入点后 OLS 指数偏移、Huber 恢复);exclude 剔除注入点且 audit 记录索引/来源;**`test_exclude_refuses_mass_flagging`**(标记超 1/3 → audit['refused'] 拒绝文本);report_only 永不调整;`test_exclude_leverage_semantics`(默认不剔杠杆点/exclude_leverage=True 才剔);robust 不被"仅杠杆"阻塞;未知策略显式失败;fatal 数据层拒绝拟合。
- `TestIRLS`:收敛 + 降权计数;`TestPiSpaceIntegration`:守卫包装 discover 流水线。

## 9.5 `test_families.py`(30 测试)、`test_autodiscovery.py`(10)、`test_scenarios.py`(9)、`test_monitoring.py`(11)、`test_uncertainty.py`(6)、`test_cli.py`(16)

- **test_families**:基底类契约(子类抛异常被兜底为失败候选、无原因失败补默认)、CandidateModel to/from_dict 往返与向前兼容、幂律族与 fit_power_law 行为一致、分段族(单段如实/变点检出/增益不足退单段/断点区间)、rational 族(单特征 MVP 路径/多特征交替/发散显式失败)、symbolic 占位、注册表与 discover_multi 排序。
- **test_autodiscovery**:assess_family 的 CV 契约(确定性折、折失败留痕、"验证未完成"与低分分离)、discover_auto 的推荐(一标准误差内偏简单、family 优先序)、后端调整(auto→幂律推荐时跑 fit_with_guard、非幂律记 not_applicable)。
- **test_scenarios**:三结构拟合与预测、共享指数结构恢复、未见场景 ValueError、场景分层 CV、n_params/AICc 选择、样本不足拒绝、标签碰撞拒绝。
- **test_monitoring**:基线标定、数据校验失败报告、包络越界/PSI 漂移告警、残差离群率/中心漂移/RMSE 恶化、未知场景失败、场景构成漂移、持续异常升级(连续≥3)、阈值未知键报错、历史有界与 snapshot/reset_history 语义。
- **test_uncertainty**(6):**docstring 记录玩具系统教训**——初版用量纲独立的 x[m]/y[s],零空间维 0 连群都不存在("玩具系统先手算零空间维数再写测试"),改成三个长度变量(零空间维 2)。`test_basic_propagation` 的注释是方法论金句:"指数的真值是框架表示的函数:π 基选出目标群 y/x、特征群 x/L,真指数是 −0.27 而非 0.73——断言对照基准拟合值(框架自己的表示),不对照手写假设(与负对照'判据建在框架表示上'同款纪律)"。其余:无声明拒绝/非法值拒绝/缺列门卫/确定性/回归 se 对照与 note 语义。
- **test_cli**(16):退出码契约(成功 0/降级 1/用法 2)、schema 副文件自动发现、schema 各字段校验错误、target 缺失或多于一个、演示模式、非量纲模式的诚实声明、引号剥离。

---

## 9.6 测试层小结

三个反复出现的断言模式构成这个测试套件的签名:

1. **"显式失败"模式**:`assert not result["success"]` + `assert "原因关键词" in result["reason"]`——几乎每个公开 API 都有一条,锁住"绝不抛异常/绝不硬输出"。
2. **"判据建在框架表示上"模式**:断言的真值目标取框架 π 基的等价形式(−0.27 而非 0.73;form_check t 而非手工曲率),避免测试者与框架各说各话。
3. **"回归测试携带事故档案"模式**:每个测试类的 docstring 指向一次真实事故(NaN 崩溃、canonicalization 残差 0.049、parity 缺口、零空间维 0 玩具、v0.4 静默缺口),修复即测试,测试即档案。
