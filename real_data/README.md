# real_data — form_check 真实工程数据台架

2026-09-06 建立，同日二期、三期扩容。为 form_check vs RESET 比较提供**大量真实工程数据**：
九个免鉴权开放源，91 个数据集，150 条候选幂律关系，761,515 行原始观测。
全部原始字节落盘(`raw/`，约 115MB)、清洗脚本可重放(`prepare_data.py`)、
来源与 sha256 记录在 `clean/manifest.json`。

## 与已有真实数据的分工

`runs/real_engineering_holdout/`(NIST Chwirut1/Nelson)与 `wb_run.py`
(世界银行宏观)是已有的两条真实数据线;本目录把它们没有的**物理工程
谱系**(流体/土木/水利/能源/气象/地震/化工/风工程/航空/汽车)补齐成批。

## 数据源与清单(一期 + 二期 + 三期)

| 源 | 数据集 | 特点 |
|---|---|---|
| UCI (12 拉 11 用) | auto_mpg 392 / automobile 205 / concrete 1030 / slump 103 / energy_eff 768 / yacht 308 / ccpp 9568 / air_quality 9358 / gas_turbine 36733 / tetouan_power 52417 / steel_energy 35040 行 | 大样本工程回归台架 |
| NIST StRD (13) | Chwirut1/2, Misra1a-d, DanWood, Kirby2, Hahn1, Roszman1, Thurber, Eckerle4, Bennett5 | 真实测量标定曲线, 官方有认证非线性模型(非幂律) |
| NIST WebBook (26) | 烷烃/烯烃/芳烃/醇/卤代/无机/制冷剂共 26 组 ΔvapH(或 ΔsubH/ΔfusH) vs T, n=10-34 | 相变焓逐点测量; 临界点趋零, **已知非幂律** |
| USGS NWIS 瞬时值 (12) | potomac/nooksack/guadalupe + 密苏里/科罗拉多/哥伦比亚/萨斯奎哈纳/特拉华/密西西比/萨克拉门托/阿查法拉亚/康涅狄格 12 站 2023 | 率定曲线 Q=C(h-h0)^n, 工程界最经典幂律 |
| Open-Meteo/ERA5 (10) | 北海/德国湾/波罗的海/堪萨斯/德州/Altamont/艾奥瓦/内蒙古/黄海/巴塔哥尼亚 | 风速切变 v100~v10 幂律正考场 + 阵风系数; **再分析数据已如实标注** |
| USGS FDSN 地震 (1) | 全球 M≥5.5, 2000-2024, 12334 事件 | felt~mag 感知超线性 / nst~mag 检测能力 / cdi~mag 烈度换算 |
| RDatasets (14) | pressure/cars/rock/trees/attenu/housing/electricity + stormer(粘度计)/whiteside(保温)/cement_hald(水化热)/petrol(炼厂)/women(异速)/airquality(臭氧)/puromycin(酶动力学) | 经典物理/工程小样本(cars: dist~v², puromycin: Michaelis-Menten) |
| OpenML (9 拉 3 用, 三期) | ailerons 13750(F-16 飞控: 舵偏~攻角)/cpu_act 8192(服务器遥测: 利用率~调用强度)/elusage 55(住宅用电~气温) | 航空控制/系统工程遥测; 归一化量级, 只做形式检验数值台架 |
| EPA fueleconomy (1, 三期) | 全美车辆库 1984-2025 燃料车主样本 48,624 行 | 油耗/尾气 CO2 vs 排量; co2~mpg 化学计量近精确反比(理论指数 -1)的大样本 |

已拉取未清洗/已跳过的: UCI 994(微燃机, 无幂律价值)、UCI 316(船用燃机,
**仿真器数据**违背真实数据纪律)、494/882(波浪能, 数值建模)、
USGS 03255000 辛辛那提(全年瞬时流量序列未发布, 0 行留痕跳过)、
WebBook 相变表 <10 行的 17 个流体页(o_xylene 8 行/nonane 9 行等, 小样本拒检线以下)、
OpenML elevators(全正列只剩 Goal 一个, 舵面历史 SaTime* 为负符号约定)、
delta_elevators/delta_ailerons(全符号列)、puma8NH(关节角/力矩全符号)、
cloud(名义列主导)、pol(f1-f48 多常数退化列)——后五者 log 空间不可用, 抓取留档不注册。
**不可达源留痕**: EDP/Engie 风机 SCADA 开放数据(opendatasoft 平台)在本机
网络 TLS 握手被断(2026-09-06 实测), 功率曲线考场待网络环境允许后补。

## 预期物理行为(考场设计)

- **幂律近似成立(H0 侧)**: ccpp PE~AT, gas_turbine TEY~CDP, cars
  dist~speed(v²), trees Volume~Girth(D²), auto_mpg mpg~weight,
  **meteo v100~v10(风工程切变律, 二期新增的正考场主菜)**,
  steel usage~CO2(化学计量), usgs 12 站 Q~h(带 h0 截距轻微失配),
  epa comb_mpg~displ(道路负载负幂律)与 co2~mpg(理论指数 -1, 三期)
- **已知非幂律(H1 侧, 应报警)**: webbook 26 组 ΔvapH~T(临界趋零),
  yacht resid_resist~froude(高速段指数漂移), concrete strength~age_day
  (固化趋饱和), puromycin rate~conc(Michaelis-Menten 饱和, 二期),
  tetouan zone~temp_c(U 形双峰, 二期), ccpp PE~AP, usgs_eq felt~mag
  (指数型), housing price~lotsize, ailerons se_pos~alpha(配平关系非
  幂律, 三期), cpu_act usr~scall(排队饱和, 三期)
- **诚实失败考场**: bennett5(y 全负), roszman1(x 全负), danwood(n=6),
  cement_hald(n=13), women(n=15)——小样本/符号问题的拒检留痕
- 全量冒烟(seed=42, 大样本抽 n≤800): **147/150 格成功**,
  fc 拒绝率 0.660 / fitted-F 0.755 / exog-HC3 0.639;
  fc 的 10 个失败格全部是 n≤11 的小样本诚实拒检(checked_false_small_n)。
  三个真实数据新素材: ① meteo v100~v10 在 n=800 下 fitted-F 10/10 站
  全拒而 fc 保留 4/10(切变律仅近似幂律+RESET 大样本过敏感的现场版);
  ② puromycin rate~conc: fc=0.011 检出 MM 曲率而 fitted-F=0.069 漏检
  (小样本功效真实数据回声); ③ stormer Time~Viscosity+Wt: 检定律
  time·Wt=k(μ+c) 含交互项, 纯幂律双特征两法齐拒(真实交互结构考场)。
  三期新分歧格: epa co2~mpg 两法齐拒(近精确反比律上的车级混合曲率),
  epa comb_mpg~displ 两法齐不拒; automobile horsepower~engine_size
  fc=0.878 不拒而 fitted-F=1.2e-4 拒(异方差散布下 fc 保守性的现场版)。
  详见 `smoke_results.csv`。

## 用法

```bash
"$PYTHON" real_data/fetch_data.py      # 幂等抓取原始数据(raw/, 约 115MB)
"$PYTHON" real_data/prepare_data.py    # 清洗 -> clean/*.csv + manifest.json
"$PYTHON" real_data/smoke_test.py --max-n 800   # 全关系冒烟
```

程序化使用:

```python
import json, pandas as pd
mf = json.load(open("real_data/clean/manifest.json", encoding="utf-8"))
for name, e in mf.items():
    df = pd.read_csv(f"real_data/clean/{name}.csv")
    for r in e["candidate_relations"]:
        ...  # r["target"], r["features"], r["note"]
```

或走实验台: `common.run_all_methods(df_sub, target, features, seed)`
(在仓库根导入;大样本建议先抽 n≤2000, 排列检验 n_boot=200 起)。

## 纪律

- 原始数据只进 `raw/`, 清洗**不插值不填缺**: 缺失/负值行丢行并在
  manifest `notes` 留痕;单位换算仅两处(CCPP 无, USGS 原样英制,
  whiteside °C→K 显式留痕)。
- 候选关系只声明物理上有讨论意义的组合, 不为凑数发明关系;
  ERA5 再分析、传感器哨兵值(-200)、保温前后混合等数据性质全部在
  manifest notes 如实标注, 不冒充直接测量。
- 真实数据无已知真值: 判定差异只呈现不判胜负, 功效/校准结论只能来自
  合成台架(`runs/full/`)。
- USGS 瞬时值/Open-Meteo 逐时值时间自相关未处理, 排列检验独立性
  假设是近似(manifest 已注)。
