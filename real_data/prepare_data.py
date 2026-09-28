"""prepare_data.py — 把 real_data/raw/ 的原始数据清洗为 form_check 可用格式。

产出:
  real_data/clean/<name>.csv     每列物理量, 数值型; 候选关系列全为正值的行保留
                                 (原始行数也记录在 manifest, 负值/缺失行在运行时过滤)
  real_data/clean/manifest.json  每个数据集: 来源 URL/描述/列单位/候选幂律关系/
                                 原始文件 sha256/注意事项

清洗纪律(与框架哲学一致):
  - 不插值、不填缺失、不改单位含义; 单位换算显式写注释并在 manifest 留痕
  - 候选关系只声明物理上有讨论意义的 (target, features); 不为凑数发明关系
  - 已知非幂律形状(蒸气压/蒸发焓/固化龄期)如实标注, 它们是 form_check 的
    正考场, 不是数据缺陷
"""
from __future__ import annotations

import hashlib
import io
import json
import os
import re
import zipfile

import lxml.html
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw")
CLEAN = os.path.join(HERE, "clean")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


ENTRIES: dict[str, dict] = {}


def register(name, df, *, source, url, description, columns, relations,
             raw_files, notes=()):
    path = os.path.join(CLEAN, f"{name}.csv")
    df.to_csv(path, index=False, encoding="utf-8")
    ENTRIES[name] = {
        "source": source,
        "source_url": url,
        "description": description,
        "columns": columns,
        "n_rows": int(len(df)),
        "candidate_relations": relations,
        "raw_files": [{"path": os.path.relpath(p, HERE), "sha256": sha256(p)}
                      for p in raw_files],
        "notes": list(notes),
    }
    print(f"  clean {name}: {len(df)} rows, "
          f"{len(relations)} candidate relation(s)")


def rel(target, features, note=""):
    return {"target": target, "features": list(features), "note": note}


# ------------------------------------------------------------------ UCI
def prep_auto_mpg():
    d = os.path.join(RAW, "uci", "9_auto_mpg")
    df = pd.read_csv(os.path.join(d, "data.csv"), na_values=["?"])
    df = df.drop(columns=["car_name", "origin", "model_year", "cylinders",
                          "acceleration"]).dropna()
    register(
        "uci_auto_mpg", df.reset_index(drop=True), source="uci",
        url="https://archive.ics.uci.edu/dataset/9/auto+mpg",
        description="1970-1982 款汽车: 油耗(mpg) vs 排量/马力/整备质量。"
                    "Physical intuitiveness: heavier/more powerful -> lower "
                    "mpg, log-log 中近似线性但整备质量与排量强共线。",
        columns={"displacement": "in**3", "horsepower": "hp",
                 "weight": "lb", "mpg": "mile/gal (target)"},
        relations=[
            rel("mpg", ["weight"], "经典: 重车费油"),
            rel("mpg", ["displacement"], "排量 vs 油耗"),
            rel("mpg", ["horsepower"], "马力 vs 油耗"),
            rel("mpg", ["displacement", "weight", "horsepower"],
                "三特征, 强共线考共线性告警"),
        ],
        raw_files=[os.path.join(d, "data.csv")],
        notes=["原始 398 行, horsepower 含 6 个 '?' 缺失, 直接丢行不插值",
               "cylinders/origin/model_year/acceleration 为离散或弱物理列, 清洗期剔除"])


def prep_concrete():
    d = os.path.join(RAW, "uci", "165_concrete_strength")
    path = os.path.join(d, "data.csv")
    df = pd.read_csv(path)
    ren = {"Cement": "cement_kg_m3", "Blast Furnace Slag": "slag_kg_m3",
           "Fly Ash": "flyash_kg_m3", "Water": "water_kg_m3",
           "Superplasticizer": "superplast_kg_m3",
           "Coarse Aggregate": "coarseagg_kg_m3",
           "Fine Aggregate": "fineagg_kg_m3", "Age": "age_day",
           "Concrete compressive strength": "strength_mpa"}
    df = df.rename(columns=ren)
    register(
        "uci_concrete_strength", df, source="uci",
        url="https://archive.ics.uci.edu/dataset/165/concrete+compressive+strength",
        description="1030 组混凝土试块抗压强度 vs 配比组分与龄期 (Yeh 1998)。"
                    "强度~龄期是固化增长曲线(亚线性、趋饱和), 幂律只在中段近似。",
        columns={c: "kg in a m3 mixture" for c in ren if c.endswith("_m3")
                 and c != "age_day"} | {"age_day": "day",
                                        "strength_mpa": "MPa (target)"},
        relations=[
            rel("strength_mpa", ["cement_kg_m3"], "Abrams/Bolomey 类单变量"),
            rel("strength_mpa", ["water_kg_m3"], "水灰比方向, 预期负指数"),
            rel("strength_mpa", ["age_day"], "固化龄期, 已知亚线性趋饱和"),
            rel("strength_mpa", ["cement_kg_m3", "water_kg_m3", "age_day"],
                "三特征多变量"),
        ],
        raw_files=[path],
        notes=["含大量 0 组分(配比本就可为 0); 框架逐特征降级路径的天然考场"])


def prep_energy_efficiency():
    d = os.path.join(RAW, "uci", "242_energy_efficiency")
    path = os.path.join(d, "data.csv")
    df = pd.read_csv(path)
    ren = {"X1": "rel_compact", "X2": "surface_area", "X3": "wall_area",
           "X4": "roof_area", "X5": "overall_height",
           "Y1": "heating_load", "Y2": "cooling_load"}
    df = df.rename(columns=ren).drop(
        columns=["X6", "X7", "X8"])
    register(
        "uci_energy_efficiency", df, source="uci",
        url="https://archive.ics.uci.edu/dataset/242/energy+efficiency",
        description="768 栋矩形建筑能耗仿真 (Tsanas & Xifara 2012): "
                    "采暖/制冷负荷 vs 形状参数。制冷负荷~glazing 在原研究中非单调。",
        columns={"rel_compact": "-", "surface_area": "m**2",
                 "wall_area": "m**2", "roof_area": "m**2",
                 "heating_load": "kWh/m**2 (target)",
                 "cooling_load": "kWh/m**2 (target)"},
        relations=[
            rel("heating_load", ["surface_area"]),
            rel("heating_load", ["roof_area"]),
            rel("cooling_load", ["surface_area"]),
            rel("heating_load", ["surface_area", "wall_area", "roof_area"],
                "三面积特征, 相互强相关"),
        ],
        raw_files=[path],
        notes=["X5 overall_height 仅 2 个水平(3.5/7 m), 做特征会阶跃, 已剔除",
               "X6/X7/X8 为离散编码列, 已剔除"])


def prep_yacht():
    d = os.path.join(RAW, "uci", "243_yacht_hydro")
    path = os.path.join(d, "yacht_hydrodynamics.data")
    df = pd.read_csv(path, sep=r"\s+", header=None, names=[
        "long_pos", "pris_coef", "ldr_ratio", "bdr_ratio", "lbl_ratio",
        "froude", "resid_resist"])
    register(
        "uci_yacht_hydro", df, source="uci",
        url="https://archive.ics.uci.edu/dataset/243/yacht+hydrodynamics",
        description="308 组帆船船模拖曳水池试验 (Tsur & Tuller, Delft 系列): "
                    "剩余阻力/单位排水量 vs Froude 数与船型参数。"
                    "高速段阻力~速度指数可达 5-6, 是 log-log 强弯曲的经典考场。",
        columns={"long_pos": "-", "pris_coef": "-", "ldr_ratio": "L/D**3",
                 "bdr_ratio": "B/T", "lbl_ratio": "L/B",
                 "froude": "V/sqrt(gL) (feature)",
                 "resid_resist": "base unit (target)"},
        relations=[
            rel("resid_resist", ["froude"], "阻力~Froude, 高速段强弯曲"),
            rel("resid_resist", ["froude", "ldr_ratio"], "双特征"),
        ],
        raw_files=[path],
        notes=["long_pos 有负值, 不进特征; 各船型参数在每条船内为常数, "
               "跨船才变化(分组共线结构)"])


def prep_ccpp():
    d = os.path.join(RAW, "uci", "294_ccpp")
    path = os.path.join(d, "data.csv")
    df = pd.read_csv(path)
    register(
        "uci_ccpp", df, source="uci",
        url="https://archive.ics.uci.edu/dataset/294/combined+cycle+power+plant",
        description="9568 条联合循环电厂满负荷净电力输出 (Tüfekci 2014, "
                    "6 年 1 小时平均): PE vs 环境温度/真空/气压/湿度。"
                    "PE~AT 近线性(汽轮机效率), 是幂律 H0 的真实代理。",
        columns={"AT": "C ambient temperature", "V": "cm Hg exhaust vacuum",
                 "AP": "mbar ambient pressure", "RH": "% relative humidity",
                 "PE": "MW net electrical output (target)"},
        relations=[
            rel("PE", ["AT"], "效率~温度, 近线性幂律(指数约 0.9-1)"),
            rel("PE", ["V"], "真空度影响"),
            rel("PE", ["AP"], "气压影响, 弱"),
            rel("PE", ["AT", "V"], "双特征"),
        ],
        raw_files=[path],
        notes=["全列正值; 大样本(n=9568)下排列检验较慢, 冒烟只抽子样"])


def prep_gas_turbine():
    d = os.path.join(RAW, "uci", "551_gas_turbine_emissions")
    path = os.path.join(d, "data.csv")
    df = pd.read_csv(path)
    register(
        "uci_gas_turbine_emissions", df, source="uci",
        url="https://archive.ics.uci.edu/dataset/551/gas+turbine+co+and+nox+emission+data+set",
        description="燃机 CO/NOx 排放与热效率遥测 (Kaya et al. 2019, "
                    "2011-2015, ~36733 行): 环境参数+燃烧参数。"
                    "TEY~CDP 近线性; CO~TIT 呈化学反应非线性。",
        columns={"AT": "C ambient", "AP": "mbar", "AH": "% absolute humidity",
                 "AFDP": "mbar air-fuel dp", "GTEP": "mbar turbine ep",
                 "TIT": "K turbine inlet temp", "TAT": "K turbine after temp",
                 "TEY": "MWh energy yield (target)", "CDP": "mbar compressor dp",
                 "CO": "mg/m3 (target)", "NOX": "mg/m3 (target)"},
        relations=[
            rel("TEY", ["CDP"], "能量产出~压气机出口压力, 近线性"),
            rel("TEY", ["GTEP"], "turbine energy pressure"),
            rel("NOX", ["TIT"], "NOx~燃烧温度, 热力型生成非线性"),
            rel("CO", ["TIT"], "CO~燃烧温度, 完全燃烧方向"),
        ],
        raw_files=[path],
        notes=["AT 全期含负值(冬季), 相关行会被运行时正值过滤剔除; CO 含 0 值",
               "giga_conflicts: 无; 大样本, 冒烟只抽子样"])


def prep_concrete_slump():
    d = os.path.join(RAW, "uci", "182_concrete_slump")
    path = os.path.join(d, "slump_test.data")
    df = pd.read_csv(path)
    ren = {"Cement": "cement", "Slag": "slag", "Fly ash": "fly_ash",
           "Water": "water", "SP": "superplast",
           "Coarse Aggr.": "coarse_agg", "Fine Aggr.": "fine_agg",
           "SLUMP(cm)": "slump_cm", "FLOW(cm)": "flow_cm",
           "Compressive Strength (28-day)(Mpa)": "strength_mpa"}
    df = df.drop(columns=["No"]).rename(columns=ren)
    register(
        "uci_concrete_slump", df, source="uci",
        url="https://archive.ics.uci.edu/dataset/182/concrete+slump+test",
        description="103 组高强混凝土坍落度试验 (Yeh 2007): 坍落/流动/"
                    "28 天抗压强度 vs 七组分配比。强度~水泥/水是配比设计"
                    "经典方向; 流动性~用水率近线性。",
        columns={**{c: "kg/m3 batch component" for c in
                    ["cement", "slag", "fly_ash", "water", "superplast",
                     "coarse_agg", "fine_agg"]},
                 "slump_cm": "cm", "flow_cm": "cm (target)",
                 "strength_mpa": "MPa (target)"},
        relations=[
            rel("strength_mpa", ["cement"], "强度~水泥用量"),
            rel("strength_mpa", ["water"], "强度~用水量, 预期负指数方向"),
            rel("flow_cm", ["water"], "流动性~用水量近线性"),
        ],
        raw_files=[path],
        notes=["slump_cm 含 0 值(高流态), 该列不进候选关系避免正值过滤失真",
               "zip 内另有 slump_test.test(40 行), 不用"])


def prep_air_quality():
    d = os.path.join(RAW, "uci", "360_air_quality")
    path = os.path.join(d, "data.csv")
    df = pd.read_csv(path, na_values=[-200])  # -200 是文档化缺失哨兵
    ren = {"CO(GT)": "co_gt", "PT08.S1(CO)": "pt08s1_co",
           "NMHC(GT)": "nmhc_gt", "C6H6(GT)": "c6h6_gt",
           "PT08.S2(NMHC)": "pt08s2_nmhc", "NOx(GT)": "nox_gt",
           "PT08.S3(NOx)": "pt08s3_nox", "NO2(GT)": "no2_gt",
           "PT08.S4(NO2)": "pt08s4_no2", "PT08.S5(O3)": "pt08s5_o3",
           "T": "temp_c", "RH": "rh_pct", "AH": "ah"}
    df = df.drop(columns=["Date", "Time"]).rename(columns=ren)
    register(
        "uci_air_quality", df, source="uci",
        url="https://archive.ics.uci.edu/dataset/360/air+quality",
        description="意大利某公路边空气污染监测站一年小时值 (De Vito 2008): "
                    "5 种金属氧化物传感器响应 vs 参比分析仪浓度 + 温湿度。"
                    "传感器~浓度是计量校准曲线的天然台架。",
        columns={"co_gt": "mg/m3 reference", "pt08s1_co": "sensor a.u.",
                 "nmhc_gt": "ppb reference", "c6h6_gt": "ug/m3 reference",
                 "pt08s2_nmhc": "sensor a.u.", "nox_gt": "ppb reference",
                 "pt08s3_nox": "sensor a.u.", "no2_gt": "ug/m3 reference",
                 "pt08s4_no2": "sensor a.u.", "pt08s5_o3": "sensor a.u.",
                 "temp_c": "C", "rh_pct": "%", "ah": "absolute humidity"},
        relations=[
            rel("pt08s1_co", ["co_gt"], "锡氧化物传感器 vs CO 参比: 校准曲线"),
            rel("pt08s2_nmhc", ["c6h6_gt"], "钨氧化物传感器 vs 苯参比"),
            rel("pt08s2_nmhc", ["c6h6_gt", "temp_c"],
                "加温度特征: 传感器温度交叉敏感"),
        ],
        raw_files=[path],
        notes=["-200 为数据方文档化的缺失哨兵, 读入即 NaN(不插值)",
               "nmhc_gt 缺失约 90%, 仅存档不进关系",
               "小时均值, 时间自相关未处理"])


def prep_tetouan():
    d = os.path.join(RAW, "uci", "849_tetouan_power")
    path = os.path.join(d, "data.csv")
    df = pd.read_csv(path)
    # Zone 2/3 列名含双空格, 归一化内部空白后再改名
    df.columns = [" ".join(c.split()) for c in df.columns]
    ren = {"Temperature": "temp_c", "Humidity": "humidity_pct",
           "Wind Speed": "wind_speed_ms",
           "general diffuse flows": "gen_diffuse",
           "diffuse flows": "diffuse",
           "Zone 1 Power Consumption": "zone1_kw",
           "Zone 2 Power Consumption": "zone2_kw",
           "Zone 3 Power Consumption": "zone3_kw"}
    df = df.drop(columns=["DateTime"]).rename(columns=ren)
    register(
        "uci_tetouan_power", df, source="uci",
        url="https://archive.ics.uci.edu/dataset/849/power+consumption+of+tetouan+city",
        description="摩洛哥 Tetouan 市三配区 2017 全年 10min 电力消费 "
                    "(52,417 行) vs 温/湿/风/辐照度。负荷~气温呈冬夏双峰 "
                    "U 形(采暖+制冷), 对幂律是连续弯曲考场。",
        columns={"temp_c": "C", "humidity_pct": "%",
                 "wind_speed_ms": "m/s", "gen_diffuse": "irradiance a.u.",
                 "diffuse": "irradiance a.u.",
                 "zone1_kw": "kW (target)", "zone2_kw": "kW (target)",
                 "zone3_kw": "kW (target)"},
        relations=[
            rel("zone1_kw", ["temp_c"], "负荷~气温: U 形双峰, 已知非单调"),
            rel("zone1_kw", ["temp_c", "humidity_pct", "wind_speed_ms"],
                "三气象特征多变量"),
            rel("zone3_kw", ["temp_c"], "第三区(商业为主)负荷~气温"),
        ],
        raw_files=[path],
        notes=["diffuse flows 原文未标单位, 列名不带单位声明",
               "10min 分辨率时间自相关; 大样本冒烟抽子样"])


def prep_steel_energy():
    d = os.path.join(RAW, "uci", "851_steel_energy")
    path = os.path.join(d, "data.csv")
    df = pd.read_csv(path, encoding="utf-8-sig")  # 去 BOM
    ren = {"Usage_kWh": "usage_kwh",
           "Lagging_Current_Reactive.Power_kVarh": "lag_reactive_kvarh",
           "Leading_Current_Reactive.Power_kVarh": "lead_reactive_kvarh",
           "CO2(tCO2)": "co2_tco2",
           "Lagging_Current_Power_Factor": "lag_pf",
           "Leading_Current_Power_Factor": "lead_pf", "NSM": "nsm"}
    df = df.drop(columns=["date", "WeekStatus", "Day_of_week", "Load_Type"]
                 ).rename(columns=ren)
    register(
        "uci_steel_energy", df, source="uci",
        url="https://archive.ics.uci.edu/dataset/851/steel+industry+energy+consumption",
        description="韩国某小型钢厂 2018 年 15min 智能电表遥测 (35,040 行): "
                    "有功/无功/功率因数/CO2。能耗~CO2 由燃烧化学计量近比例"
                    "耦合; 无功~负载是感应电机饱和曲线。",
        columns={"usage_kwh": "kWh (target)",
                 "lag_reactive_kvarh": "kVarh", "lead_reactive_kvarh": "kVarh",
                 "co2_tco2": "t CO2", "lag_pf": "%", "lead_pf": "%",
                 "nsm": "seconds from midnight"},
        relations=[
            rel("usage_kwh", ["co2_tco2"], "能耗~CO2: 化学计量近比例(含0行被滤)"),
            rel("usage_kwh", ["lag_reactive_kvarh"], "感应电机无功~负载"),
            rel("usage_kwh", ["co2_tco2", "lag_reactive_kvarh"], "双特征"),
        ],
        raw_files=[path],
        notes=["轻载行 CO2=0(电弧炉间歇), 正值过滤后保留重载段——留痕不掩盖",
               "nsm 呈日内周期, 不进关系(非物理幂律候选)"])


# ------------------------------------------------------------- NIST StRD
NIST_DESC = {
    "Chwirut1": "超声标定: 超声响应 y vs 金属距离 x (n=214, NIST)",
    "Chwirut2": "超声标定短版 (n=54)",
    "Misra1a": "牙科根管种植体/骨界面 (Misra 1978), n=14",
    "Misra1b": "同 Misra1a 模型的第二组, n=40",
    "Misra1c": "同 Misra1a 模型的第三组, n=24",
    "Misra1d": "同 Misra1a 模型的第四组, n=44",
    "DanWood": "能源装置: 能耗 kWh y vs 运行天数 x (Daniel & Wood), n=6",
    "Kirby2": "微荧光分析 (Kirby 1976), n=151",
    "Hahn1": "铜镍合金热膨胀系数 y vs 温度 x (Hahn 1970s), n=236",
    "Roszman1": "电子量子产额 y vs 频率 x (Roszman), n=176",
    "Thurber": "电子器件电流 y vs 电压 x (Thurber), n=37",
    "Eckerle4": "圆形干涉滤波器透射率 y vs 波长 x (Eckerle), n=33",
    "Bennett5": "磁滞回线类标定 (Bennett 1954?), n=154",
}

# log-空间幂律可测性(数据符号/样本量的诚实声明, 不改数据迁就检验)
NIST_LOG_NOTES = {
    "Misra1a": ["n=14: 低于 run_all_methods 的小样本门槛, fc 拒检(留痕)",
                "官方 n=14, 四个 Misra 集 .dat 均为 14 行(以文件为准)"],
    "Misra1b": ["n=14, 同上"],
    "Misra1c": ["n=14, 同上"],
    "Misra1d": ["n=14, 同上"],
    "DanWood": ["n=6: 低于任何检验的可用样本量, 仅存档"],
    "Bennett5": ["y 全负(-34.8~-31.8): log 不可用, 保留为框架'诚实失败'考场"],
    "Roszman1": ["x 全负(-4869~-464): log 不可用, 同上"],
    "Thurber": ["x 含负值(-3.07~2.20): 正值过滤后 n=11, fc 小样本拒检留痕"],
}


def prep_nist_strd():
    for name, desc in NIST_DESC.items():
        path = os.path.join(RAW, "nist_strd", f"{name}.dat")
        lines = open(path, encoding="ascii").read().splitlines()
        idx = next(i for i, l in enumerate(lines)
                   if l.startswith("Data:") and "y" in l
                   and "Response" not in l)
        header_tokens = lines[idx].split(":")[1].split()
        arr = np.loadtxt(io.StringIO("\n".join(lines[idx + 1:])))
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        cols = (header_tokens + [f"c{i}" for i in range(arr.shape[1])])
        df = pd.DataFrame(arr, columns=cols[:arr.shape[1]])
        feats = [c for c in df.columns if c != "y"]
        relations = [rel("y", feats[:1], "NIST 认证非线性模型的变量, 幂律仅作对照")]
        if len(feats) >= 2:
            relations.append(rel("y", feats, "全特征"))
        register(
            f"nist_{name.lower()}", df, source="nist_strd",
            url=("https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/"
                 f"{name}.dat"),
            description=desc,
            columns={"y": "response (target)", **{c: "predictor"
                                                  for c in feats}},
            relations=relations,
            raw_files=[path],
            notes=["官方认证模型是非线性最小二乘形式, 不是幂律; "
                   "这些集用于考 form_check 在小 n 真实标定数据上的行为",
                   f"变量声明: {header_tokens}",
                   *NIST_LOG_NOTES.get(name, [])])


# ------------------------------------------------------------- WebBook
def prep_webbook():
    wb_dir = os.path.join(RAW, "webbook")
    # ΔvapH/ΔfusH/ΔsubH 等 vs Temperature(K) 逐点表:
    # 表头第一列是量名(带单位), 第二列是 'Temperature (K)'。
    qty_file = {"ΔvapH": "dvaph", "ΔsubH": "dsubh", "ΔfusH": "dfush",
                "ΔfusS": "dfuss"}
    for fname in sorted(os.listdir(wb_dir)):
        if not fname.endswith(".html"):
            continue
        stem = fname[: -len(".html")]
        fluid = stem.removesuffix("_phase").removesuffix("_fluid")
        doc = lxml.html.parse(os.path.join(wb_dir, fname)).getroot()
        for table in doc.findall(".//table"):
            heads = [[c.text_content().strip() for c in r.findall(".//th")]
                     for r in table.findall(".//tr")]
            if not heads or len(heads[0]) < 2:
                continue
            h0 = heads[0]
            if "Temperature (K)" not in h0:
                continue
            qty_raw = h0[0]
            m = re.match(r"^(Δ\w+)\s*\((.+)\)$", qty_raw)
            if not m:
                continue
            qty, unit = m.group(1), m.group(2)
            slug = qty_file.get(qty)
            if slug is None:
                continue
            rows = []
            for tr in table.findall(".//tr")[1:]:
                cells = [c.text_content().strip()
                         for c in tr.findall(".//td")]
                if len(cells) < 2:
                    continue
                try:
                    qv = float(cells[0].split("±")[0])
                    tv = float(cells[1].split("±")[0].split(" to ")[0])
                except ValueError:
                    continue
                if tv > 0 and qv > 0:
                    rows.append((tv, qv))
            if len(rows) < 10:
                continue
            df = pd.DataFrame(rows, columns=["temp_k", f"{slug}_{unit}"])
            col = f"{slug}_{unit}"
            register(
                f"webbook_{fluid}_{slug}", df, source="webbook",
                url=("https://webbook.nist.gov/cgi/cbook.cgi"
                     f"?ID=&Units=SI&Mask=4#Phase-Change"),
                description=f"{fluid}: 相变 {qty} vs 温度逐点测量值 "
                            f"(NIST Chemistry WebBook/TRC 汇编)。"
                            "ΔvapH(T) 在临界点趋零(Watson 型), 对 T 的幂律"
                            "只在中段近似——已知非幂律形状。",
                columns={"temp_k": "K", col: f"{unit} (target)"},
                relations=[rel(col, ["temp_k"],
                               "蒸发焓~温度: 已知非幂律(临界趋零)")],
                raw_files=[os.path.join(wb_dir, fname)],
                notes=["同温度多来源测量点全部保留(重复测量是真实数据的一部分)",
                       "含 '±' 不确定度的取中心值"])


# ------------------------------------------------------------- USGS
def prep_usgs():
    for fname in sorted(os.listdir(os.path.join(RAW, "usgs"))):
        if not fname.startswith("iv_") or not fname.endswith(".json"):
            continue
        path = os.path.join(RAW, "usgs", fname)
        data = json.load(open(path, encoding="utf-8"))
        series = {}
        for s in data["value"]["timeSeries"]:
            code = s["variable"]["variableCode"][0]["value"]
            recs = {}
            for v in s["values"][0]["value"]:
                try:
                    recs[v["dateTime"]] = float(v["value"])
                except ValueError:
                    pass
            series[code] = recs
        if "00060" not in series or "00065" not in series:
            print(f"  skip {fname}: missing Q or h series")
            continue
        q, h = series["00060"], series["00065"]
        common_ts = sorted(set(q) & set(h))
        df = pd.DataFrame({
            "gage_height_ft": [h[t] for t in common_ts],
            "discharge_cfs": [q[t] for t in common_ts],
        })
        site = data["value"]["timeSeries"][0]["sourceInfo"]["siteCode"][0]["value"]
        name = data["value"]["timeSeries"][0]["sourceInfo"]["siteName"]
        df = df[(df["discharge_cfs"] > 0) & (df["gage_height_ft"] > 0)]
        if len(df) < 10:
            print(f"  skip {fname}: 仅 {len(df)} 行可用(Q 序列未发布或全被滤), "
                  "不注册空数据集")
            continue
        register(
            f"usgs_{site}_rating", df.reset_index(drop=True), source="usgs",
            url=("https://waterservices.usgs.gov/nwis/iv/?sites=" f"{site}"
                 "&startDT=2023-01-01&endDT=2023-12-31"
                 "&parameterCd=00060,00065&format=json"),
            description=f"USGS 站 {site} ({name}) 2023 年瞬时值: "
                        "流量 Q vs 水位 h——教科书率定曲线 Q=C(h-h0)^n, "
                        "n≈1.5-2.5, 是工程界最经典的幂律之一。"
                        "h0>0 使纯幂律在低水位端系统性失配。",
            columns={"gage_height_ft": "ft (feature)",
                     "discharge_cfs": "ft**3/s (target)"},
            relations=[rel("discharge_cfs", ["gage_height_ft"],
                           "率定曲线幂律正例(带 h0 截距的轻微失配)")],
            raw_files=[path],
            notes=["Q/h 按时间戳精确配对, 负值/缺失(冰封等)已丢行",
                   "时间自相关未处理: 排列检验独立性假设在此是近似"])


# ------------------------------------------------------------- Open-Meteo
def prep_open_meteo():
    om_dir = os.path.join(RAW, "open_meteo")
    if not os.path.isdir(om_dir):
        return
    for fname in sorted(os.listdir(om_dir)):
        if not fname.endswith(".json"):
            continue
        site = fname.split("_2023")[0]
        path = os.path.join(om_dir, fname)
        js = json.load(open(path, encoding="utf-8"))
        hourly = js.get("hourly", {})
        if "time" not in hourly:
            print(f"  skip {fname}: no hourly block")
            continue
        df = pd.DataFrame({
            "wind_speed_10m": hourly["wind_speed_10m"],
            "wind_speed_100m": hourly["wind_speed_100m"],
            "wind_gusts_10m": hourly["wind_gusts_10m"],
        })
        lat, lon = js.get("latitude"), js.get("longitude")
        register(
            f"meteo_{site}", df, source="open_meteo",
            url=("https://archive-api.open-meteo.com/v1/archive"
                 f"?latitude={lat}&longitude={lon}"
                 "&start_date=2023-01-01&end_date=2023-12-31&timezone=UTC"),
            description=f"Open-Meteo/ERA5 再分析 {site} ({lat}, {lon}) "
                        "2023 逐小时风速。风工程教科书幂律: "
                        "v(z2)=v(z1)·(z2/z1)^α(α≈0.11 海上~0.2 陆上), "
                        "阵风系数 gust/v10 近常数。",
            columns={"wind_speed_10m": "km/h", "wind_speed_100m": "km/h",
                     "wind_gusts_10m": "km/h (target)"},
            relations=[
                rel("wind_speed_100m", ["wind_speed_10m"],
                    "风速切变幂律 v100~v10 (H0 正考场)"),
                rel("wind_gusts_10m", ["wind_speed_10m"],
                    "阵风~平均风速, 阵风系数近常数"),
            ],
            raw_files=[path],
            notes=["ERA5 是观测同化再分析, 非现场测风塔; 剪切指数系统性"
                   "偏低是已知偏差, 不作现场标定用",
                   "静风 0 值行被运行时正值过滤; 单位 km/h 未换算",
                   "逐小时时间自相关, 排列检验独立性假设是近似"])


# ------------------------------------------------------------- USGS 地震
def prep_usgs_eq():
    eq_dir = os.path.join(RAW, "usgs_eq")
    if not os.path.isdir(eq_dir):
        return
    path = os.path.join(eq_dir, "global_m55_2000_2024.geojson")
    if not os.path.exists(path):
        return
    js = json.load(open(path, encoding="utf-8"))
    rows = []
    for f in js["features"]:
        p = f["properties"]
        if p.get("type") != "earthquake":
            continue
        rows.append((p.get("mag"), (f.get("geometry") or {}).get(
            "coordinates", [None, None, None])[2],
            p.get("felt"), p.get("cdi"), p.get("nst"), p.get("rms")))
    df = pd.DataFrame(rows, columns=["mag", "depth_km", "felt", "cdi",
                                     "nst", "rms"])
    register(
        "usgs_eq_global", df, source="usgs_eq",
        url=("https://earthquake.usgs.gov/fdsnws/event/1/query"
             "?format=geojson&starttime=2000-01-01&endtime=2024-12-31"
             "&minmagnitude=5.5&orderby=time-asc"),
        description=f"USGS FDSN 全球 M≥5.5 地震目录 2000-2024 "
                    f"({js['metadata']['count']} 事件): 震级/深度/有感"
                    "报告数/社区烈度/定位台站数。有感报告~震级是感知"
                    "人数×烈度放大的超线性标度; 台站数~震级是检测能力标度。",
        columns={"mag": "magnitude (feature)", "depth_km": "km",
                 "felt": "DYFI report count (target)", "cdi": "MMI-like (target)",
                 "nst": "station count (target)", "rms": "s"},
        relations=[
            rel("felt", ["mag"], "有感报告数~震级: 人口暴露×烈度超线性"),
            rel("nst", ["mag"], "定位台站数~震级: 检测能力标度"),
            rel("cdi", ["mag"], "社区烈度~震级: 能量-烈度换算(指数型)"),
        ],
        raw_files=[path],
        notes=["felt/cdi 仅有人报告的事件非空(NaN 行保留, 运行时按关系 dropna)",
               "目录完整性随年代/区域变化, 检测阈值下截断",
               "事件间独立假设对余震序列是近似"])


# ------------------------------------------------------------- RDatasets
def prep_rdatasets():
    common_url = ("https://vincentarelbundock.github.io/Rdatasets/csv/")

    def rd(name, path, description, columns, relations, notes=()):
        p = os.path.join(RAW, "rdatasets", name + ".csv")
        df = pd.read_csv(p)
        if df.columns[0] in ("rownames", ""):
            df = df.drop(columns=df.columns[0])
        register(f"rdataset_{name}", df, source="rdatasets",
                 url=common_url + path, description=description,
                 columns=columns, relations=relations, raw_files=[p],
                 notes=notes)

    rd("pressure", "datasets/pressure.csv",
       "汞蒸气压 vs 温度 (Weast 手册, n=19): log P 对 1/T 的 Clausius-"
       "Clapeyron 形状, 对 T 的幂律显著弯曲。小 n 名场面。",
       {"temperature": "degC", "pressure": "mm Hg (target)"},
       [rel("pressure", ["temperature"], "蒸气压曲线, 非幂律")],
       notes=["n=19, 小样本校准考场"])

    rd("cars", "datasets/cars.csv",
       "1920 年代汽车: 制动距离 vs 车速 (Ezekiel 1930, n=50)。物理上 "
       "dist ≈ v²/2μg(动能项), 二次律是真结构。",
       {"speed": "mph", "dist": "ft (target)"},
       [rel("dist", ["speed"], "制动距离~速度, 理论指数 2")])

    rd("rock", "datasets/rock.csv",
       "石油工程岩芯: 渗透率 perm vs 面积 area/周长 peri (n=48)。",
       {"area": "px", "peri": "px", "perm": "millidarcy (target)"},
       [rel("perm", ["area"]), rel("perm", ["peri"]),
        rel("perm", ["area", "peri"], "双特征共线")],
       notes=["area/peri 为图像像素单位"])

    rd("trees", "datasets/trees.csv",
       "黑樱桃树: 材积 vs 胸径/树高 (n=31)。几何上 V ∝ D²·H, "
       "接近真实幂律结构(正对照素材)。",
       {"Girth": "in (直径)", "Height": "ft",
        "Volume": "ft**3 (target)"},
       [rel("Volume", ["Girth"], "材积~胸径, 理论指数≈2"),
        rel("Volume", ["Height"], "材积~树高, 理论指数≈1"),
        rel("Volume", ["Girth", "Height"], "双特征, 接近精确幂律")])

    rd("attenu", "datasets/attenu.csv",
       "加州大地震强震记录 (Joyner-Boore 1981, n=182): 峰值地面加速度 vs "
       "震级/震中距。GMPE 类关系的真实样本, 距离衰减近似幂律。",
       {"mag": "Richter", "dist": "km", "accel": "g (target)"},
       [rel("accel", ["dist"], "加速度~震中距, 衰减"),
        rel("accel", ["mag", "dist"], "双特征 GMPE 类")],
       notes=["event/station 为分类列, 清洗期剔除"])

    rd("housing", "Ecdat/Housing.csv",
       "House-sold-in-Windsor: 成交价 vs 地块/结构 (Anglin & Gencay 1996, "
       "n=546)。价格~地块面积是 hedonic 幂律候选。",
       {"price": "CAD (target)", "lotsize": "ft**2",
        "bedrooms": "count", "bathrms": "count", "stories": "count"},
       [rel("price", ["lotsize"]),
        rel("price", ["lotsize", "bedrooms", "bathrms", "stories"],
            "多特征 hedonic")],
       notes=["driveway/recroom 等哑变量列已保留但不进关系"])

    rd("electricity", "Ecdat/Electricity.csv",
       "美国电力公司横截面 (n=158): 总成本 vs 产出与要素价格, "
       "成本函数 cost~q 是经典规模弹性考场。",
       {"cost": "total cost (target)", "q": "output (kWh)",
        "pl": "wage", "pk": "capital price", "pf": "fuel price"},
       [rel("cost", ["q"], "成本~产出: 规模经济考题"),
        rel("cost", ["q", "pf"], "双特征")],
       notes=["要素份额 sl/sk/sf 列已保留但不进关系"])

    # ---- 二期(2026-09-06): 工程/物理经典补充 ----
    rd("stormer", "MASS/stormer.csv",
       "Stormer 旋转粘度计 (n=23): 流出时间 vs 液体粘度与驱动砝码重。"
       "粘度计量学经典小样本; time~viscosity 近过原点线性(幂律 p≈1)。",
       {"Viscosity": "poise", "Wt": "g weight", "Time": "s (target)"},
       [rel("Time", ["Viscosity"], "流出时间~粘度"),
        rel("Time", ["Viscosity", "Wt"], "双特征: 砝码重改变驱动力矩")],
       notes=["Williams 标定模型是 time·Wt = k·(viscosity + c), 幂律作对照"])

    # whiteside 需显式单位换算, 不走通用 rd()
    p_ws = os.path.join(RAW, "rdatasets", "whiteside.csv")
    df_ws = pd.read_csv(p_ws)
    if df_ws.columns[0] in ("rownames", ""):
        df_ws = df_ws.drop(columns=df_ws.columns[0])
    df_ws["temp_k"] = df_ws["Temp"] + 273.15
    df_ws = df_ws.drop(columns=["Temp"])
    register(
        "rdataset_whiteside", df_ws, source="rdatasets",
        url="https://vincentarelbundock.github.io/Rdatasets/csv/MASS/whiteside.csv",
        description="谢菲尔德住宅保温前后燃气消耗 vs 室外温度 (n=56, "
                    "Whiteside 1965)。采暖需求~温差的工程能效关系; "
                    "绝对温标下近线性带截距。",
        columns={"temp_k": "K (converted from C)",
                 "Gas": "1000 ft3/week (target)", "Insul": "Before/After"},
        relations=[rel("Gas", ["temp_k"], "采暖负荷~室外绝对温度, 近线性带截距")],
        raw_files=[p_ws],
        notes=["Temp 原列单位 C 且含负值, 显式换算 temp_k=C+273.15 "
               "(本台架第二处单位换算, 留痕)",
               "Insul 保温前/后混合在同一回归, 列保留供分组核查"])

    rd("cement_hald", "MASS/cement.csv",
       "Hald 水泥数据 (n=13): 四种熟料组分 vs 水化热。多元回归最著名的"
       "强共线小样本之一。",
       {"x1": "% aluminate", "x2": "% silicate",
        "x3": "% ferrite", "x4": "% dicalcium", "y": "cal/g heat (target)"},
       [rel("y", ["x1"], "水化热~铝酸钙"),
        rel("y", ["x1", "x2", "x3", "x4"], "四组分全特征, 强共线考场")],
       notes=["n=13 低于多数检验的可用样本量, 小样本拒检留痕"])

    rd("petrol", "MASS/petrol.csv",
       "Prater 汽油炼厂数据 (n=32): 汽油收率 Y vs 原油比重/蒸气压/"
       "ASTM 馏程。Daniel-Wood 同源, 分组结构(No 列)是场景变量。",
       {"SG": "specific gravity", "VP": "lbf/in2 vapor pressure",
        "V10": "F ASTM 10% point", "EP": "F end point",
        "Y": "% yield (target)"},
       [rel("Y", ["EP"], "收率~终馏点"),
        rel("Y", ["SG", "VP", "V10", "EP"], "四物性特征")],
       notes=["No 列为炼厂分组(A-E), 保留供 fit_scenarios 分组实验"])

    rd("women", "datasets/women.csv",
       "美国 30-39 岁女性平均身高体重 (n=15): 体重~身高^p 的异速律"
       "小样本(工程人因/生物力学口径)。",
       {"height": "in", "weight": "lb (target)"},
       [rel("weight", ["height"], "体重~身高, 理论指数 2-3")],
       notes=["n=15 小样本; 逐身高平均后的理想化数据, 方差被压缩"])

    rd("airquality", "datasets/airquality.csv",
       "纽约 1973 年 5-9 月日均值 (n=153): 臭氧 vs 风/温/辐射。"
       "风速对地面臭氧是稀释扩散负幂律; 温度是光化学正相关。",
       {"Ozone": "ppb (target)", "Solar.R": "langley",
        "Wind": "mph", "Temp": "F"},
       [rel("Ozone", ["Wind"], "臭氧~风速: 稀释扩散负指数"),
        rel("Ozone", ["Temp"], "臭氧~温度: 光化学"),
        rel("Ozone", ["Wind", "Temp", "Solar.R"], "三特征多变量")],
       notes=["Ozone/Solar.R 含 NA, 运行时 dropna 不插值",
              "Month/Day 列保留不进关系"])

    rd("puromycin", "datasets/Puromycin.csv",
       "嘌呤霉素酶促反应 (n=23): 反应速率 vs 底物浓度。Michaelis-Menten "
       "饱和曲线是幂律的教科书反例(低浓度近线性、高浓度趋平台)。",
       {"conc": "ppm", "rate": "counts/min (target)", "state": "treated factor"},
       [rel("rate", ["conc"], "Michaelis-Menten 饱和, 已知非幂律")],
       notes=["state(treated/untreated)混合在同一回归, 列保留供分组"])


# ------------------------------------------------------------- OpenML
def _read_arff(path):
    """稠密数值 ARFF → (列名, DataFrame)。名义/字符串值解析为 NaN(留痕剔除)。"""
    with open(path, encoding="utf-8", errors="replace") as fh:
        txt = fh.read()
    attrs = re.findall(r"@attribute\s+'?\"?([^'\"\s]+)'?\"?\s", txt,
                       re.IGNORECASE)
    di = txt.lower().index("@data")
    rows = []
    for raw_ln in txt[di + 5:].strip().splitlines():
        ln = raw_ln.strip()
        if not ln or ln.startswith("%"):
            continue
        vals = []
        for tok in ln.split(","):
            v = tok.strip().strip("'\"")
            try:
                vals.append(float(v))
            except ValueError:  # 名义值(seeded/season 等)
                vals.append(np.nan)
        rows.append(vals)
    return pd.DataFrame(np.array(rows), columns=attrs)


def prep_openml():
    om_dir = os.path.join(RAW, "openml")
    if not os.path.isdir(om_dir):
        return

    # elevators: F-16 纵向控制 (Jaroszewicz & Simovici 2002 域)。
    # 抓取留档但不注册: 全正列只剩 Goal 一个(舵面历史 SaTime* 为负符号
    # 约定), 无成对正关系, log 空间不可用(见 README)。
    # ailerons: 同域横向通道。alpha(攻角)与 Se(舵面位置)全正。
    df = _read_arff(os.path.join(om_dir, "ailerons.arff"))
    df = df.rename(columns={"Se": "se_pos",
                            "SeTime1": "se_t1"})
    register(
        "openml_ailerons", df[["curPitch", "alpha", "se_pos", "se_t1",
                               "goal"]],
        source="openml",
        url="https://www.openml.org/d/296",
        description="F-16 飞控横向通道 13750 条记录 (同 elevators 域): "
                    "攻角 alpha 与舵面位置 se_pos 为正量级。"
                    "se_pos~alpha 是纵向配平的真实气动/控制关系。",
        columns={"curPitch": "pitch (mostly positive)", "alpha": "AoA (rad-ish)",
                 "se_pos": "control surface position (target)",
                 "se_t1": "position t-1", "goal": "roll goal (signed, archive)"},
        relations=[
            rel("se_pos", ["alpha"], "舵面偏度~攻角: 配平关系"),
            rel("se_pos", ["se_t1"], "舵位一阶自相关(伺服速率限制)"),
            rel("se_pos", ["alpha", "se_t1"], "双特征"),
        ],
        raw_files=[os.path.join(om_dir, "ailerons.arff")],
        notes=["归一化量级, 单位不可恢复(同 elevators)",
               "goal 列全负(符号约定), 仅存档不进关系"])

    # cpu_act: 服务器系统遥测 (8192 条 5s 采样)。
    df = _read_arff(os.path.join(om_dir, "cpu_act.arff"))
    df = df[["lread", "scall", "rchar", "freemem", "runqsz", "usr"]].rename(
        columns={"usr": "usr_pct"})
    register(
        "openml_cpu_act", df, source="openml",
        url="https://www.openml.org/d/197",
        description="Unix 服务器系统遥测 8192 条 (每 5s, CMU statlib 系): "
                    "CPU 用户态利用率 vs 系统调用/IO 流量/内存/运行队列。"
                    "利用率~负载是排队论型饱和曲线, 幂律只在中段近似。",
        columns={"lread": "reads/s", "scall": "syscalls/s",
                 "rchar": "chars read/s", "freemem": "pages",
                 "runqsz": "run queue size", "usr_pct": "% CPU user (target)"},
        relations=[
            rel("usr_pct", ["scall"], "利用率~系统调用强度: 饱和曲线"),
            rel("usr_pct", ["rchar"], "利用率~读流量"),
            rel("usr_pct", ["scall", "freemem", "runqsz"],
                "多特征: 调用/内存/运行队列"),
        ],
        raw_files=[os.path.join(om_dir, "cpu_act.arff")],
        notes=["usr=0 的空闲行约 3%, 运行时正值过滤",
               "lread 含 0 值(8%), 该列不进关系"])

    # elusage: 55 个月度电力用量 vs 气温 (CMU statlib)。
    df = _read_arff(os.path.join(om_dir, "elusage.arff"))
    df = df.rename(columns={"average_temperature": "temp_f", "month": "month",
                            "average_electricity_usage": "usage"})
    register(
        "openml_elusage", df, source="openml",
        url="https://www.openml.org/d/228",
        description="某住宅 55 个月平均电力用量 vs 月均气温 (CMU statlib "
                    "elusage): 用电~气温呈冬夏双峰 U 形(采暖+制冷), "
                    "对幂律是连续弯曲考场。",
        columns={"temp_f": "F", "month": "month index", "usage": "kWh (target)"},
        relations=[
            rel("usage", ["temp_f"], "用电~气温: U 形非单调"),
            rel("usage", ["temp_f", "month"], "加月份协变量(季节)"),
        ],
        raw_files=[os.path.join(om_dir, "elusage.arff")],
        notes=["n=55 小样本; 用量单位原文未明示, 按 kWh 记",
               ("elevators/delta_elevators/delta_ailerons/puma8NH/cloud/pol"
                " 已抓取留档但不注册: 全符号关节量/舵面历史负符号约定/"
                "退化归一化特征(f1-f48 多常数列)/名义列主导, "
                "log 空间不可用(见 README)")])


# ------------------------------------------------------------- EPA 油耗
def prep_epa_fueleconomy():
    zpath = os.path.join(RAW, "epa_fueleconomy", "vehicles.csv.zip")
    if not os.path.exists(zpath):
        return
    zf = zipfile.ZipFile(zpath)
    df = pd.read_csv(zf.open("vehicles.csv"), low_memory=False)
    # 只保留燃料车主行(fuelType1 非空): EV/PHEV 无排量口径
    df = df[df["fuelType1"].notna()]
    df = df[df["displ"] > 0]  # displ=0 是占位错误行
    df = df[["year", "displ", "cylinders", "city08", "highway08",
             "comb08", "co2TailpipeGpm"]].rename(columns={
                 "displ": "displ_l", "city08": "city_mpg",
                 "highway08": "hwy_mpg", "comb08": "comb_mpg",
                 "co2TailpipeGpm": "co2_g_per_mi"})
    register(
        "epa_vehicles", df.reset_index(drop=True), source="epa_fueleconomy",
        url="https://www.fueleconomy.gov/feg/download.shtml",
        description=f"EPA/DOE 全美车辆数据库 1984-2025 燃料车主样本 "
                    f"({len(df)} 行): 油耗经济性/尾气 CO2 vs 排量与缸数。"
                    "comb_mpg~displ 是道路负载(滚阻+气动)负幂律; "
                    "co2~mpg 由燃烧化学计量近精确反比(指数≈-1)——"
                    "少有的已知指数真实大样本。",
        columns={"year": "model year", "displ_l": "L",
                 "cylinders": "count", "city_mpg": "MPG",
                 "hwy_mpg": "MPG", "comb_mpg": "MPG (target)",
                 "co2_g_per_mi": "g/mi (target)"},
        relations=[
            rel("comb_mpg", ["displ_l"], "油耗经济性~排量: 负幂律"),
            rel("co2_g_per_mi", ["displ_l"], "尾气 CO2~排量: 正幂律"),
            rel("co2_g_per_mi", ["comb_mpg"],
                "CO2~MPG: 化学计量近精确反比, 理论指数 -1"),
            rel("comb_mpg", ["displ_l", "cylinders"],
                "排量+缸数(缸数离散 2-16)"),
        ],
        raw_files=[zpath],
        notes=["EV/PHEV/双燃料行(无 fuelType1 主燃料排量口径)剔除",
               "displ=0 的占位错误行剔除留痕; 少数 co2=0 行运行时正值过滤",
               "同一动力总成多年复测, 行间非独立(同型重复)"])


# ------------------------------------------------------------- UCI 汽车
def prep_automobile():
    d = os.path.join(RAW, "uci", "10_automobile")
    path = os.path.join(d, "data.csv")
    df = pd.read_csv(path, na_values=["?"])
    ren = {"wheel-base": "wheel_base_in", "length": "length_in",
           "width": "width_in", "height": "height_in",
           "curb-weight": "curb_weight_lb", "engine-size": "engine_size_cc",
           "compression-ratio": "compression_ratio", "peak-rpm": "peak_rpm",
           "city-mpg": "city_mpg", "highway-mpg": "hwy_mpg",
           "price": "price_usd"}
    df = df[["wheel-base", "length", "width", "height", "curb-weight",
             "engine-size", "bore", "stroke", "compression-ratio",
             "horsepower", "peak-rpm", "city-mpg", "highway-mpg",
             "price"]].rename(columns=ren)
    register(
        "uci_automobile", df, source="uci",
        url="https://archive.ics.uci.edu/dataset/10/automobile",
        description="1985 年款进口轿车 205 辆 (Schlimmer): 价格/马力/油耗 vs "
                    "尺寸与整备质量。马力~排量近比例(进气量上限); "
                    "油耗~车重负指数(道路负载)。",
        columns={"wheel_base_in": "in", "length_in": "in", "width_in": "in",
                 "height_in": "in", "curb_weight_lb": "lb",
                 "engine_size_cc": "cc", "bore": "in", "stroke": "in",
                 "compression_ratio": "-", "horsepower": "hp",
                 "peak_rpm": "rpm", "city_mpg": "MPG (target)",
                 "hwy_mpg": "MPG", "price_usd": "USD (target)"},
        relations=[
            rel("price_usd", ["curb_weight_lb"], "价格~整备质量"),
            rel("price_usd", ["engine_size_cc"], "价格~排量"),
            rel("city_mpg", ["curb_weight_lb"],
                "城市油耗~车重: 负指数(滚阻+惯量)"),
            rel("horsepower", ["engine_size_cc"], "马力~排量: 近比例"),
        ],
        raw_files=[path],
        notes=["'?' 缺失读为 NaN(原格式), 运行时按关系 dropna 不插值",
               "make/车身形式等分类列清洗期剔除"])


def main():
    os.makedirs(CLEAN, exist_ok=True)
    prep_auto_mpg()
    prep_concrete()
    prep_concrete_slump()
    prep_energy_efficiency()
    prep_yacht()
    prep_ccpp()
    prep_air_quality()
    prep_tetouan()
    prep_steel_energy()
    prep_gas_turbine()
    prep_nist_strd()
    prep_webbook()
    prep_usgs()
    prep_rdatasets()
    prep_open_meteo()
    prep_usgs_eq()
    prep_openml()
    prep_epa_fueleconomy()
    prep_automobile()
    out = os.path.join(CLEAN, "manifest.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(ENTRIES, f, ensure_ascii=False, indent=2)
    n_rel = sum(len(e["candidate_relations"]) for e in ENTRIES.values())
    n_rows = sum(e["n_rows"] for e in ENTRIES.values())
    print(f"\ndone: {len(ENTRIES)} datasets, {n_rel} candidate relations, "
          f"{n_rows} rows -> {os.path.relpath(out, HERE)}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
