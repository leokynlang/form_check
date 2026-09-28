"""fetch_data.py — 为 form_check 真实数据测试批量抓取免鉴权工程数据源。

全部免鉴权、可复现(原始字节落盘 + sha256):
  1. UCI ML Repository (archive.ics.uci.edu)      — 大样本工程回归台架
  2. NIST StRD 非线性回归 (itl.nist.gov)          — 真实测量的标定/物理曲线
  3. NIST Chemistry WebBook (webbook.nist.gov)    — 流体物性 vs 温度(已知非幂律形状)
     3b. WebBook Mask=4 相变页扩容(~40 流体, 2026-09-06 二期)
  4. USGS NWIS 瞬时值 (waterservices.usgs.gov)   — 河流率定曲线 Q-h(教科书幂律)
     4b. 大河测站扩容(10 站, 2026-09-06 二期)
  5. RDatasets 镜像 (vincentarelbundock.github.io)— 经典物理/工程小样本
  6. Open-Meteo 历史存档 (archive-api.open-meteo.com, ERA5 再分析) — 风速切变/阵风
  7. USGS FDSN 地震目录 (earthquake.usgs.gov)     — 震级 vs 有感报告/台站数
  8. OpenML (www.openml.org)                      — F-16 飞控/机器人动力学/系统遥测
  9. EPA fueleconomy.gov                          — 全美车辆工程与油耗数据库(~48k 行)

用法: python fetch_data.py            (幂等: 已存在的 raw 文件跳过)
写入: real_data/raw/<source>/...      (仅原始数据; 清洗在 prepare_data.py)
"""
from __future__ import annotations

import hashlib
import io
import json
import os
import sys
import time
import zipfile

import requests

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw")

UA = {"User-Agent": "Mozilla/5.0 (form_check data fetch; research use)"}


def fetch_bytes(url, dest, tries=4, timeout=90, headers=None):
    """下载 url 到 dest(字节级缓存: 已存在则跳过)。返回 (bytes, cached)。"""
    if os.path.exists(dest):
        with open(dest, "rb") as f:
            return f.read(), True
    hdr = headers or UA
    last = None
    for attempt in range(tries):
        try:
            r = requests.get(url, headers=hdr, timeout=timeout)
            r.raise_for_status()
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            with open(dest, "wb") as f:
                f.write(r.content)
            return r.content, False
        except Exception as exc:  # noqa: BLE001 — 网络抖动重试
            last = exc
            time.sleep(2.0 * (attempt + 1))
    raise RuntimeError(f"fetch failed: {url}: {last!r}")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------- 1. UCI
UCI_IDS = {
    9:   "auto_mpg",             # 汽车: mpg vs 排量/马力/车重 (398 行)
    165: "concrete_strength",    # 混凝土: 抗压强度 vs 配比+龄期 (1030 行)
    182: "concrete_slump",       # 混凝土坍落度试验: 坍落/流动/强度 vs 配比 (103 行)
    242: "energy_efficiency",    # 建筑形状参数 → 采暖/制冷负荷 (768 行)
    243: "yacht_hydro",          # 船体水动力: 剩余阻力 vs Froude 数 (308 行)
    294: "ccpp",                 # 联合循环电厂: 净电力 vs 环境变量 (9568 行)
    360: "air_quality",          # 意大利气体传感器阵列 vs 参考分析仪 (9358 小时)
    551: "gas_turbine_emissions",# 燃机 CO/NOx vs 环境+燃烧参数 (~36k 行/年)
    849: "tetouan_power",        # 摩洛哥 Tetouan 三区电力消费 vs 气象 (52k 行)
    851: "steel_energy",         # 韩国小型钢厂能耗遥测 (35k 行, 15min)
    994: "micro_gas_turbine",    # 微型燃机两年 10min 遥测 (~7k 行)
    10:  "automobile",           # 1985 进口轿车: 价格/马力/整备质量 (205 行)
}


def fetch_uci():
    out = {}
    api = "https://archive.ics.uci.edu/api/dataset?id={id}"
    static = "https://archive.ics.uci.edu/static/public/{id}/{file}"
    for uid, short in UCI_IDS.items():
        meta_path = os.path.join(RAW, "uci", f"{uid}_meta.json")
        try:
            text, _ = fetch_bytes(api.format(id=uid), meta_path)
            data = json.loads(text)["data"]
        except Exception as exc:  # noqa: BLE001
            print(f"  ! UCI {uid} meta failed: {exc}")
            continue
        dest_dir = os.path.join(RAW, "uci", f"{uid}_{short}")
        os.makedirs(dest_dir, exist_ok=True)
        got = []
        try:
            if data.get("data_url"):
                b, _ = fetch_bytes(data["data_url"],
                                   os.path.join(dest_dir, "data.csv"))
                got.append(("data.csv", b))
        except Exception:
            pass
        if not got:
            # 兜底: 名字 slug 的 zip(Yacht/微型燃机走此路径)
            slug = data["repository_url"].rstrip("/").split("/")[-1]
            try:
                b, _ = fetch_bytes(static.format(id=uid, file=f"{slug}.zip"),
                                   os.path.join(dest_dir, "bundle.zip"))
                zf = zipfile.ZipFile(io.BytesIO(b))
                zf.extractall(dest_dir)
                got.append(("bundle.zip", b))
            except Exception as exc:  # noqa: BLE001
                print(f"  ! UCI {uid} download failed: {exc}")
                continue
        out[short] = {"uci_id": uid, "name": data["name"],
                      "abstract": data.get("abstract", ""),
                      "n_instances": data.get("num_instances"),
                      "variables": [
                          {"name": v.get("name"), "role": v.get("role"),
                           "units": v.get("units"), "type": v.get("type")}
                          for v in (data.get("variables") or [])],
                      "dir": os.path.relpath(dest_dir, HERE)}
        print(f"  ok UCI {uid} {short}: {data['name']} "
              f"({data.get('num_instances')} inst)")
    return out


# ------------------------------------------------- 2. NIST StRD 非线性回归
# 真实测量工程数据(超声标定/热膨胀/光学干涉/能量效率...), .dat 头部含说明。
NIST_NLS = ["Chwirut1", "Chwirut2", "Misra1a", "Misra1b", "Misra1c",
            "Misra1d", "DanWood", "Kirby2", "Hahn1", "Roszman1",
            "Thurber", "Eckerle4", "Bennett5"]


def fetch_nist_strd():
    base = ("https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/"
            "{name}.dat")
    out = {}
    for name in NIST_NLS:
        dest = os.path.join(RAW, "nist_strd", f"{name}.dat")
        try:
            fetch_bytes(base.format(name=name), dest)
            out[name] = sha256(dest)
        except Exception as exc:  # noqa: BLE001
            print(f"  ! NIST {name}: {exc}")
            continue
        print(f"  ok NIST {name}")
    return out


# ------------------------------------------- 3. NIST WebBook 流体物性
# Mask=1024 是 "Fluid Properties"(粘度/热导率 vs 温度表), HTML 表格。
# C7732185=水, C124389=CO2, C7727379=氮气, C6418=甲醇, C67175=乙醇。
WEBBOOK = {
    "water":      "C7732185",
    "co2":        "C124389",
    "nitrogen":   "C7727379",
    "methanol":   "C6418",
    "ethanol":    "C67175",
}
WEBBOOK_URL = ("https://webbook.nist.gov/cgi/cbook.cgi"
               "?ID={cid}&Units=SI&Mask=1024")


def fetch_webbook():
    out = {}
    for short, cid in WEBBOOK.items():
        dest = os.path.join(RAW, "webbook", f"{short}_fluid.html")
        try:
            fetch_bytes(WEBBOOK_URL.format(cid=cid), dest)
            out[short] = {"cid": cid, "sha256": sha256(dest)}
            print(f"  ok WebBook {short} ({cid})")
            time.sleep(0.4)
        except Exception as exc:  # noqa: BLE001
            print(f"  ! WebBook {short}: {exc}")
    return out


# 3b. 相变页扩容(Mask=4): ΔvapH/ΔsubH vs T 逐点表, 小 n 已知非幂律考场。
#     CID = "C" + CAS 号去连字符。含一期 5 流体的相变补抓(其 Mask=1024
#     流体物性页不含相变表)。
WEBBOOK_PHASE = {
    "methane": "C74828",        "ethane": "C74840",
    "propane": "C74986",        "butane": "C106978",
    "pentane": "C109660",       "nonane": "C111842",
    "decane": "C124184",        "isooctane": "C540841",
    "ethylene": "C74851",       "propylene": "C115071",
    "o_xylene": "C95476",       "ethylbenzene": "C100414",
    "1_propanol": "C71238",     "1_butanol": "C71363",
    "ethylene_glycol": "C107211", "aniline": "C62533",
    "phenol": "C108952",        "ammonia": "C7664417",
    "chlorine": "C7782505",     "sulfur_dioxide": "C7446095",
    "hydrogen_sulfide": "C7783064", "oxygen": "C7782447",
    "carbon_monoxide": "C630080", "carbon_disulfide": "C75150",
    "carbon_tetrachloride": "C56235", "chloroform": "C67663",
    "methylene_chloride": "C75092", "acetylene": "C74862",
    "methyl_chloride": "C74873", "hydrogen": "C1333740",
    "helium": "C7440597",       "neon": "C7440019",
    "nitrous_oxide": "C10024972", "r134a": "C811972",
    # 一期 5 流体的相变页补抓
    "water": "C7732185",        "co2": "C124389",
    "nitrogen": "C7727379",     "methanol": "C67561",
    "ethanol": "C67175",
}
WEBBOOK_PHASE_URL = ("https://webbook.nist.gov/cgi/cbook.cgi"
                     "?ID={cid}&Units=SI&Mask=4")


def fetch_webbook_phase():
    out = {}
    for short, cid in WEBBOOK_PHASE.items():
        dest = os.path.join(RAW, "webbook", f"{short}_phase.html")
        try:
            fetch_bytes(WEBBOOK_PHASE_URL.format(cid=cid), dest)
            out[short] = {"cid": cid, "sha256": sha256(dest)}
            print(f"  ok WebBook-phase {short} ({cid})")
            time.sleep(0.4)
        except Exception as exc:  # noqa: BLE001
            print(f"  ! WebBook-phase {short}: {exc}")
    return out


# ------------------------------------------------------- 4. USGS 日值
# Q(00060 流量) + h(00065 水位) 日值, 教科书率定曲线 Q = C(h-h0)^n。
USGS_SITES = {
    "01646500": ("potomac",   "2010-01-01", "2024-12-31"),
    "12484500": ("nooksack",  "2010-01-01", "2024-12-31"),
    "08158000": ("guadalupe", "2010-01-01", "2024-12-31"),
}
USGS_URL = ("https://waterservices.usgs.gov/nwis/dv/"
            "?sites={site}&startDT={start}&endDT={end}"
            "&parameterCd=00060,00065&statCd=00003&format=json")


def fetch_usgs():
    out = {}
    for site, (short, start, end) in USGS_SITES.items():
        dest = os.path.join(RAW, "usgs", f"{short}_{site}.json")
        try:
            fetch_bytes(USGS_URL.format(site=site, start=start, end=end),
                        dest)
            out[short] = {"site": site, "start": start, "end": end,
                          "sha256": sha256(dest)}
            print(f"  ok USGS {short} ({site})")
        except Exception as exc:  # noqa: BLE001
            print(f"  ! USGS {site}: {exc}")
    return out


# 率定曲线用瞬时值(00065 水位无日值发布, 走 iv, Q+h 按时间戳配对)。
USGS_IV = {
    "01646500": "potomac",
    "12484500": "nooksack",
    "08158000": "guadalupe",
}
USGS_IV_URL = ("https://waterservices.usgs.gov/nwis/iv/"
               "?sites={site}&startDT={start}&endDT={end}"
               "&parameterCd=00060,00065&format=json")


def fetch_usgs_iv():
    out = {}
    for site, short in USGS_IV.items():
        dest = os.path.join(RAW, "usgs", f"iv_{short}_{site}_2023.json")
        try:
            fetch_bytes(USGS_IV_URL.format(site=site,
                                           start="2023-01-01",
                                           end="2023-12-31"),
                        dest, tries=3, timeout=300)
            out[short] = {"site": site, "start": "2023-01-01",
                          "end": "2023-12-31", "sha256": sha256(dest)}
            print(f"  ok USGS-iv {short} ({site})")
        except Exception as exc:  # noqa: BLE001
            print(f"  ! USGS-iv {site}: {exc}")
    return out


# 4b. 大河测站扩容(2026-09-06 二期): 10 站覆盖密苏里/俄亥俄/科罗拉多/
#     哥伦比亚/萨斯奎哈纳/特拉华/密西西比/萨克拉门托/阿查法拉亚/康涅狄格。
#     部分站可能无 00065 瞬时序列(如 07289000 Vicksburg), 清洗期校验并跳过。
USGS_IV2 = {
    "06934500": "missouri_hermann",
    "03255000": "ohio_cincinnati",
    "09380000": "colorado_leesferry",
    "14105700": "columbia_thedalles",
    "01570500": "susquehanna_harrisburg",
    "01463500": "delaware_trenton",
    "07010000": "mississippi_memphis",
    "11447650": "sacramento_freeport",
    "07374000": "atchafalaya_simmesport",
    "01184000": "connecticut_thompsonville",
}


def fetch_usgs_iv2():
    out = {}
    for site, short in USGS_IV2.items():
        dest = os.path.join(RAW, "usgs", f"iv_{short}_{site}_2023.json")
        try:
            fetch_bytes(USGS_IV_URL.format(site=site,
                                           start="2023-01-01",
                                           end="2023-12-31"),
                        dest, tries=3, timeout=300)
            out[short] = {"site": site, "start": "2023-01-01",
                          "end": "2023-12-31", "sha256": sha256(dest)}
            print(f"  ok USGS-iv2 {short} ({site})")
            time.sleep(0.5)
        except Exception as exc:  # noqa: BLE001
            print(f"  ! USGS-iv2 {site}: {exc}")
    return out


# ------------------------------------------- 6. Open-Meteo 历史存档(ERA5)
# 免 key 免费档用于非商业研究。风速切变 v100 = v10·(100/10)^α 是风工程
# 教科书幂律(α≈0.11 海上-0.2 陆上); 阵风系数 gust/v10 近常数。
# 注意: ERA5 是再分析(观测同化模型), 非直接测风, manifest 如实标注。
OPEN_METEO_SITES = {
    "north_sea":      (56.5, 3.0),    # 北海海上风场区
    "german_bight":   (54.0, 7.5),    # 德国湾海上
    "baltic_ruegen":  (54.6, 13.4),   # 波罗的海吕根岛外海
    "kansas_plains":  (38.5, -98.5),  # 大平原风电带
    "texas_panhandle": (35.2, -101.8),# 德州潘汉德尔风电带
    "altamont_pass":  (37.73, -121.66),# 加州 Altamont 山口老风场
    "iowa_fields":    (42.0, -93.5),  # 艾奥瓦平原
    "inner_mongolia": (44.0, 113.5),  # 内蒙古风电基地
    "qingdao_offshore": (35.9, 120.8),# 黄海近海
    "patagonia":      (-47.0, -68.5), # 巴塔哥尼亚(全球风资源极值区)
}
OPEN_METEO_URL = ("https://archive-api.open-meteo.com/v1/archive"
                  "?latitude={lat}&longitude={lon}"
                  "&start_date=2023-01-01&end_date=2023-12-31"
                  "&hourly=wind_speed_10m,wind_speed_100m,wind_gusts_10m"
                  "&timezone=UTC")


def fetch_open_meteo():
    out = {}
    for short, (lat, lon) in OPEN_METEO_SITES.items():
        dest = os.path.join(RAW, "open_meteo", f"{short}_2023.json")
        try:
            fetch_bytes(OPEN_METEO_URL.format(lat=lat, lon=lon), dest)
            out[short] = {"lat": lat, "lon": lon, "sha256": sha256(dest)}
            print(f"  ok Open-Meteo {short}")
            time.sleep(0.5)
        except Exception as exc:  # noqa: BLE001
            print(f"  ! Open-Meteo {short}: {exc}")
    return out


# ------------------------------------------- 7. USGS FDSN 地震目录
# M≥5.5 全球事件 2000-2024 (count=12336)。GeoJSON 才带 felt/cdi/mmi。
USGS_EQ_URL = ("https://earthquake.usgs.gov/fdsnws/event/1/query"
               "?format=geojson&starttime=2000-01-01&endtime=2024-12-31"
               "&minmagnitude=5.5&orderby=time-asc")


def fetch_usgs_eq():
    dest = os.path.join(RAW, "usgs_eq", "global_m55_2000_2024.geojson")
    try:
        fetch_bytes(USGS_EQ_URL, dest, timeout=300)
        n = json.loads(open(dest, encoding="utf-8").read())["metadata"]["count"]
        print(f"  ok USGS-EQ global M>=5.5 2000-2024: {n} events")
        return {"global_m55_2000_2024": {"url": USGS_EQ_URL,
                                         "n_events": n,
                                         "sha256": sha256(dest)}}
    except Exception as exc:  # noqa: BLE001
        print(f"  ! USGS-EQ: {exc}")
        return {}


# ---------------------------------------------------- 5. RDatasets 经典
# 列: 数据集(包/名), 描述。
RDATASETS = {
    "pressure": ("datasets/pressure.csv",
                 "汞蒸气压 vs 温度(19 行, Weast 手册)"),
    "cars":     ("datasets/cars.csv",
                 "1920s 汽车速度 vs 制动距离(50 行)"),
    "rock":     ("datasets/rock.csv",
                 "石油工程: 岩芯渗透率 vs 面积/周长(48 行)"),
    "trees":    ("datasets/trees.csv",
                 "黑樱桃树: 材积 vs 胸径/树高(31 行, 材积~直径^2)"),
    "attenu":   ("datasets/attenu.csv",
                 "加州大地震工程: 峰值加速度 vs 震级/震中距(182 行)"),
    "dnipro":   ("Ecdat/Dnidaho.csv", None),  # 占位, 下面会过滤
}
RDATASETS = {k: v for k, v in RDATASETS.items() if v[1] is not None}
RDATASETS["housing"] = ("Ecdat/Housing.csv",
                        "房屋: 销售价 vs 地块/结构属性(546 行)")
RDATASETS["electricity"] = ("Ecdat/Electricity.csv",
                            "电力公司: 边际成本 vs 产出/燃料(158 行)")
# 二期(2026-09-06): 工程/物理经典补充
RDATASETS["stormer"] = ("MASS/stormer.csv",
                        "Stormer 旋转粘度计: 流出时间 vs 粘度/砝码重(23 行)")
RDATASETS["whiteside"] = ("MASS/whiteside.csv",
                          "住宅保温前后: 燃气消耗 vs 室外温度(56 行)")
RDATASETS["cement_hald"] = ("MASS/cement.csv",
                            "水泥熟料四组分 vs 水化热(Hald 数据, 13 行)")
RDATASETS["petrol"] = ("MASS/petrol.csv",
                        "Prater 炼厂: 汽油收率 vs 原油物性(32 行)")
RDATASETS["women"] = ("datasets/women.csv",
                      "美国女性平均身高体重(15 行): 体重~身高^p 异速律")
RDATASETS["airquality"] = ("datasets/airquality.csv",
                           "纽约空气质量: 臭氧 vs 风/温/辐射(153 行)")
RDATASETS["puromycin"] = ("datasets/Puromycin.csv",
                          "酶促反应速率 vs 底物浓度(23 行): Michaelis-Menten")


def fetch_rdatasets():
    base = ("https://vincentarelbundock.github.io/Rdatasets/csv/{path}")
    out = {}
    for short, (path, desc) in RDATASETS.items():
        dest = os.path.join(RAW, "rdatasets", f"{short}.csv")
        try:
            fetch_bytes(base.format(path=path), dest)
            out[short] = {"path": path, "desc": desc,
                          "sha256": sha256(dest)}
            print(f"  ok RDatasets {short}")
        except Exception as exc:  # noqa: BLE001
            print(f"  ! RDatasets {short}: {exc}")
    return out


# ------------------------------------------- 8. OpenML 航空/机器人/系统
# F-16 飞控系(Jaroszewicz: elevators/ailerons/delta_*)、PUMA 机械臂动力学、
# 服务器系统遥测、经典小样本(elusage/cloud)。ARFF 静态直链(openml.org 对
# 浏览器式 UA 返回 412, 用 curl 风格 UA 实测 200)。
OPENML_SETS = {
    "elevators": (216, "https://openml.org/data/v1/download/3653/elevators.arff"),
    "ailerons": (296, "https://openml.org/data/v1/download/52060/Ailerons.arff"),
    "delta_elevators": (198, "https://openml.org/data/v1/download/52348/delta_elevators.arff"),
    "delta_ailerons": (803, "https://openml.org/data/v1/download/53337/delta_ailerons.arff"),
    "puma8NH": (225, "https://openml.org/data/v1/download/3662/puma8NH.arff"),
    "cpu_act": (197, "https://openml.org/data/v1/download/3634/cpu_act.arff"),
    "elusage": (228, "https://openml.org/data/v1/download/3665/elusage.arff"),
    "cloud": (210, "https://openml.org/data/v1/download/3647/cloud.arff"),
    "pol": (201, "https://openml.org/data/v1/download/3638/pol.arff"),
}
OPENML_UA = {"User-Agent": "curl/8.0"}  # openml 反爬: Mozilla 式 UA 被 412


def fetch_openml():
    out = {}
    for short, (did, url) in OPENML_SETS.items():
        dest = os.path.join(RAW, "openml", f"{short}.arff")
        try:
            fetch_bytes(url, dest, headers=OPENML_UA)
            out[short] = {"did": did, "url": url, "sha256": sha256(dest)}
            print(f"  ok OpenML {short} (did={did})")
        except Exception as exc:  # noqa: BLE001
            print(f"  ! OpenML {short}: {exc}")
    return out


# ------------------------------------------- 9. EPA fueleconomy.gov
# 全美车辆数据库(1984-至今, ~48k 行): 排量/缸数/油耗/尾气 CO2。
# 官方年度快照 zip, 免鉴权。
EPA_URL = ("https://www.fueleconomy.gov/feg/epadata/vehicles.csv.zip")


def fetch_epa():
    dest = os.path.join(RAW, "epa_fueleconomy", "vehicles.csv.zip")
    try:
        fetch_bytes(EPA_URL, dest)
        zf = zipfile.ZipFile(dest)
        names = zf.namelist()
        print(f"  ok EPA fueleconomy: {names} "
              f"({os.path.getsize(dest) // 1024} KB)")
        return {"vehicles": {"url": EPA_URL, "zip_entries": names,
                             "sha256": sha256(dest)}}
    except Exception as exc:  # noqa: BLE001
        print(f"  ! EPA fueleconomy: {exc}")
        return {}


def main():
    os.makedirs(RAW, exist_ok=True)
    manifest = {}
    print("== 1/9 UCI ==")
    manifest["uci"] = fetch_uci()
    print("== 2/9 NIST StRD ==")
    manifest["nist_strd"] = fetch_nist_strd()
    print("== 3/9 NIST WebBook ==")
    manifest["webbook"] = fetch_webbook()
    manifest["webbook_phase"] = fetch_webbook_phase()
    print("== 4/9 USGS ==")
    manifest["usgs"] = fetch_usgs()
    manifest["usgs_iv"] = fetch_usgs_iv()
    manifest["usgs_iv2"] = fetch_usgs_iv2()
    print("== 5/9 RDatasets ==")
    manifest["rdatasets"] = fetch_rdatasets()
    print("== 6/9 Open-Meteo (ERA5) ==")
    manifest["open_meteo"] = fetch_open_meteo()
    print("== 7/9 USGS Earthquakes ==")
    manifest["usgs_eq"] = fetch_usgs_eq()
    print("== 8/9 OpenML ==")
    manifest["openml"] = fetch_openml()
    print("== 9/9 EPA fueleconomy ==")
    manifest["epa_fueleconomy"] = fetch_epa()

    manifest_path = os.path.join(RAW, "fetch_manifest.json")
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, ensure_ascii=False, indent=2)
    n = sum(len(v) for v in manifest.values())
    print(f"done: {n} 项 -> {os.path.relpath(manifest_path, HERE)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
