"""wb_blindspot.py — 真实数据上的 fc 盲区探针(交互/稀释方向)。

机理(scenarios.py):fc 的检验统计量 = 唯一一列"平均中心化二次项"的系数,
对三类失配天然低敏感:
  INTERACT — lx1*lx2 交叉项(与二次列无关)
  DILUTION — 双特征只有一个有曲率,平均列信号减半
  OPPOSING — 两特征曲率反向,平均列信号精确抵消
本探针在世行农业/能源数据上找同方向的"RESET 拒绝而 fc 沉默"的分歧格:
  1) yield ~ fert                (单特征对照:曲率应被双方看到)
  2) yield ~ fert + energy_pc    (双特征:交互/稀释疑点)
  3) ghg_pc ~ gdp_pc + energy_pc (双特征,能源强度路径)
真实数据无真值,只报告分歧及其与合成盲区的方向一致性,不作功效结论。
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import requests

import common

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "runs", "wb_real")
RAW = os.path.join(OUT, "raw")
YEARS = ["2018", "2019", "2022"]
SEEDS = [11, 22, 33, 44, 55]

IND = {
    "AG.YLD.CREL.KG": "yield",           # 谷物单产 kg/ha
    "AG.CON.FERT.ZS": "fert",            # 施肥量(可得的口径)
    "EG.USE.PCAP.KG.OE": "energy_pc",    # 人均能耗 kg油当量
    "NY.GDP.PCAP.PP.CD": "gdp_pc_ppp",
    "EN.GHG.ALL.MT.CE.AR5": "ghg_mt",
    "SP.POP.TOTL": "pop",
}
REL = [
    ("yield__fert", "yield", ["fert"], "单特征对照"),
    ("yield__fert_energy", "yield", ["fert", "energy_pc"], "双特征:交互/稀释疑点"),
    ("ghgpc__gdppc_energy", "ghg_pc", ["gdp_pc_ppp", "energy_pc"], "双特征:能源强度路径"),
]


def _get(url, cache):
    path = os.path.join(RAW, cache)
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    last = None
    for a in range(6):
        try:
            r = requests.get(url, timeout=60)
            r.raise_for_status()
            data = r.json()
            with open(path, "w", encoding="utf-8") as f:
                json.dump(data, f, ensure_ascii=False)
            return data
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2 * (a + 1))
    raise RuntimeError(f"{url}: {last!r}")


def main():
    os.makedirs(RAW, exist_ok=True)
    meta = _get("https://api.worldbank.org/v2/country?format=json&per_page=400",
                "countries.json")
    real = {c["id"] for c in meta[1] if c["region"]["value"] != "Aggregates"}
    rows = []
    for year in YEARS:
        cols = {}
        for code, name in IND.items():
            d = _get(f"https://api.worldbank.org/v2/country/all/indicator/{code}"
                     f"?date={year}&format=json&per_page=20000",
                     f"{code}_{year}.json")
            if len(d) < 2 or not d[1]:
                print(f"  ! 无数据: {code} {year}")
                continue
            cols[name] = {r["countryiso3code"]: float(r["value"])
                          for r in d[1] if r["value"] is not None}
        df = pd.DataFrame({k: pd.Series(v) for k, v in cols.items()})
        df = df[df.index.isin(real)].reset_index().rename(columns={"index": "iso3"})
        if "ghg_mt" in df and "pop" in df:
            df["ghg_pc"] = df["ghg_mt"] * 1e6 / df["pop"]
        for task, tgt, feats, desc in REL:
            sub = df[["iso3", tgt] + feats].dropna()
            sub = sub[(sub[tgt] > 0) & (sub[feats] > 0).all(axis=1)]
            n = len(sub)
            if n < 40:
                print(f"{year} {task}: n={n} 不足跳过")
                continue
            for seed in SEEDS:
                row = {"year": int(year), "task": task, "desc": desc, "n": n,
                       "seed": seed}
                row.update(common.run_all_methods(sub, tgt, feats, seed=seed))
                rows.append(row)
            r0 = next(r for r in rows if r["task"] == task and r["year"] == int(year))
            print(f"{year} {task}: n={n} p_fc={r0['p_fc']:.4f} "
                  f"pF={r0['p_fitted_f']:.2e} pHC3={r0['p_fitted_hc3']:.2e} "
                  f"pEX={r0['p_exog_hc3']:.2e}")
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT, "blindspot.csv"), index=False, encoding="utf-8-sig")
    med = res.groupby(["year", "task"])[
        ["p_fc", "p_fitted_f", "p_fitted_hc3", "p_exog_hc3"]].median()
    print("\n=== 中位 p ===")
    print(med.to_string())


if __name__ == "__main__":
    sys.exit(main())
