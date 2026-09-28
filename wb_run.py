"""wb_run.py — 世界银行真实宏观数据上的 form_check vs RESET 对比。

数据源: World Bank Open Data API v2 (免鉴权, https://api.worldbank.org)。
联合国基金会(UN Foundation)本身无公开数据 API;世界银行是联合国系统相关
机构中数据开放最完整的来源,故以世行数据作为"联合国体系"真实数据代表。

设计: 跨国横截面(剔除世行聚合体, 只留主权经济体), 若干候选幂律关系,
每个关系 × 多个年份, 每格 5 个排列种子看经验 p 稳定性。
没有合成真值, 故本实验只报告四种方法的判定(p 值/拒绝/失败/条件数),
不做功效/FPR 结论 — 那需要已知 H0/H1, 真实数据不提供。

用法: python wb_run.py   (幂等: 已缓存 JSON 不重拉)
写入: runs/wb_real/ (raw 缓存 / results.csv / report.md)
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
YEARS = ["2019", "2021", "2022", "2023"]      # 跳过 2020 疫情年异常,保留前后对照
SEEDS = [11, 22, 33, 44, 55]
ALPHA = common.ALPHA

INDICATORS = {
    "NY.GDP.MKTP.PP.CD": "gdp_ppp",        # GDP, PPP (现价国际元)
    "NY.GDP.PCAP.PP.CD": "gdp_pc_ppp",     # 人均 GDP, PPP
    "SP.POP.TOTL": "pop",                  # 总人口
    "SP.DYN.LE00.IN": "life_exp",          # 出生预期寿命
    "SP.DYN.IMRT.IN": "inf_mort",          # 婴儿死亡率(每千活产)
    "EN.GHG.CO2.MT.CE.AR5": "co2_mt",      # CO2 排放(百万吨, AR5 口径)
    "SH.XPD.CHEX.PC.CD": "health_exp_pc",  # 人均卫生支出
    "IT.NET.USER.ZS": "net_users_pct",     # 互联网使用率 %
}

# 预注册关系: (task, target, features, 说明)
RELATIONS = [
    ("lifeexp__gdppc", "life_exp", ["gdp_pc_ppp"],
     "预期寿命 ~ 人均GDP(经典凹关系, 预期非幂律)"),
    ("infmort__gdppc", "inf_mort", ["gdp_pc_ppp"],
     "婴儿死亡率 ~ 人均GDP( Preston 曲线另一支)"),
    ("gdp__pop", "gdp_ppp", ["pop"],
     "总量GDP ~ 总人口(规模关系)"),
    ("co2pc__gdppc", "co2_pc", ["gdp_pc_ppp"],
     "人均CO2 ~ 人均GDP(环境库兹涅茨类)"),
    ("netusers__gdppc", "net_users_pct", ["gdp_pc_ppp"],
     "互联网普及率 ~ 人均GDP(S 型中段)"),
    ("lifeexp__gdppc_health", "life_exp", ["gdp_pc_ppp", "health_exp_pc"],
     "预期寿命 ~ 人均GDP + 人均卫生支出(双特征)"),
]


def _get(url, cache_name):
    os.makedirs(RAW, exist_ok=True)
    path = os.path.join(RAW, cache_name)
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    last = None
    for attempt in range(6):
        try:
            r = requests.get(url, timeout=60)
            r.raise_for_status()
            data = r.json()
            break
        except Exception as exc:  # noqa: BLE001 — 网络抖动重试
            last = exc
            time.sleep(2.0 * (attempt + 1))
    else:
        raise RuntimeError(f"fetch failed after retries: {url}: {last!r}")
    # 世行 API 返回 [meta, payload]
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False)
    time.sleep(0.3)
    return data


def load_country_set():
    """主权经济体白名单: 剔除 'Aggregates' 聚合体。"""
    data = _get("https://api.worldbank.org/v2/country?format=json&per_page=400",
                "countries.json")
    rows = data[1]
    real = {c["id"] for c in rows
            if c["region"]["value"] != "Aggregates"}
    return real


def load_indicator(code, year):
    data = _get(f"https://api.worldbank.org/v2/country/all/indicator/{code}"
                f"?date={year}&format=json&per_page=20000",
                f"{code}_{year}.json")
    out = {}
    if len(data) < 2 or not data[1]:
        print(f"  ! 指标无数据/已归档: {code} {year}")
        return out
    for rec in data[1]:
        if rec["value"] is not None:
            out[rec["countryiso3code"]] = float(rec["value"])
    return out


def build_dataset(year, real_countries):
    cols = {}
    for code, name in INDICATORS.items():
        cols[name] = load_indicator(code, year)
    df = pd.DataFrame({k: pd.Series(v) for k, v in cols.items()})
    df.index.name = "iso3"
    df = df[df.index.isin(real_countries)]
    # 派生: 人均 CO2(千吨/人 → 吨/人)
    if "co2_mt" in df and "pop" in df:
        df["co2_pc"] = df["co2_mt"] * 1e6 / df["pop"]
    return df.reset_index()


def main():
    os.makedirs(OUT, exist_ok=True)
    real = load_country_set()
    print(f"主权经济体: {len(real)} 个")

    results = []
    for year in YEARS:
        df = build_dataset(year, real)
        for task, target, features, desc in RELATIONS:
            sub = df[["iso3", target] + features].dropna()
            sub = sub[(sub[target] > 0) & (sub[features] > 0).all(axis=1)]
            n = len(sub)
            base = {"year": int(year), "task": task, "desc": desc,
                    "target": target, "features": ",".join(features), "n": n}
            if n < max(30, 5 * (len(features) + 3)):
                base.update({"status": "insufficient_n"})
                results.append(base)
                continue
            for seed in SEEDS:
                row = dict(base)
                row["seed"] = seed
                row["status"] = "ok"
                try:
                    row.update(common.run_all_methods(sub, target, features,
                                                      seed=seed))
                except Exception as exc:  # noqa: BLE001 — 真实数据留痕不抛
                    row["status"] = f"runner_error: {exc!r}"
                results.append(row)
            print(f"{year} {task}: n={n}")

    res = pd.DataFrame(results)
    res.to_csv(os.path.join(OUT, "results.csv"), index=False, encoding="utf-8-sig")

    # ---- 汇总: 每格取种子中位 p 与一致拒绝率 ----
    ok = res[res["status"] == "ok"].copy()
    method_keys = ["fc", "fitted_f", "fitted_hc3", "exog_hc3"]
    lines = ["# form_check vs RESET — 世界银行真实宏观数据", "",
             f"- 数据: World Bank Open Data API, 跨国横截面, 年份 {YEARS}, "
             f"剔除聚合体, 每格 {len(SEEDS)} 个排列种子(form_check 经验 p 的稳定性)",
             f"- 统一 α={ALPHA}; 方法 A=fc(排列经验 p), B=RESET fitted-F, "
             "C=RESET fitted-HC3, D=RESET exog-HC3", "",
             "| 年份 | 关系 | n | fc p (中位) | fitted-F p | fitted-HC3 p | exog-HC3 p | 一致结论 |",
             "|---|---|---|---|---|---|---|---|"]
    for (year, task), g in ok.groupby(["year", "task"]):
        med = {k: g[f"p_{k}"].median() for k in method_keys}
        rejs = {k: g[f"rej_{k}"].mean() for k in method_keys}
        fails = [f"{k}:{g[f'fail_{k}'].notna().sum()}" for k in method_keys
                 if g[f"fail_{k}"].notna().any()]
        verdict = " / ".join(
            f"{ 'REJ' if rejs[k] >= 0.5 else 'no-rej' }" for k in method_keys)
        ps = " / ".join(
            "—" if not np.isfinite(med[k]) else f"{med[k]:.4f}" for k in method_keys)
        extra = f" fail[{'; '.join(fails)}]" if fails else ""
        lines.append(f"| {year} | {task} | {int(g['n'].iloc[0])} | "
                     f"{ps.split(' / ')[0]} | {ps.split(' / ')[1]} | "
                     f"{ps.split(' / ')[2]} | {ps.split(' / ')[3]} | "
                     f"{verdict}{extra} |")
    lines += ["", "## 备注", "- 真实数据无已知真值, 本表只呈现判定差异, 不含功效/校准结论。",
              "- p_fc 为排列经验 p(5 种子中位数), RESET 三个 p 值为确定性解析值。"]
    with open(os.path.join(OUT, "report.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("done ->", OUT)
    print(ok.groupby(["year", "task"])[[f"p_{k}" for k in method_keys]]
          .median().to_string())


if __name__ == "__main__":
    sys.exit(main())
