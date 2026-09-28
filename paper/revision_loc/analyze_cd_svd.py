# -*- coding: utf-8 -*-
"""analyze_cd_svd.py — v5 重跑汇总:C/D Gram(存档) vs SVD(新) 逐套件对比。

输出 runs/cd_svd/summary_cd_svd.json + 控制台报告:
  - S1/S2 H0 尺寸(C/D,按 n,合并噪声层;旧 vs 新 + Wilson 区间)
  - S1 结构 H1 功效(n=200, rel=1%,逐场景,旧 vs 新)
  - S2 H1 功效(n=200, 合并 rel≥0.5%,总整体;旧 vs 新)
  - 哨兵统计与失败类型变化
"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = r"D:\form_check"
CD = os.path.join(HERE, "runs", "cd_svd")


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"),) * 2
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)


def rate_cell(rej):
    rej = pd.Series(rej).fillna(False).astype(bool)
    k, n = int(rej.sum()), len(rej)
    lo, hi = wilson(k, n)
    return {"rate": round(k / n, 4) if n else None, "lo": round(lo, 4),
            "hi": round(hi, 4), "n": n}


def s1_report():
    old = pd.read_csv(os.path.join(HERE, "runs", "full_500", "raw.csv"),
                      low_memory=False)
    new = pd.read_csv(os.path.join(CD, "s1_cd_svd.csv"), low_memory=False)
    out = {"h0_size": {}, "h1_power_n200_rel1": {}, "sentinel": {}, "fails": {}}
    out["sentinel"] = {
        "n": len(new),
        "sentinel_false": int((~new.sentinel_ok.astype(bool)).sum()),
    }
    # H0 尺寸(仅经典 H0 场景;噪声层 H1 的拒绝是假阳性,另行给出)
    h0_scen = ["H0_BASELINE", "H0_MULTIVAR", "H0_CORR"]
    noise_scen = ["H1_HEAVYTAIL", "H1_OUTLIER", "H1_OUTLIER_GRP",
                  "H1_VAR_GROWING", "H1_AR1"]
    gauss = (old.rel > 0) if "rel" in old else None
    for tag, src, rejcol_old, rejcol_new in (
            ("C_old", old, "rej_fitted_hc3", None),
            ("C_svd", new, None, "rej_C_svd"),
            ("D_old", old, "rej_exog_hc3", None),
            ("D_svd", new, None, "rej_D_svd")):
        d = src[~src.fail_fitted_hc3.notna() if rejcol_old ==
                "rej_fitted_hc3" else slice(None)] if False else src
        for n in (50, 200, 1000):
            if tag.startswith("C"):
                col_o, col_n = "rej_fitted_hc3", "rej_C_svd"
            else:
                col_o, col_n = "rej_exog_hc3", "rej_D_svd"
            m_old = old[(old.scenario.isin(h0_scen)) & (old.n == n)
                        & (old.rel > 0)]
            m_new = new[(new.scenario.isin(h0_scen)) & (new.n == n)
                        & (new.rel > 0)]
            key = f"{'C' if tag.startswith('C') else 'D'}_n{n}"
            if key not in out["h0_size"]:
                out["h0_size"][key] = {
                    "old": rate_cell(m_old[col_o].astype(bool)),
                    "svd": rate_cell(m_new[col_n].astype(bool))}
    # 结构 H1 功效 n=200 rel=1%
    scen_all = [s for s in old.scenario.unique()
                if s not in h0_scen + noise_scen]
    for scen in scen_all:
        mo = old[(old.scenario == scen) & (old.n == 200) & (old.rel == 0.01)]
        mn = new[(new.scenario == scen) & (new.n == 200) & (new.rel == 0.01)]
        out["h1_power_n200_rel1"][scen] = {
            "C_old": rate_cell(mo.rej_fitted_hc3.astype(bool)),
            "C_svd": rate_cell(mn.rej_C_svd.astype(bool)),
            "D_old": rate_cell(mo.rej_exog_hc3.astype(bool)),
            "D_svd": rate_cell(mn.rej_D_svd.astype(bool))}
    # 噪声层(应为假阳性,须保持安静)
    for scen in noise_scen:
        mo = old[(old.scenario == scen) & (old.n == 200) & (old.rel == 0.01)]
        mn = new[(new.scenario == scen) & (new.n == 200) & (new.rel == 0.01)]
        out["fails"][scen] = {
            "C_old": rate_cell(mo.rej_fitted_hc3.astype(bool))["rate"],
            "C_svd": rate_cell(mn.rej_C_svd.astype(bool))["rate"],
            "D_old": rate_cell(mo.rej_exog_hc3.astype(bool))["rate"],
            "D_svd": rate_cell(mn.rej_D_svd.astype(bool))["rate"]}
    return out


def s2_report():
    shards = [pd.read_csv(os.path.join(HERE, "runs", f"srsd_full_r{i}",
                                       "raw_srsd.csv"), low_memory=False)
              for i in range(5)]
    old = pd.concat(shards, ignore_index=True)
    new = pd.read_csv(os.path.join(CD, "s2_cd_svd.csv"), low_memory=False)
    H0 = {"i.12.1", "i.29.4", "i.34.27", "i.38.12", "i.39.10", "i.39.22",
          "i.43.31", "i.47.23", "ii.27.16", "ii.27.18", "ii.8.31",
          "iii.12.43", "iii.15.14"}
    new["klass"] = np.where(new.task.str[-6:].isin(
        [t[-6:] for t in H0]) | new.task.isin([f"feynman-{t}" for t in H0]),
        "H0", "H1")
    out = {"h0_size": {}, "h1_power_n200": {}, "sentinel": {}}
    out["sentinel"] = {"n": len(new),
                       "sentinel_false": int((~new.sentinel_ok.astype(bool)).sum())}
    for n in (50, 200, 1000):
        mo = old[(old.klass == "H0") & (old.n == n) & (old.rel > 0)]
        mn = new[(new["klass"] == "H0") & (new["n"] == n)
                 & (new["rel"] > 0)]
        out["h0_size"][f"C_n{n}"] = {
            "old": rate_cell(mo.rej_fitted_hc3.astype(bool)),
            "svd": rate_cell(mn.rej_C_svd.astype(bool))}
        out["h0_size"][f"D_n{n}"] = {
            "old": rate_cell(mo.rej_exog_hc3.astype(bool)),
            "svd": rate_cell(mn.rej_D_svd.astype(bool))}
    mo = old[(old.klass == "H1") & (old.n == 200) & (old.rel >= 0.005)]
    mn = new[(new["klass"] == "H1") & (new["n"] == 200)
             & (new["rel"] >= 0.005)]
    out["h1_power_n200"]["C"] = {"old": rate_cell(mo.rej_fitted_hc3.astype(bool)),
                                 "svd": rate_cell(mn.rej_C_svd.astype(bool))}
    out["h1_power_n200"]["D"] = {"old": rate_cell(mo.rej_exog_hc3.astype(bool)),
                                 "svd": rate_cell(mn.rej_D_svd.astype(bool))}
    return out


def s3_report():
    new = pd.read_csv(os.path.join(CD, "s3_cd_svd.csv"), low_memory=False)
    old = pd.read_csv(r"D:\dimreflex\experiments\form_check_vs_reset"
                      r"\runs\feynman_families\raw_families.csv",
                      low_memory=False)
    out = {"per_family_size_h0": {}, "sentinel": {}}
    out["sentinel"] = {"n": len(new),
                       "sentinel_false": int((~new.sentinel_ok.astype(bool)).sum())}
    key_cols = ["eq_name", "rel", "hetero", "range_scale", "seed"]
    old2 = old.set_index(key_cols)
    idx = pd.MultiIndex.from_frame(new[key_cols])
    new["rej_C_old"] = old2.loc[idx, "rej_fitted_hc3"].values
    new["rej_D_old"] = old2.loc[idx, "rej_exog_hc3"].values
    h0 = new[new.eq_name.isin(
        old[old.klass == "H0"].eq_name.unique())]
    out["per_family_size_h0"]["C_old"] = rate_cell(h0.rej_C_old.astype(bool))
    out["per_family_size_h0"]["C_svd"] = rate_cell(h0.rej_C_svd.astype(bool))
    out["per_family_size_h0"]["D_old"] = rate_cell(h0.rej_D_old.astype(bool))
    out["per_family_size_h0"]["D_svd"] = rate_cell(h0.rej_D_svd.astype(bool))
    new.drop(columns=["rej_C_old", "rej_D_old"]).to_csv(
        os.path.join(CD, "s3_cd_svd.csv"), index=False)
    return out


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "s1"
    summary = {}
    if which == "s1":
        summary["S1"] = s1_report()
    elif which == "s2":
        summary["S2"] = s2_report()
    elif which == "s3":
        summary["S3"] = s3_report()
    path = os.path.join(CD, f"summary_{which}.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=1, ensure_ascii=False)
    print(json.dumps(summary, indent=1, ensure_ascii=False))
    print("->", path)


if __name__ == "__main__":
    main()
