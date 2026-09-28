# -*- coding: utf-8 -*-
"""e3_counts_svd.py — E 台架 designed-regime 计数(SVD 重算口径)。"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import pandas as pd

e = pd.read_csv(r"D:\form_check\runs\cd_svd\e_smoke_results_svd.csv")
e = e[e.status == "ok"]
H0_SIDE = ["ccpp", "gas_turbine", "rdataset_cars", "rdataset_trees",
           "auto_mpg", "meteo", "steel", "usgs_", "epa"]
H1_SIDE = ["webbook", "yacht", "concrete", "puromycin", "tetouan",
           "usgs_eq", "housing", "ailerons", "cpu_act"]


def side(name):
    n = str(name).lower()
    if any(k in n for k in H1_SIDE):
        return "H1"
    if any(k in n for k in H0_SIDE):
        return "H0"
    return None


e["side"] = e.dataset.map(side)
sub = e[e.side.notna()]
g = sub.groupby("side")[["rej_fc", "rej_fitted_f", "rej_fitted_hc3",
                         "rej_exog_hc3"]].agg(["sum", "count"])
print("designed-regime 计数(README 考场口径):")
print(g.to_string())
rc = e[e.dataset.str.contains("usgs_1") & ~e.dataset.str.contains("eq")]
print(f"\nrating-curve 格数: {len(rc)} | A 拒 {int(rc.rej_fc.sum())} | "
      f"B 拒 {int(rc.rej_fitted_f.sum())} | C(svd) 拒 {int(rc.rej_fitted_hc3.sum())} | "
      f"D(svd) 拒 {int(rc.rej_exog_hc3.sum())}")
ws = e[e.dataset.str.contains("meteo")]
print(f"meteo 格数: {len(ws)} | A 拒 {int(ws.rej_fc.sum())} | "
      f"B 拒 {int(ws.rej_fitted_f.sum())}")
pu = e[e.dataset.str.contains("puromycin")]
print("puromycin:", pu[["dataset", "target", "n_used", "p_fc", "p_fitted_f",
                        "p_fitted_hc3", "p_exog_hc3"]].to_string(index=False))
EOF_MARK = None
