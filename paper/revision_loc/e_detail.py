# -*- coding: utf-8 -*-
"""e_detail.py — rating-curve D 未拒站 + meteo A 保留格明细。"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import pandas as pd

e = pd.read_csv(r"D:\form_check\runs\cd_svd\e_smoke_results_svd.csv")
e = e[e.status == "ok"]
rc = e[e.dataset.str.contains("_rating")]
miss = rc[~rc.rej_exog_hc3.astype(bool)]
cols = ["dataset", "target", "n_used", "p_fc", "p_fitted_f",
        "p_fitted_hc3", "p_exog_hc3"]
print("D 未拒站:")
print(miss[cols].to_string(index=False))
ws = e[e.dataset.str.contains("meteo")]
print("\nmeteo 目标分布:", ws["target"].value_counts().to_dict())
keep = ws[~ws.rej_fc.astype(bool)]
print("A 保留格:")
print(keep[cols].to_string(index=False))
