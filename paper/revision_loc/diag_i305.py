# -*- coding: utf-8 -*-
"""diag_i305.py — i.30.5 哨兵失配诊断:存档 NaN?重算 NaN?还是真失配?"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, r"D:\form_check")
import numpy as np
import pandas as pd

patch = pd.read_csv(r"D:\form_check\runs\cd_svd\s2_i305_patch.csv",
                    low_memory=False)
old = pd.concat([pd.read_csv(rf"D:\form_check\runs\srsd_full_r{i}\raw_srsd.csv",
                             low_memory=False,
                             usecols=["task", "n", "rel", "seed", "p_fitted_f",
                                      "stat_fitted_f"])
                 for i in range(5)], ignore_index=True)
old = old[old.task == "feynman-i.30.5"].reset_index(drop=True)
n_grid, rel_grid = [50, 200, 1000], [0.0, 0.005, 0.01, 0.03, 0.05]
old["_bs"] = [200_000 + int(s) * 97 + n_grid.index(int(n)) * 13
              + rel_grid.index(float(r)) * 5_000_000
              for s, n, r in zip(old.seed, old.n, old.rel)]
m = patch.reset_index(drop=True).join(old[["_bs", "p_fitted_f",
                                           "stat_fitted_f"]])
ok = m.sentinel_ok.astype(bool).values
pf = m.p_fitted_f.values
print("失败行中存档 p NaN 比例:",
      round(float(np.mean(np.isnan(pf[~ok]))), 3))
print("失败行中重算 stat NaN 比例:",
      round(float(np.mean(~np.isfinite(m.stat_F_recheck.values[~ok]))), 3))
both = ~ok & np.isfinite(pf) & np.isfinite(m.stat_F_recheck.values)
print("真失配行数:", int(both.sum()))
import collections
fails = m[~ok]
print("失败行按 n:", collections.Counter(fails.loc[:, "n"].tolist()).most_common()
      if "n" in fails else "n 列不在 patch(见 rerun 脚本列集)")
print("失败行数:", len(fails), "| 列:", list(fails.columns)[:12])
