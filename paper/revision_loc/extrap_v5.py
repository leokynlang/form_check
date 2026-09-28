# -*- coding: utf-8 -*-
"""extrap_v5.py — Table 8 外推证据复核(31 任务池,i.30.5 剔除后)。"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import pandas as pd

old = pd.concat([pd.read_csv(rf"D:\form_check\runs\srsd_full_r{i}\raw_srsd.csv",
                             low_memory=False,
                             usecols=["task", "klass", "n", "rel", "seed",
                                      "extrap_med_relerr", "p_fc", "p_fitted_f"])
                 for i in range(5)], ignore_index=True)
old = old[old.task != "feynman-i.30.5"]
sub = old[old.rel >= 0.005]
sub = old[old.extrap_med_relerr.notna()]
base = float((sub.extrap_med_relerr > 0.10).mean())
print(f"外推基率(31 任务池): {base:.3f}  (v4 报 .198 为 32 任务池)")
for c, lbl in (("p_fc", "A"), ("p_fitted_f", "B")):
    alarm = (sub[c] < .05) & sub[c].notna()
    severe = sub.extrap_med_relerr > 0.10
    tp = alarm & severe
    prec = float(tp.sum() / alarm.sum())
    rec = float(tp.sum() / severe.sum())
    print(f"{lbl}: precision {prec:.2f}  recall {rec:.2f}  "
          f"(n_alarms {int(alarm.sum())})")
