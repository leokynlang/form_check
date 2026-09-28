# -*- coding: utf-8 -*-
"""calibrated_power_v5.py — Table 3(v5):S2 校准后功效,C/D 全链路 SVD。

与 v4 口径完全一致(srsd_analyze.py::tau = H0 池 p 的 5% 分位,
功效 = H1 池 p <= tau 的比例),唯一区别是 C/D 的 p 来自 SVD 路径。
v4 的 tau_C/tau_D 是在伪影 p 分布上取的(42-56% 质量堆积在 0,tau=0),
因此其 C/D『校准功效』沿伪影尺度计算,应整体作废重述。
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, r"D:\form_check")
import json
import numpy as np
import pandas as pd

HERE = r"D:\form_check"
H0 = ["feynman-" + t for t in
      ["i.12.1", "i.29.4", "i.34.27", "i.38.12", "i.39.10", "i.39.22",
       "i.43.31", "i.47.23", "ii.27.16", "ii.27.18", "ii.8.31",
       "iii.12.43", "iii.15.14"]]


def main():
    new = pd.read_csv(rf"{HERE}\runs\cd_svd\s2_cd_svd.csv", low_memory=False)
    # A/B 用存档 p(未受伪影影响);C/D 用 SVD p
    old_full = pd.concat(
        [pd.read_csv(rf"{HERE}\runs\srsd_full_r{i}\raw_srsd.csv",
                     low_memory=False,
                     usecols=["task", "n", "rel", "seed", "p_fc", "p_fitted_f"])
         for i in range(5)], ignore_index=True)
    old_full = old_full[old_full.task != "feynman-i.30.5"]
    new = new.merge(old_full, on=["task", "n", "rel", "seed"], how="left")
    tasks_all = [t for t in pd.concat([old_full.task, new.task]).unique()]
    v4 = {"A": {50: .41, 200: .51, 1000: .59},
          "B": {50: .19, 200: .28, 1000: .42},
          "C": {50: .05, 200: .28, 1000: .41},
          "D": {50: .19, 200: .37, 1000: .65}}
    print(f"{'方法':20s} {'n=50':>22s} {'n=200':>22s} {'n=1000':>22s}")
    out = {}
    for name, pcol, is_svd in (("A (permutation)", "p_fc", False),
                               ("B (RESET fitted F)", "p_fitted_f", False),
                               ("C (HC3 fitted, SVD)", "p_C_svd", True),
                               ("D (HC3 exog, SVD)", "p_D_svd", True)):
        row = {}
        for n in (50, 200, 1000):
            h0 = new[(new.task.isin(H0)) & (new["n"] == n) & (new["rel"] > 0)]
            h1 = new[(~new.task.isin(H0)) & (new["n"] == n)
                     & (new["rel"] >= 0.005)]
            tau = float(np.nanquantile(h0[pcol].values, 0.05))
            p = h1[pcol].values
            p = p[np.isfinite(p)]
            pw = float(np.mean(p <= tau))
            row[n] = (round(pw, 3), tau)
        out[name] = row
        cells = []
        for n in (50, 200, 1000):
            pw, tau = row[n]
            old = v4[name.split()[0]]
            cells.append(f"{pw:.3f} (τ={tau:.3g}; v4 {old[n]:.2f})")
        print(f"{name:20s} {cells[0]:>22s} {cells[1]:>22s} {cells[2]:>22s}")
    with open(rf"{HERE}\runs\cd_svd\v5_table3_calibrated.json", "w",
              encoding="utf-8") as f:
        json.dump({k: {str(n): v for n, v in row.items()}
                   for k, row in out.items()}, f, indent=1)


if __name__ == "__main__":
    main()
