# -*- coding: utf-8 -*-
"""task_bootstrap_v5.py — Table 2(v5)的任务级区间:C/D 用 SVD 新值。

口径与 v4 Appendix C 一致:自助重抽 12 个 H0 方程(10,000 次),
报告池化尺寸的 bootstrap 均值区间;另给 LOEO 范围。i.30.5 不在 H0 池,
不受剔除影响。
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, r"D:\form_check")
import numpy as np
import pandas as pd

HERE = r"D:\form_check"
H0 = ["feynman-" + t for t in
      ["i.12.1", "i.29.4", "i.34.27", "i.38.12", "i.39.10", "i.39.22",
       "i.43.31", "i.47.23", "ii.27.16", "ii.27.18", "ii.8.31",
       "iii.12.43", "iii.15.14"]]


def main():
    new = pd.read_csv(rf"{HERE}\runs\cd_svd\s2_cd_svd.csv", low_memory=False)
    new = new[new.task.isin(H0) & (new["rel"] > 0)]
    per_task = new.groupby("task").agg(
        n=("rel", "size"),
        c=("rej_C_svd", "mean"), d=("rej_D_svd", "mean"),
        c_old=("rej_C_svd", "size"))
    # 旧值同池(用于并排)
    old = pd.concat(
        [pd.read_csv(rf"{HERE}\runs\srsd_full_r{i}\raw_srsd.csv",
                     low_memory=False) for i in range(5)], ignore_index=True)
    old = old[old.task.isin(H0) & (old.rel > 0)]
    per_old = old.groupby("task").agg(c=("rej_fitted_hc3", "mean"),
                                      d=("rej_exog_hc3", "mean"))
    print("=== 逐任务 H0 拒绝率(n=200 合并噪声;旧 Gram -> 新 SVD) ===")
    for t in H0:
        print(f"{t:22s} C {per_old.loc[t,'c']:.3f} -> {per_task.loc[t,'c']:.3f}   "
              f"D {per_old.loc[t,'d']:.3f} -> {per_task.loc[t,'d']:.3f}")
    rng = np.random.default_rng(20260907)
    tasks = np.array(H0)
    reps = 10000
    for n in (200,):
        sub = new[new["n"] == n]
        rates = sub.groupby("task").apply(
            lambda g: pd.Series({"c": g.rej_C_svd.mean(),
                                 "d": g.rej_D_svd.mean()}), include_groups=False)
        subo = old[old.n == n]
        rateso = subo.groupby("task").apply(
            lambda g: pd.Series({"c": g.rej_fitted_hc3.mean(),
                                 "d": g.rej_exog_hc3.mean()}), include_groups=False)
        for meth in ("c", "d"):
            vals = rates[meth].reindex(tasks).values
            valso = rateso[meth].reindex(tasks).values
            bs_new, bs_old = [], []
            for _ in range(reps):
                pick = rng.integers(0, len(tasks), len(tasks))
                bs_new.append(vals[pick].mean())
                bs_old.append(valso[pick].mean())
            lo_n, hi_n = np.percentile(bs_new, [2.5, 97.5])
            lo_o, hi_o = np.percentile(bs_old, [2.5, 97.5])
            print(f"n={n} {'C' if meth=='c' else 'D'}: "
                  f"task-bootstrap 旧 {np.mean(valso):.3f} [{lo_o:.3f},{hi_o:.3f}]"
                  f" -> 新 {np.mean(vals):.3f} [{lo_n:.3f},{hi_n:.3f}]")
        # LOEO
        for meth in ("c", "d"):
            vals = rates[meth].reindex(tasks).values
            valso = rateso[meth].reindex(tasks).values
            weights = sub.groupby("task").size().reindex(tasks).values
            loo_n = [np.average(np.delete(vals, i), weights=np.delete(weights, i))
                     for i in range(len(tasks))]
            loo_o = [np.average(np.delete(valso, i), weights=np.delete(weights, i))
                     for i in range(len(tasks))]
            print(f"n={n} {'C' if meth=='c' else 'D'} LOEO: "
                  f"旧 [{min(loo_o):.3f},{max(loo_o):.3f}] -> "
                  f"新 [{min(loo_n):.3f},{max(loo_n):.3f}]")


if __name__ == "__main__":
    main()
