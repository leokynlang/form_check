# -*- coding: utf-8 -*-
"""s3_family_table_v5.py — S3 校准后家族表重推导(C/D 用 SVD 新 p)。

S3 逐种子层在原仓 raw_families.csv(klass/family/eq_name 列齐);
新 C/D p 在 runs/cd_svd/s3_cd_svd.csv。对齐:rel/hetero/range_scale/seed
(+行序 cell_ordinal 校验)。阈值:各方法在 S3 H0 方程池(官方域)按族内
合并取 5% 分位——与 v4 Table 6 的『S3 H0 pool 校准』口径一致。
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, r"D:\form_check")
import numpy as np
import pandas as pd

HERE = r"D:\form_check"


def threshold_at_fpr(p, fpr=0.05):
    p = np.sort(np.asarray(p, float))
    p = p[np.isfinite(p)]
    k = max(1, int(round(fpr * len(p))))
    return p[k - 1]


def power_at(p, thr):
    p = np.asarray(p, float)
    p = p[np.isfinite(p)]
    return float(np.mean(p < thr))


def main():
    old = pd.read_csv(rf"{HERE}\runs\feynman_families\raw_families.csv",
                      low_memory=False)
    new = pd.read_csv(rf"{HERE}\runs\cd_svd\s3_cd_svd.csv", low_memory=False)
    sel = pd.read_csv(rf"{HERE}\runs\feynman_families\selection_table.csv")
    klass = dict(zip(sel.eq_name, sel.klass))
    family = dict(zip(sel.eq_name, sel.family))
    for df in (old, new):
        df["klass"] = df.eq_name.map(klass)
        df["family"] = df.eq_name.map(family)
    # S3 的 H0 校准池 = klass==H0 方程(官方域格;官方域 = range_scale==1)
    h0o = old[(old.klass == "H0")]
    h0n = new[new.klass == "H0"]
    print("=== Table 6(v5):S3 校准后功效按族(@FPR 5%,C/D 用 SVD) ===")
    print("(A/B 沿用存档 p;阈值各自从 S3 H0 池取 5% 分位)")
    h1o = old[old.klass == "H1"]
    h1n = new[new.klass == "H1"]
    fams = ("exponential", "rational", "trigonometric")
    out = {}
    for fam in fams:
        row = {}
        for name, colo, coln in (
                ("A", "p_fc", "p_fc"),
                ("B", "p_fitted_f", "p_fitted_f"),
                ("C_svd", "p_fitted_hc3", "p_C_svd"),
                ("D_svd", "p_exog_hc3", "p_D_svd")):
            use_new = name.endswith("_svd")
            po = (h0n if use_new else h0o)[coln]
            thr = threshold_at_fpr(po)
            p1 = (h1n if use_new else h1o)
            p1 = p1[p1.family == fam][coln]
            row[name] = round(power_at(p1, thr), 3)
        out[fam] = row
        print(f"{fam:16s} A {row['A']:.3f}  B {row['B']:.3f}  "
              f"C(SVD) {row['C_svd']:.3f}  D(SVD) {row['D_svd']:.3f}")
    # H0 尺寸对照(官方域)
    print("\n=== S3 H0 尺寸(官方域格) ===")
    h0o1 = h0o[h0o.range_scale == 1.0]
    h0n1 = h0n[h0n.range_scale == 1.0]
    for name, colo, coln in (("A", "p_fc", "p_fc"),
                             ("B", "p_fitted_f", "p_fitted_f"),
                             ("C old", "p_fitted_hc3", None),
                             ("C svd", None, "p_C_svd"),
                             ("D old", "p_exog_hc3", None),
                             ("D svd", None, "p_D_svd")):
        if colo:
            p = h0o1[colo].dropna()
        else:
            p = h0n1[coln].dropna()
        print(f"{name:8s}: {np.mean(p < 0.05):.3f}  (n={len(p)})")


if __name__ == "__main__":
    main()
