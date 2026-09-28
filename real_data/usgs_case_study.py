"""Real-data case study: USGS 11447650 rating curve.

Question: does the four-gate panel certify power-law fits that are silent
in-domain but demonstrably wrong at the far field, on REAL data (no
simulation)? Truth here is not closed-form: it is held-out real USGS
measurements above the 90th percentile of training gage height. The power
law Q = a*G^b is the deployed rating-curve approximation; the true physical
relation Q = a*(G-g0)^b contains a datum offset g0, so fitting without the
offset plants a misspecification whose far-field signature is real.
"""
import sys, os
sys.path.insert(0, "D:/form_check")
import numpy as np, pandas as pd
import common

df = pd.read_csv("D:/form_check/real_data/clean/usgs_11447650_rating.csv")
df = df[(df.gage_height_ft > 0) & (df.discharge_cfs > 0)].dropna()
print("rows:", len(df), "G range:", df.gage_height_ft.min(), df.gage_height_ft.max())

df = df.sort_values("gage_height_ft").reset_index(drop=True)
cut = df.gage_height_ft.quantile(0.90)
train = df[df.gage_height_ft <= cut]
hold = df[df.gage_height_ft > cut]
print("train:", len(train), "holdout:", len(hold), "cut G:", round(cut, 2))

res = []
for seed in range(50):
    tr = train.sample(n=min(200, len(train)), replace=True, random_state=seed)
    sub = tr[["discharge_cfs", "gage_height_ft"]].rename(
        columns={"discharge_cfs": "y", "gage_height_ft": "x"})
    row = common.run_all_methods(sub, "y", ["x"], seed)
    lx, ly = np.log(sub[["x"]].values), np.log(sub["y"].values)
    X = np.column_stack([np.ones(len(lx)), lx])
    b = np.linalg.lstsq(X, ly, rcond=None)[0]
    lxh = np.log(hold[["gage_height_ft"]].values)
    pred = np.exp(np.column_stack([np.ones(len(lxh)), lxh]) @ b)
    rel = np.abs(pred - hold.discharge_cfs.values) / hold.discharge_cfs.values
    row.update(seed=seed, med_far=float(np.median(rel)),
               p90_far=float(np.quantile(rel, 0.9)))
    res.append(row)

r = pd.DataFrame(res)
out = "D:/form_check/real_data/usgs_case_study_rows.csv"
r.to_csv(out, index=False)
P = ["p_fc", "p_fitted_f", "p_fitted_hc3", "p_exog_hc3"]
ok = r[P].notna().all(axis=1)
rej = (r[P] < 0.05).any(axis=1)
silent = ok & ~rej
print("\nn fits:", len(r), " evaluable:", int(ok.sum()))
print("joint silence:", int(silent.sum()), "=", round(silent.sum() / ok.sum(), 3))
print("median far err | silent:",
      round(float(r.loc[silent, "med_far"].median()), 4) if silent.any() else "n/a")
print("median far err | flagged:",
      round(float(r.loc[~silent, "med_far"].median()), 4))
print("P(med_far>10% | silent):",
      round(float((r.loc[silent, "med_far"] > 0.10).mean()), 3) if silent.any() else "n/a")
print("P(med_far>10% | flagged):",
      round(float((r.loc[~silent, "med_far"] > 0.10).mean()), 3))
print("written:", out)
