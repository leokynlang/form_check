import sys
sys.path.insert(0, "D:/form_check")
import numpy as np, pandas as pd
import common

df = pd.read_csv("D:/form_check/real_data/clean/rdataset_pressure.csv")
df = df[(df.temperature > 0) & (df.pressure > 0)].reset_index(drop=True)
train_full = df[df.temperature <= 200]
hold = df[df.temperature > 200]
print("train pts:", len(train_full), " holdout pts:", len(hold))
print("holdout T:", sorted(hold.temperature.values))

res = []
for seed in range(200):
    tr = train_full.sample(n=12, replace=True, random_state=seed)
    sub = tr.rename(columns={"pressure": "y", "temperature": "x"})
    row = common.run_all_methods(sub, "y", ["x"], seed)
    lx, ly = np.log(sub[["x"]].values), np.log(sub["y"].values)
    X = np.column_stack([np.ones(len(lx)), lx])
    b = np.linalg.lstsq(X, ly, rcond=None)[0]
    lxh = np.log(hold[["temperature"]].values)
    pred = np.exp(np.column_stack([np.ones(len(lxh)), lxh]) @ b)
    rel = np.abs(pred - hold.pressure.values) / hold.pressure.values
    row.update(seed=seed, med_far=float(np.median(rel)), max_far=float(rel.max()))
    res.append(row)

r = pd.DataFrame(res)
out = "D:/form_check/real_data/pressure_case_study_rows.csv"
r.to_csv(out, index=False)
P = ["p_fc", "p_fitted_f", "p_fitted_hc3", "p_exog_hc3"]
ok = r[P].notna().all(axis=1)
rej = (r[P] < 0.05).any(axis=1)
silent = ok & ~rej
print("n fits:", len(r), "evaluable:", int(ok.sum()), "declines:", int((~ok).sum()))
print("JOINT SILENCE:", int(silent.sum()), "=", round(silent.sum()/ok.sum(), 3))
if silent.any():
    print("far-field err (silent): median", round(float(r.loc[silent,"med_far"].median()),2),
          " p90", round(float(r.loc[silent,"med_far"].quantile(.9)),2),
          " P(med_far>50%):", round(float((r.loc[silent,"med_far"]>0.5).mean()),2))
print("far-field err (flagged): median", round(float(r.loc[~silent,"med_far"].median()),2))
print("enrichment:", round(float((r.loc[silent,"med_far"]>0.5).mean()/max((r.loc[~silent,"med_far"]>0.5).mean(),1e-9)),2))
print("written:", out)
