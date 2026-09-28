import sys
sys.path.insert(0, "D:/form_check")
import numpy as np, pandas as pd
import common

df = pd.read_csv("D:/form_check/real_data/clean/uci_concrete_strength.csv")
df = df[(df.age_day > 0) & (df.strength_mpa > 0)]
train_pool = df[df.age_day <= 28]
hold = df[df.age_day >= 90]
print("train pool:", len(train_pool), "holdout:", len(hold),
      "holdout ages:", sorted(hold.age_day.unique()))

res = []
for seed in range(300):
    tr = train_pool.sample(n=100, replace=True, random_state=seed)
    sub = tr[["strength_mpa", "age_day"]].rename(columns={"strength_mpa": "y", "age_day": "x"})
    row = common.run_all_methods(sub, "y", ["x"], seed)
    lx, ly = np.log(sub[["x"]].values), np.log(sub["y"].values)
    X = np.column_stack([np.ones(len(lx)), lx])
    b = np.linalg.lstsq(X, ly, rcond=None)[0]
    lxh = np.log(hold[["age_day"]].values)
    pred = np.exp(np.column_stack([np.ones(len(lxh)), lxh]) @ b)
    rel = (pred - hold.strength_mpa.values) / hold.strength_mpa.values  # signed: power law over-predicts
    row.update(seed=seed, med_far=float(np.median(np.abs(rel))),
               med_bias=float(np.median(rel)))
    res.append(row)

r = pd.DataFrame(res)
out = "D:/form_check/real_data/concrete_case_study_rows.csv"
r.to_csv(out, index=False)
P = ["p_fc", "p_fitted_f", "p_fitted_hc3", "p_exog_hc3"]
ok = r[P].notna().all(axis=1)
rej = (r[P] < 0.05).any(axis=1)
silent = ok & ~rej
print("n fits:", len(r), "evaluable:", int(ok.sum()), "declines:", int((~ok).sum()))
print("JOINT SILENCE:", int(silent.sum()), "=", round(silent.sum()/ok.sum(), 3))
if silent.any():
    s = r.loc[silent]
    print("silent med_far: median", round(float(s.med_far.median()),3),
          "| P(med_far>10%):", round(float((s.med_far>0.10).mean()),3),
          "| median signed bias:", round(float(s.med_bias.median()),3))
f = r.loc[~silent]
print("flagged med_far median:", round(float(f.med_far.median()),3),
      "| P(med_far>10%):", round(float((f.med_far>0.10).mean()),3))
print("written:", out)
