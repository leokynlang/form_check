"""joint_silence.py — "beautiful result, wrong function" joint-silence audit.

Question (user's reframing, 2026-09-09):
NOT "is form_check better?" but: when the fitted functional form is WRONG,
how often does EVERY falsification check stay quiet — i.e. the pipeline
certifies a misfit model as if it were right?

Event of interest, per dataset: no gate rejects at its threshold.
Reported per bench, per n, at two thresholds:
  alpha    = nominal 0.05 (the deployed default)
  calibrated = H0-pool 5th-percentile threshold (matched 5% budget)

Outputs runs/joint_silence/:
  joint_silence_synthetic.csv  per (scenario, n, rel>0) + pool rows
  joint_silence_srsd.csv       per (task, n), SVD-path C/D
  conditional_extrapolation.csv  SRSD: extrap error among silent vs flagged
  summary.json

Disclosures baked in:
- rel = 0 excluded: zero-noise cells emit no p-values (fail-registered,
  degenerate fits) — no certificate exists there, silence is undefined.
- feynman-i.30.5 withdrawn (v5 App B); its replay-side runner_error rows
  dropped via sentinel_ok.
- Noise-layer H1s are wrong-error-model, not wrong-form: excluded from the
  structural pool (semantic separation per §3.3).
- Observational-data cells (World Bank, engineering bench) are NOT pooled
  here: no ground truth for "the form is wrong".
"""
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "runs", "joint_silence")
os.makedirs(OUT, exist_ok=True)

PNAMES = ["p_fc", "p_fitted_f", "p_fitted_hc3", "p_exog_hc3"]
GATES = ["fc", "fitted_f", "fitted_hc3", "exog_hc3"]
ALPHA = 0.05

def reject(ps, tau):
    """p < threshold with NaN (failed evaluation) -> False (counted silent?
    No: NaN p means the gate DECLINED, registered as fail_* — the framework
    does not certify. We count declined as NOT silent (an honest decline is
    not a beautiful certificate)."""
    r = (ps < tau).fillna(False).values
    declined = ps.isna().values
    return r | declined, declined

# ================================================================= S1 synthetic
s1 = pd.read_parquet(os.path.join(HERE, "runs", "full", "raw.parquet"))
s1_struct = s1[(s1.klass == "H1") & (~s1.family.str.startswith("noise"))
               & (s1.rel > 0)].copy()
s1_h0 = s1[(s1.klass == "H0") & (s1.rel > 0)].copy()

rows = []
for (n, rel), h1g in s1_struct.groupby(["n", "rel"]):
    h0g = s1_h0[(s1_h0.n == n) & (s1_h0.rel == rel)]
    tau = {p: float(np.nanquantile(h0g[p], 0.05)) for p in PNAMES}
    rej = {g: reject(h1g[p], tau[p])[0] for p, g in zip(PNAMES, GATES)}
    any_rej = np.zeros(len(h1g), dtype=bool)
    any_decl = np.zeros(len(h1g), dtype=bool)
    for p, g in zip(PNAMES, GATES):
        r, d = reject(h1g[p], tau[p])
        any_rej |= r
        any_decl |= d
    row = {"n": n, "rel": rel, "N": len(h1g),
           "N_declined": int(any_decl.sum()),
           **{f"tpr_{g}": float(rej[g].mean()) for g in GATES},
           "tpr_any": float(any_rej.mean()),
           "joint_silence": float((~any_rej).mean())}
    rows.append({"scenario": "STRUCT_pool", **row})
    for scen, g2 in h1g.groupby("scenario"):
        r_any = np.zeros(len(g2), dtype=bool)
        for p, g in zip(PNAMES, GATES):
            r_any |= reject(g2[p], tau[p])[0]
        rows.append({"scenario": scen, "n": n, "rel": rel, "N": len(g2),
                     "N_declined": int(g2[PNAMES].isna().any(axis=1).sum()),
                     **{f"tpr_{g}": float(reject(g2[p], tau[p])[0].mean())
                        for p, g in zip(PNAMES, GATES)},
                     "tpr_any": float(r_any.mean()),
                     "joint_silence": float((~r_any).mean())})
syn = pd.DataFrame(rows)

# ============================================================ S2 SRSD (SVD rerun)
s2 = pd.read_csv(os.path.join(HERE, "runs", "cd_svd", "s2_cd_svd.csv"),
                 low_memory=False)
s2 = s2[s2.sentinel_ok == True]  # noqa: E712
s2 = s2[s2.task != "feynman-i.30.5"]
base = pd.read_parquet(os.path.join(HERE, "runs", "srsd_full", "raw_srsd.parquet"))
base = base[["task", "n", "rel", "seed", "p_fc", "p_fitted_f", "klass",
             "data_fail", "extrap_med_relerr", "extrap_far_relerr"]].copy()
m = base.merge(s2[["task", "n", "rel", "seed", "p_C_svd", "p_D_svd"]],
               on=["task", "n", "rel", "seed"], how="inner")

PAIRS = {"fc": "p_fc", "fitted_f": "p_fitted_f",
         "fitted_hc3": "p_C_svd", "exog_hc3": "p_D_svd"}

def silent_flag(g, tau):
    """True = every gate quiet AND no gate declined (a clean certificate).
    tau keyed by gate name."""
    not_rej = np.ones(len(g), dtype=bool)
    no_decl = np.ones(len(g), dtype=bool)
    for gname, p in PAIRS.items():
        r, d = reject(g[p], tau[gname])
        not_rej &= ~r
        no_decl &= ~d
    return not_rej & no_decl

def tau_for(h0g, key):
    return float(np.nanquantile(h0g[key], 0.05))

h0 = m[m.klass == "H0"]
h1 = m[m.klass == "H1"]

# ---- per (task, n) ----
rows2 = []
for (task, n), g in h1.groupby(["task", "n"]):
    h0g = h0[h0.n == n]
    tau_a = {k: ALPHA for k in PAIRS}
    tau_c = {k: tau_for(h0g, v) for k, v in PAIRS.items()}
    sil_a = silent_flag(g, tau_a)
    sil_c = silent_flag(g, tau_c)
    far = g["extrap_far_relerr"].values
    row = {"task": task, "n": n, "N": len(g),
           **{f"tpr_{k}": float((g[v] < ALPHA).fillna(False).mean())
              for k, v in PAIRS.items()},
           "tpr_any_alpha": float((~sil_a).mean()),
           "joint_silence_alpha": float(sil_a.mean()),
           "tpr_any_cal": float((~sil_c).mean()),
           "joint_silence_cal": float(sil_c.mean()),
           "med_farerr_silent_alpha": (float(np.median(far[sil_a]))
                                       if sil_a.sum() >= 5 else np.nan),
           "med_farerr_flagged_alpha": (float(np.median(far[~sil_a]))
                                        if (~sil_a).sum() >= 5 else np.nan)}
    rows2.append(row)
srsd_task = pd.DataFrame(rows2)

# ---- pool per n ----
rows3 = []
for n, g in h1.groupby("n"):
    h0g = h0[h0.n == n]
    tau_a = {k: ALPHA for k in PAIRS}
    tau_c = {k: tau_for(h0g, v) for k, v in PAIRS.items()}
    sil_a = silent_flag(g, tau_a)
    sil_c = silent_flag(g, tau_c)
    rows3.append({"n": n, "N": len(g),
                  **{f"tpr_{k}": float((g[v] < ALPHA).fillna(False).mean())
                     for k, v in PAIRS.items()},
                  "tpr_any_alpha": float((~sil_a).mean()),
                  "joint_silence_alpha": float(sil_a.mean()),
                  "tpr_any_cal": float((~sil_c).mean()),
                  "joint_silence_cal": float(sil_c.mean())})
srsd_pool = pd.DataFrame(rows3)

# ---- conditional extrapolation (pool, per n): silent vs flagged ----
cond = []
for n, g in h1.groupby("n"):
    h0g = h0[h0.n == n]
    tau_a = {k: ALPHA for k in PAIRS}
    tau_c = {k: tau_for(h0g, v) for k, v in PAIRS.items()}
    for lbl, sil in [("alpha", silent_flag(g, tau_a)),
                     ("cal", silent_flag(g, tau_c))]:
        far = g["extrap_far_relerr"].values
        med_all = g["extrap_med_relerr"].values
        ok = ~np.isnan(far)                  # data_fail rows: no extrapolation
        sil_ok = sil & ok
        flg_ok = (~sil) & ok
        cond.append({
            "n": n, "threshold": lbl, "N_silent": int(sil.sum()),
            "N_silent_ok": int(sil_ok.sum()),
            "med_farerr_silent": (float(np.median(far[sil_ok]))
                                  if sil_ok.sum() >= 5 else np.nan),
            "p90_farerr_silent": (float(np.quantile(far[sil_ok], 0.9))
                                  if sil_ok.sum() >= 5 else np.nan),
            "frac_silent_far_gt10pct": (float((far[sil_ok] > 0.10).mean())
                                        if sil_ok.sum() >= 5 else np.nan),
            "med_farerr_flagged": (float(np.median(far[flg_ok]))
                                   if flg_ok.sum() >= 5 else np.nan),
            "med_mederr_silent": (float(np.median(med_all[sil_ok]))
                                  if sil_ok.sum() >= 5 else np.nan),
            "med_mederr_flagged": (float(np.median(med_all[flg_ok]))
                                   if flg_ok.sum() >= 5 else np.nan)})
cond = pd.DataFrame(cond)

# ==================================================================== write
syn.to_csv(os.path.join(OUT, "joint_silence_synthetic.csv"), index=False)
srsd_task.to_csv(os.path.join(OUT, "joint_silence_srsd_per_task.csv"), index=False)
srsd_pool.to_csv(os.path.join(OUT, "joint_silence_srsd_pool.csv"), index=False)
cond.to_csv(os.path.join(OUT, "conditional_extrapolation.csv"), index=False)

summary = {
    "question": "P(all gates quiet AND none declined | wrong functional form)",
    "thresholds": {"alpha": "nominal .05",
                   "calibrated": "H0-pool 5th-pctile p threshold per gate"},
    "disclosures": [
        "rel=0 excluded (no certificates emitted; fail-registered)",
        "feynman-i.30.5 withdrawn (v5 App B); sentinel_ok filter applied",
        "noise-layer H1 excluded (wrong error model, not wrong form)",
        "declined evaluations (NaN p, fail_*) counted as NOT silent"],
    "synthetic_pool": syn[syn.scenario == "STRUCT_pool"].to_dict("records"),
    "srsd_pool": srsd_pool.to_dict("records"),
    "conditional_extrapolation": cond.to_dict("records"),
}
with open(os.path.join(OUT, "summary.json"), "w", encoding="utf-8") as f:
    json.dump(summary, f, ensure_ascii=False, indent=1, default=float)

pd.set_option("display.width", 250)
pool = syn[syn.scenario == "STRUCT_pool"]
print("=== SYNTHETIC pool (calibrated 5%) ===")
print(pool.to_string(index=False))
print()
print("=== SYNTHETIC: scenarios with nonzero joint silence (any n, rel) ===")
nz = syn[(syn.scenario != "STRUCT_pool") & (syn.joint_silence > 0)]
print(nz.to_string(index=False))
print()
print("=== SRSD pool ===")
print(srsd_pool.to_string(index=False))
print()
print("=== SRSD per task: worst silence (n=50, alpha) ===")
t50 = srsd_task[srsd_task.n == 50].sort_values("joint_silence_alpha",
                                              ascending=False)
print(t50.head(10).to_string(index=False))
print()
print("=== SRSD: silent cells' extrapolation error vs flagged ===")
print(cond.to_string(index=False))
print()
bad = srsd_task[(srsd_task.n == 50)
                & (srsd_task.med_farerr_silent_alpha > 0.02)]
print("=== DANGEROUS SILENCE: silent-certified tasks w/ median far-extrap err > 2% (n=50) ===")
print(bad[["task", "N", "joint_silence_alpha", "joint_silence_cal",
           "med_farerr_silent_alpha"]].to_string(index=False))
