"""Omnibus combination check for the four-gate panel (S2 physics bench).

Question (reviewer-shaped): is the blind-spot pattern of the joint-silence
analysis a failure of AGGREGATION (no panel statistic was formed) or of
DIRECTION COVERAGE (the blind directions are orthogonal to every tested
column set, so no monotone combination of the four marginals can see them)?

Method: combine the four per-row gate p-values with the two standard
omnibus rules --
  Fisher   X = -2 * sum(log p), 8 df, nominal, and H0-pool-calibrated
           (5% empirical size per n stratum);
  Tippett  min-p with the exact 4-test correction 1-(1-pmin)^4.
Recompute joint silence under each rule and compare against the paper's
any-gate rule on identical rows.

Reads the same sources as joint_silence.py:
  runs/srsd_full/raw_srsd.parquet   (A/B p-values, klass, extrapolation)
  runs/cd_svd/s2_cd_svd.csv         (SVD-path C/D p-values)
Writes:
  runs/joint_silence/omnibus_combination.csv
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "runs", "joint_silence")
PAIRS = {"A_fc": "p_fc", "B_fitted_f": "p_fitted_f",
         "fitted_hc3": "p_C_svd", "exog_hc3": "p_D_svd"}
COLS = list(PAIRS.values())


def load():
    base = pd.read_parquet(os.path.join(HERE, "runs", "srsd_full", "raw_srsd.parquet"))
    base = base[["task", "n", "rel", "seed", "p_fc", "p_fitted_f", "klass"]].copy()
    s2 = pd.read_csv(os.path.join(HERE, "runs", "cd_svd", "s2_cd_svd.csv"), low_memory=False)
    s2 = s2[s2.sentinel_ok == True]  # noqa: E712
    s2 = s2[s2.task != "feynman-i.30.5"]
    m = base.merge(s2[["task", "n", "rel", "seed", "p_C_svd", "p_D_svd"]],
                   on=["task", "n", "rel", "seed"], how="inner")
    return m[m.rel > 0]


def pmat(sub):
    P = np.column_stack([pd.to_numeric(sub[c], errors="coerce").values for c in COLS])
    ok = ~np.isnan(P).any(axis=1)
    return ok, np.where(np.isnan(P), 1.0, P)


def silence(sub, rule, tau=0.05):
    ok, P = pmat(sub)
    if rule == "any":
        rej = (P < tau).any(axis=1)
    elif rule == "fisher":
        X = -2.0 * np.log(np.clip(P, 1e-300, 1.0)).sum(axis=1)
        rej = chi2.sf(X, 8) < tau
    elif rule == "fisher_cal":
        raise ValueError("use tau from fisher_tau()")
    elif rule == "tippett":
        pmin = P.min(axis=1)
        rej = (1.0 - (1.0 - pmin) ** 4) < tau
    else:
        raise ValueError(rule)
    return int((ok & ~rej).sum()), int(ok.sum())


def fisher_tau(h0sub):
    """H0-pool 95th percentile of the Fisher X statistic (5% empirical size)."""
    _, P = pmat(h0sub)
    X = -2.0 * np.log(np.clip(P, 1e-300, 1.0)).sum(axis=1)
    return float(np.quantile(X, 0.95))


def fisher_silence(sub, xtau):
    ok, P = pmat(sub)
    X = -2.0 * np.log(np.clip(P, 1e-300, 1.0)).sum(axis=1)
    rej = X > xtau
    return int((ok & ~rej).sum()), int(ok.sum())


def main():
    mn = load()
    rows = []

    # ---- pooled level, by n stratum -------------------------------------
    for n, g in mn.groupby("n"):
        xtau = fisher_tau(g[g.klass == "H0"])
        rec = {"unit": "pool", "n": n, "fisher_xtau_cal": round(xtau, 2),
               "chi2_8_nominal_tau": round(chi2.ppf(0.95, 8), 2)}
        for klass in ["H0", "H1"]:
            sub = g[g.klass == klass]
            for rule in ["any", "fisher", "tippett"]:
                s, d = silence(sub, rule)
                rec[f"{klass}_sil_{rule}"] = round(s / d, 4)
            s, d = fisher_silence(sub, xtau)
            rec[f"{klass}_sil_fishercal"] = round(s / d, 4)
        rec["LR_any"] = round(rec["H0_sil_any"] / rec["H1_sil_any"], 2)
        rec["LR_fishercal"] = round(rec["H0_sil_fishercal"] / rec["H1_sil_fishercal"], 2)
        rows.append(rec)

    # ---- per task x n, H1 side ------------------------------------------
    for (task, n), g in mn[mn.klass == "H1"].groupby(["task", "n"]):
        g0 = mn[(mn.klass == "H0") & (mn.n == n)]
        xtau = fisher_tau(g0)
        rec = {"unit": task, "n": n, "fisher_xtau_cal": round(xtau, 2)}
        for rule in ["any", "fisher", "tippett"]:
            s, d = silence(g, rule)
            rec[f"H1_sil_{rule}"] = round(s / d, 4)
        s, d = fisher_silence(g, xtau)
        rec["H1_sil_fishercal"] = round(s / d, 4)
        rows.append(rec)

    out = pd.DataFrame(rows)
    os.makedirs(OUT, exist_ok=True)
    out.to_csv(os.path.join(OUT, "omnibus_combination.csv"), index=False)

    pool = out[out.unit == "pool"]
    print("=== Pooled S2 (noise-present): joint silence under combination rules ===")
    print(pool[["n", "H0_sil_any", "H1_sil_any", "LR_any",
                "H1_sil_fisher", "H1_sil_tippett",
                "H0_sil_fishercal", "H1_sil_fishercal", "LR_fishercal"]].to_string(index=False))
    print("\nWritten:", os.path.join(OUT, "omnibus_combination.csv"))


if __name__ == "__main__":
    main()
