"""joint_silence_full_toolbox.py — silence of the ENTIRE falsification toolbox.

Strongest form of the user's question (2026-09-09): when the fitted form is
wrong, how often does NOTHING in the toolbox object — not the four spec tests
(A-D) AND not the nested-rival model-comparison gates (CV-quad, AICc-quad)?
gate_cmp ran fc + CV + AICc on the same seeds, so the joint event is
computable per row.

Union event (toolbox objects) = rej_fc | cv_fire | aicc_fire
  (fc's own rej_fc column is the runner's registered alpha-.05 decision)
Silence = none of the three fires.

Thresholds for cv/aicc: calibrated d*/m* from the same bench x n H0 pool
(read from runs/gate_cmp/summary.json), and the naive 0 thresholds for
reference.

Output: runs/joint_silence/toolbox_silence.csv + appended print
"""
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "runs", "joint_silence")

g = pd.read_csv(os.path.join(HERE, "runs", "gate_cmp", "results.csv"))
summ = json.load(open(os.path.join(HERE, "runs", "gate_cmp", "summary.json")))

# calibrated thresholds per bench x n
TH = {}
for cellkey, v in summ["cells"].items():
    bench, n = cellkey.rsplit("_n", 1)
    TH[(bench, int(n))] = (v["naive"]["d_star"], v["naive"]["m_star"])

g["cv_fire_cal"] = [row.D_cv > TH[(row.bench, int(row.n))][0]
                    if np.isfinite(row.D_cv) else False
                    for row in g.itertuples()]
g["aicc_fire_cal"] = [row.delta_aicc > TH[(row.bench, int(row.n))][1]
                      if np.isfinite(row.delta_aicc) else False
                      for row in g.itertuples()]
g["fail_any"] = g[["fail_fc", "fail_mc"]].notna().any(axis=1)

struct = g[(g.klass == "H1") & (~g.group.str.startswith("noise"))]
h0 = g[g.klass == "H0"]

def silence_stats(gg, fire_cols):
    fire_any = gg[fire_cols].any(axis=1)
    decl = gg["fail_any"].values
    sil = (~fire_any) & (~decl)
    return float((~fire_any).mean()), float(sil.mean()), int(decl.sum())

rows = []
for (bench, n), gg in struct.groupby(["bench", "n"]):
    for lbl, cols in [
        ("alpha_naive", ["rej_fc", "d0_fire", "m0_fire"]),
        ("calibrated", ["rej_fc", "cv_fire_cal", "aicc_fire_cal"])]:
        s_raw, s_clean, ndecl = silence_stats(gg, cols)
        # per-task
        per = []
        for task, gt in gg.groupby("task"):
            f = gt[cols].any(axis=1)
            per.append({"task": task, "silence": float((~f).mean())})
        pdf = pd.DataFrame(per)
        rows.append({
            "bench": bench, "n": n, "threshold": lbl, "N": len(gg),
            "N_declined": ndecl,
            "toolbox_silence_pool": s_clean,
            "toolbox_silence_incl_declined": s_raw,
            "mean_task_silence": float(pdf["silence"].mean()),
            "max_task_silence": float(pdf["silence"].max()),
            "worst_task": pdf.loc[pdf["silence"].idxmax(), "task"],
        })
out = pd.DataFrame(rows)

# fc vs toolbox gap: how much does adding the model-comparison gates buy?
piv = out.pivot_table(index=["bench", "n"], columns="threshold",
                      values="toolbox_silence_pool")

# H0 reference: silence on true nulls should be ~95% by construction
h0_rows = []
for (bench, n), gg in h0.groupby(["bench", "n"]):
    for lbl, cols in [("alpha_naive", ["rej_fc", "d0_fire", "m0_fire"]),
                      ("calibrated", ["rej_fc", "cv_fire_cal", "aicc_fire_cal"])]:
        f = gg[cols].any(axis=1)
        h0_rows.append({"bench": bench, "n": n, "threshold": lbl,
                        "toolbox_silence_H0": float((~f).mean())})
h0df = pd.DataFrame(h0_rows)

out.to_csv(os.path.join(OUT, "toolbox_silence.csv"), index=False)
h0df.to_csv(os.path.join(OUT, "toolbox_silence_h0_reference.csv"), index=False)

pd.set_option("display.width", 250)
print("=== FULL TOOLBOX silence on wrong-form data (fc OR CV-quad OR AICc-quad) ===")
print(out.to_string(index=False))
print()
print("=== H0 reference (silence on true power law; ~0.95 = calibrated) ===")
print(h0df.to_string(index=False))
