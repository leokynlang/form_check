"""Read-only aggregation of cached S2 results; write ONLY beside this script.

No data generation, fitting, permutation testing, or modification of original
artifacts. Figure layout adapts paper/review3_analysis.py D2/D3/D4, but all
panels here use the explicitly documented sentinel-valid pools and decisions.
Run: python audit_results.py
"""
from __future__ import annotations

import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
assert OUT == ROOT / "paper" / "research_R2", "Unexpected audit output directory"
KEYS = ["task", "n", "rel", "seed"]
NS = [50, 200, 1000]
METHODS = {
    "A": ("p_fc", "rej_fc", "fail_fc"),
    "B": ("p_fitted_f", "rej_fitted_f", "fail_fitted_f"),
    "C": ("p_C_svd", "rej_C_svd", "fail_C_svd"),
    "D": ("p_D_svd", "rej_D_svd", "fail_D_svd"),
}
INPUTS = {
    "raw_s2": ROOT / "runs/srsd_full/raw_srsd.parquet",
    "quality_s2": ROOT / "runs/srsd_full/quality_recompute.parquet",
    "svd_s2": ROOT / "runs/cd_svd/s2_cd_svd.csv",
    "s2_config": ROOT / "runs/srsd_full/config.json",
}
REFERENCE_SCRIPTS = [ROOT / "paper/review3_analysis.py",
                     ROOT / "paper/v5_round2_analysis.py",
                     ROOT / "svd_main_results.py", ROOT / "run_srsd.py",
                     ROOT / "srsd_fitting_quality.py"]
REPS = 10000


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def ratio(a, b):
    return float(a / b) if b else None


def present(s):
    """A typed failure string denotes decline; NaN/empty denotes no failure."""
    return s.notna() & s.astype("string").str.strip().ne("").fillna(False)


def bootstrap(values, rng, reps=REPS):
    values = np.asarray(values, dtype=float)
    assert np.isfinite(values).all() and len(values)
    draws = rng.integers(0, len(values), size=(reps, len(values)))
    means = values[draws].mean(axis=1)
    return np.quantile(means, [.025, .975]).tolist()


def wilson(k, n):
    if not n:
        return None
    z = 1.959963984540054
    p = k / n
    den = 1 + z * z / n
    mid = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [float(mid - half), float(mid + half)]


def load():
    raw = pd.read_parquet(INPUTS["raw_s2"])
    quality = pd.read_parquet(INPUTS["quality_s2"])
    svd = pd.read_csv(INPUTS["svd_s2"], low_memory=False)
    m = raw.merge(quality[KEYS + ["r2"]], on=KEYS, how="left",
                  validate="one_to_one", indicator="quality_merge")
    svd_cols = ["sentinel_ok"] + [x for method in ("C", "D") for x in METHODS[method]]
    m = m.merge(svd[KEYS + svd_cols], on=KEYS, how="left",
                validate="one_to_one", indicator="svd_merge")
    assert m.quality_merge.eq("both").all() and m.svd_merge.eq("both").all()
    for method, (p, rej, fail) in METHODS.items():
        m[f"valid_{method}"] = np.isfinite(m[p]) & ~present(m[fail])
        m[f"alarm_{method}"] = m[f"valid_{method}"] & m[rej].fillna(False).astype(bool)
    m["panel_valid"] = m[[f"valid_{method}" for method in METHODS]].all(axis=1)
    m["silence"] = m.panel_valid & ~m[[f"alarm_{method}" for method in METHODS]].any(axis=1)
    m["screen"] = np.isfinite(m.r2) & (m.r2 >= .99)
    m["certificate"] = m.screen & m.silence
    m["danger"] = np.isfinite(m.extrap_far_relerr) & (m.extrap_far_relerr > .10)
    return m


def describe_pool(df):
    return {"N": len(df), "tasks": sorted(df.task.unique().tolist()),
            "task_count": int(df.task.nunique()),
            "by_n": {str(n): int((df.n == n).sum()) for n in NS}}


def power_results(valid):
    h0 = valid[(valid.klass == "H0") & (valid.rel > 0)]
    h1 = valid[(valid.klass == "H1") & (valid.rel >= .005)]
    assert h0.task.nunique() == 13 and h1.task.nunique() == 18
    tasks = sorted(h1.task.unique())
    out = {"H0_pool": describe_pool(h0), "H1_pool": describe_pool(h1),
           "denominator_policy": "All H1 rows retained. Nonfinite p or typed decline is not a rejection; decline is separately counted, never certified as acceptance.",
           "calibration": "Per method and n: numpy quantile(method=linear) at .05 of finite, non-declined H0 p. Reject iff valid p <= threshold. Thresholds are not retuned within task bootstrap.",
           "bootstrap": {"repetitions": REPS, "paired_D_A_seed": 20260908,
                         "pointwise_power_seed": 20260922,
                         "unit": "equation, equally weighted; paired resampling for D-A"},
           "by_n": {}, "per_task": {}, "paired_D_A": {}}
    power_rng = np.random.default_rng(20260922)
    paired_rng = np.random.default_rng(20260908)
    for n in NS:
        hn, sn = h0[h0.n == n], h1[h1.n == n]
        out["by_n"][str(n)] = {}
        per = {}
        for method, (p, _, fail) in METHODS.items():
            calibration = hn.loc[hn[f"valid_{method}"], p]
            threshold = float(np.quantile(calibration, .05, method="linear"))
            reject = sn[f"valid_{method}"] & (sn[p] <= threshold)
            calibration_reject = hn[f"valid_{method}"] & (hn[p] <= threshold)
            per[method] = reject.groupby(sn.task).mean().reindex(tasks).to_numpy()
            declines = ~sn[f"valid_{method}"]
            out["by_n"][str(n)][method] = {
                "threshold": threshold, "calibration_N_all": len(hn),
                "calibration_N_valid": len(calibration),
                "calibration_FPR_all_rows": float(calibration_reject.mean()),
                "N_all": len(sn), "N_valid": int((~declines).sum()),
                "declines": int(declines.sum()), "decline_rate": float(declines.mean()),
                "decline_reasons": {str(k): int(v) for k, v in sn.loc[declines, fail].fillna("nonfinite_p_without_typed_reason").value_counts().items()},
                "rejections": int(reject.sum()), "power": float(reject.mean()),
                "power_valid_only_sensitivity": ratio(int(reject.sum()), int((~declines).sum())),
                "task_bootstrap_ci95": bootstrap(per[method], power_rng),
            }
        out["per_task"][str(n)] = {t: {method: float(per[method][i]) for method in METHODS}
                                     for i, t in enumerate(tasks)}
        diff = per["D"] - per["A"]
        out["paired_D_A"][str(n)] = {
            "mean_D": float(per["D"].mean()), "mean_A": float(per["A"].mean()),
            "delta": float(diff.mean()), "ci95": bootstrap(diff, paired_rng),
            "tasks": len(tasks), "D_wins": int((diff > 0).sum()),
            "ties": int((diff == 0).sum()), "D_loses": int((diff < 0).sum()),
        }
    return out


def h0_bootstrap(valid):
    h0 = valid[(valid.klass == "H0") & (valid.rel > 0)]
    tasks = sorted(h0.task.unique())
    rng = np.random.default_rng(20260920)
    out = {"seed": 20260920, "repetitions": REPS,
           "order": "n=50,200,1000; C SVD, C legacy, D SVD, D legacy (matches review3 RNG traversal)",
           "unit": "13 equations; 2000 rows/equation/n; native alpha=.05, no empirical threshold recalibration", "by_n": {}}
    for n in NS:
        sn = h0[h0.n == n]
        out["by_n"][str(n)] = {}
        for method, legacy in [("C", "rej_fitted_hc3"), ("D", "rej_exog_hc3")]:
            for path, column in [("SVD", f"alarm_{method}"), ("legacy", legacy)]:
                rates = sn[column].fillna(False).groupby(sn.task).mean().reindex(tasks).to_numpy()
                ci = bootstrap(rates, rng)
                loo = [float(np.delete(rates, i).mean()) for i in range(len(rates))]
                out["by_n"][str(n)][f"{method}_{path}"] = {
                    "N": len(sn), "task_count": len(tasks),
                    "point": float(rates.mean()), "ci95": ci,
                    "LOEO": [min(loo), max(loo)],
                    "fit_level_Wilson95": wilson(int(sn[column].fillna(False).sum()), len(sn)),
                }
    return out


def event_row(sn, threshold=.10, danger_column=None):
    danger = sn[danger_column] if danger_column else sn.extrap_far_relerr > threshold
    A = int(sn.screen.sum())
    B = int((sn.screen & danger).sum())
    C = int((sn.certificate & danger).sum())
    D = int(danger.sum())
    assert len(sn) >= A >= B >= C and len(sn) >= D >= B
    benign = ~danger
    return {"N_H1": len(sn), "task_count": int(sn.task.nunique()),
            "A_screened_wrong": A, "B_screened_dangerous": B,
            "C_released_dangerous": C, "D_all_dangerous": D,
            "gate_only_silent_dangerous": int((sn.silence & danger).sum()),
            "panel_declines": int((~sn.panel_valid).sum()),
            "C_over_A": ratio(C, A), "C_over_B": ratio(C, B), "C_over_D": ratio(C, D),
            "danger_share_B_over_A": ratio(B, A),
            "screened_nondangerous_wrong": int((sn.screen & benign).sum()),
            "released_nondangerous_wrong": int((sn.certificate & benign).sum()),
            "release_given_screened_nondangerous_wrong": ratio(int((sn.certificate & benign).sum()), int((sn.screen & benign).sum()))}


def extrap_results(valid):
    # Failure of a gate is retained as decline; only quantities needed to define
    # the population (r2 and far-field reference error) must be finite here.
    all_ex = valid[(valid.rel > 0) & np.isfinite(valid.r2) & np.isfinite(valid.extrap_far_relerr)].copy()
    h1 = all_ex[all_ex.klass == "H1"].copy()
    assert h1.task.nunique() == 17
    out = {"pool": describe_pool(h1), "all_classes_pool": describe_pool(all_ex),
           "definitions": {"screen": "finite r2 >= .99", "danger": "finite extrap_far_relerr > threshold",
                           "certificate": "screen AND all four gate p finite AND no typed decline AND no native-alpha rejection",
                           "declines": "retained in population denominator but never accepted"},
           "table15": {str(n): event_row(h1[h1.n == n]) for n in NS},
           "pooled": event_row(h1),
           "threshold_sensitivity": {str(t): {str(n): event_row(h1[h1.n == n], t) for n in NS}
                                     for t in [.02, .05, .10, .20, .50]}}
    index = ["task", "rel", "seed"]
    pivot = h1.pivot(index=index, columns="n", values="extrap_far_relerr")
    full = pivot[NS].dropna()
    transitions = {}
    for lo, hi in [(50, 200), (200, 1000), (50, 1000)]:
        a, b = full[lo] > .1, full[hi] > .1
        transitions[f"{lo}_to_{hi}"] = {
            "both": int((a & b).sum()), "only_lower_n": int((a & ~b).sum()),
            "only_higher_n": int((~a & b).sum()), "neither": int((~a & ~b).sum())}
    fixed = h1.merge((full[1000] > .1).rename("danger_fixed"),
                     left_on=index, right_index=True, how="inner", validate="many_to_one")
    out["fixed_label_n1000"] = {
        "complete_triples": len(full), "available_triples": len(pivot),
        "label_transitions": transitions,
        "by_n": {str(n): event_row(fixed[fixed.n == n], danger_column="danger_fixed") for n in NS},
        "interpretation": "Descriptive same-index comparison only. Holding the n=1000 label fixed does not isolate a causal screen effect from a gate effect or establish common nested samples. The first-feature probe coordinate also depends on the observed training maximum and can change across n."}
    true = all_ex.klass.eq("H0")
    danger_wrong = all_ex.klass.eq("H1") & all_ex.danger
    benign_wrong = all_ex.klass.eq("H1") & ~all_ex.danger
    states = {}
    for label, mask in [("true_H0", true), ("dangerous_wrong_H1", danger_wrong),
                        ("nondangerous_wrong_H1", benign_wrong)]:
        states[label] = {"N": int(mask.sum()), "certificates": int((mask & all_ex.certificate).sum()),
                         "P_certificate_given_state": ratio(int((mask & all_ex.certificate).sum()), int(mask.sum()))}
    p0 = states["true_H0"]["P_certificate_given_state"]
    pdanger = states["dangerous_wrong_H1"]["P_certificate_given_state"]
    out["certificate"] = {
        "states": states,
        "LR_dangerous_wrong_vs_true": pdanger / p0,
        "inverse_LR_true_vs_dangerous_wrong": p0 / pdanger,
        "empirical_three_state_P_dangerous_wrong_given_certificate": ratio(int((danger_wrong & all_ex.certificate).sum()), int(all_ex.certificate.sum())),
        "empirical_three_state_certificates": int(all_ex.certificate.sum()),
        "binary_illustration_equal_prior_P_dangerous_wrong_given_certificate": pdanger / (pdanger + p0),
        "binary_empirical_restricted_to_true_or_dangerous_wrong": ratio(int((danger_wrong & all_ex.certificate).sum()), int(((true | danger_wrong) & all_ex.certificate).sum())),
        "warning": "Equal-prior posterior is only a hypothetical two-state population excluding nondangerous H1. Do not insert prior P(H1)=19/32 into likelihood P(certificate|dangerous H1): the conditioned classes differ. Empirical full-pool posterior uses all three states."}
    table7 = {"pool": describe_pool(all_ex), "danger_definition": "far-field error>.10 across H0 and H1, no fit screen",
              "N_dangerous": int(all_ex.danger.sum()), "base_rate": float(all_ex.danger.mean()), "methods": {}}
    for method in METHODS:
        alarm = all_ex[f"alarm_{method}"]
        tp = int((alarm & all_ex.danger).sum())
        table7["methods"][method] = {
            "alarms": int(alarm.sum()), "true_positives": tp,
            "precision": ratio(tp, int(alarm.sum())), "recall": ratio(tp, int(all_ex.danger.sum())),
            "declines": int((~all_ex[f"valid_{method}"]).sum())}
    a, b, d = all_ex.alarm_A, all_ex.alarm_B, all_ex.danger
    table7["true_positive_overlap_A_B"] = {
        "A_only": int((a & ~b & d).sum()), "B_only": int((~a & b & d).sum()),
        "both": int((a & b & d).sum())}
    out["table7"] = table7
    out["old_results_status"] = "Old Table15, threshold sweep, fixed-label decomposition, Table7 and D4 must not be carried over unchanged: old pools included i.30.5. This audit supplies replacements. Raw danger flips are descriptive, not identified estimation-error mechanisms."
    return out, all_ex


def make_figures(power, extrap, all_ex):
    figdir = OUT / "figures"
    figdir.mkdir(exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    paths = {}
    tasks = list(power["per_task"]["200"])
    matrix = np.array([[power["per_task"]["200"][t][m] for m in METHODS] for t in tasks])
    fig, ax = plt.subplots(figsize=(7.2, 9.4))
    im = ax.imshow(matrix, aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
    ax.set_xticks(range(4), ["A", "B", "C (SVD)", "D (SVD)"])
    ax.set_yticks(range(len(tasks)), [t.replace("feynman-", "") for t in tasks])
    for i in range(len(tasks)):
        for j in range(4):
            ax.text(j, i, f"{matrix[i,j]:.3f}", ha="center", va="center",
                    color="white" if matrix[i,j] > .6 else "black", fontsize=8)
    ax.set_title("S2 calibrated detection: 18 valid H1 equations\nn = 200; 2,000 cached datasets per equation", pad=12)
    fig.colorbar(im, ax=ax, fraction=.045, label="Rejections / all H1 datasets")
    fig.text(.5, .018, "i.30.5 excluded; i.26.2 retained (gate analysis). H0-pool thresholds.\nDeclines remain in the denominator and are not acceptances.", ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .06, 1, 1))
    p = figdir / "review3_fig_D2.png"; fig.savefig(p, dpi=170); plt.close(fig); paths["D2"] = str(p)
    fig, ax = plt.subplots(figsize=(8, 5.4))
    for method, marker in zip(METHODS, ["o", "s", "^", "v"]):
        cells = [power["by_n"][str(n)][method] for n in NS]
        y = np.array([c["power"] for c in cells]); ci = np.array([c["task_bootstrap_ci95"] for c in cells])
        ax.errorbar(NS, y, yerr=[y-ci[:,0],ci[:,1]-y], marker=marker, capsize=4,
                    label=method + (" (SVD)" if method in "CD" else ""))
    ax.set_xscale("log"); ax.set_xticks(NS, [str(n) for n in NS]); ax.set_ylim(0, 1)
    ax.set_xlabel("Dataset sample size n"); ax.set_ylabel("Rejections / all H1 datasets")
    ax.set_title("S2: 18 valid H1 equations\nH0-pool calibrated detection, equation-bootstrap 95% intervals")
    ax.grid(alpha=.2); ax.legend(loc="lower right")
    fig.text(.5, .015, "10,000 resamples; thresholds held fixed. Pointwise, not simultaneous, intervals.\nAll-row denominator; declines counted separately. No new model fits.", ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .075, 1, 1))
    p = figdir / "review3_fig_D3.png"; fig.savefig(p, dpi=170); plt.close(fig); paths["D3"] = str(p)
    t7 = extrap["table7"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    xx = np.arange(4)
    axes[0].bar(xx-.18, [t7["methods"][m]["precision"] for m in METHODS], .36, label="Precision")
    axes[0].bar(xx+.18, [t7["methods"][m]["recall"] for m in METHODS], .36, label="Recall")
    axes[0].set_xticks(xx, ["A", "B", "C (SVD)", "D (SVD)"])
    axes[0].axhline(t7["base_rate"], ls="--", color="black", lw=1, label=f"Base rate {t7['base_rate']:.3f}")
    axes[0].set_ylim(0,1); axes[0].legend(fontsize=8); axes[0].set_title("Native-alpha alarms for far-field error > 10%")
    err = np.sort(all_ex.extrap_far_relerr.to_numpy())
    axes[1].plot(err, np.arange(1,len(err)+1)/len(err), lw=1.5)
    axes[1].axvline(.1, ls="--", color="black", lw=1)
    axes[1].set_xscale("symlog", linthresh=.001); axes[1].set_ylim(0,1)
    axes[1].set_xlim(0, float(err[-1])*1.05)
    axes[1].set_xlabel("Far-field relative error (symlog scale)")
    axes[1].set_ylabel("Empirical cumulative fraction"); axes[1].set_title("All finite errors; no fit-quality screen")
    fig.suptitle("S2 valid extrapolation pool: 13 H0 + 17 H1 equations", fontsize=12)
    fig.text(.5,.012,"i.30.5 excluded (sentinel failure); i.26.2 excluded (missing extrapolation). Declines are not alarms.",ha="center",fontsize=9)
    fig.tight_layout(rect=(0,.045,1,.96))
    p=figdir / "review3_fig_D4.png";fig.savefig(p,dpi=170);plt.close(fig);paths["D4"]=str(p)
    return {k: {"path": v, "sha256": sha256(v)} for k,v in paths.items()}


def main():
    config = json.loads(INPUTS["s2_config"].read_text(encoding="utf-8"))
    assert config["far"] == 3.0
    m = load()
    valid = m[(m.sentinel_ok == True) & (m.task != "feynman-i.30.5")].copy()
    power = power_results(valid)
    extrap, all_ex = extrap_results(valid)
    exceptions = {}
    for task in ["feynman-i.30.5", "feynman-i.26.2"]:
        t = m[(m.task == task) & (m.rel > 0)]
        exceptions[task] = {"nonzero_noise_rows": len(t), "sentinel_valid": int(t.sentinel_ok.eq(True).sum()),
                            "finite_r2": int(np.isfinite(t.r2).sum()),
                            "finite_extrap": int(np.isfinite(t.extrap_far_relerr).sum()),
                            "all_gates_valid": int(t.panel_valid.sum()),
                            "rule": "Excluded from every inferential pool: replay sentinel failed." if task.endswith("i.30.5") else "Retained for gate-only Table3; excluded only from extrapolation/certificate/Table7 pools because reference extrapolation is missing."}
    result = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "audit_script": {"path": str(Path(__file__).resolve()), "sha256": sha256(__file__)},
        "inputs": {k: {"path": str(v), "sha256": sha256(v), "bytes": v.stat().st_size} for k,v in INPUTS.items()},
        "reference_scripts": [{"path":str(p),"sha256":sha256(p)} for p in REFERENCE_SCRIPTS],
        "runtime": {"python":platform.python_version(),"numpy":np.__version__,"pandas":pd.__version__,"matplotlib":matplotlib.__version__},
        "provenance": "Reaggregation of cached results only, no new datasets or model fits. Input hashes identify exact artifacts, not independent validation of their generation.",
        "join": {"keys":KEYS,"validation":"one_to_one","matched_rows":len(m),"unmatched_rows":0},
        "exclusions": exceptions,
        "far_field_definition": {
            "source": "run_srsd.py:31,55-73; runs/srsd_full/config.json far=3.0",
            "cached_column": "raw_srsd.parquet extrap_far_relerr; quality_recompute supplies only fit-quality fields",
            "probe_grid": "12 log-space points along first feature: 6 from observed log-min-ln(3) to log-min, 6 from observed log-max to log-max+ln(3); other log-features held at their training medians",
            "reported_far_metric": "abs(exp(log_prediction-log_noiseless_closed_form_truth)-1) at the LAST grid point only: first feature=3*observed training maximum",
            "not": "Not one decade, not an error maximum across the grid, not the grid-median error, not a guarantee over an extrapolation domain",
            "sampling_warning": "Training extrema and hence the evaluation coordinate can change across n and seed; cross-n label changes do not isolate estimation error alone",
            "same_distribution_test_n": config["n_test"],
            "seed_rule": config["seed_rule"],
        },
        "table3": power,
        "h0_bootstrap": h0_bootstrap(valid),
        "extrapolation": extrap,
        "limits": ["Calibration FPR is fitted on H0 reference pools, not a new independent test of deployed FPR.",
                   "Task bootstrap holds H0 thresholds fixed and does not include threshold-estimation uncertainty.",
                   "A missing/declined diagnostic never creates an accepted certificate.",
                   "A danger-label change across n is descriptive; it does not identify a causal mechanism."]}
    result["figures"] = make_figures(power, extrap, all_ex)
    evidence = OUT / "evidence"; evidence.mkdir(exist_ok=True)
    destination = evidence / "results.json"
    destination.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    print(json.dumps({"results":str(destination),"table15":extrap["table15"],
                      "pooled":extrap["pooled"],"certificate":extrap["certificate"],
                      "table7":extrap["table7"],"figures":result["figures"]},indent=2,ensure_ascii=False))


if __name__ == "__main__":
    main()
