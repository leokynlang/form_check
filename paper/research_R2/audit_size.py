"""Audit cached H0 size, without importing benchmark code or fitting models.

Writes only evidence/size.json and figures/review3_fig_D1.png beside this
script. A/B use original cached p-values; C/D use cached SVD p-values.
"""
from __future__ import annotations

import hashlib
import json
import math
import platform
from datetime import datetime, timezone
from pathlib import Path
from statistics import NormalDist

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent
SOURCE = Path(__file__).resolve().parents[2]
FILES = {
    "s1_raw": SOURCE / "runs/full_500/raw.parquet",
    "s1_svd": SOURCE / "runs/cd_svd/s1_cd_svd.csv",
    "s2_raw": SOURCE / "runs/srsd_full/raw_srsd.parquet",
    "s2_svd": SOURCE / "runs/cd_svd/s2_cd_svd.csv",
    "s1_old_summary": SOURCE / "runs/full_500/summary.json",
    "s1_svd_summary": SOURCE / "runs/cd_svd/summary_s1.json",
    "s1_runner": SOURCE / "run_experiment.py",
    "s1_old_aggregator": SOURCE / "analyze.py",
}
H0 = {"H0_BASELINE", "H0_MULTIVAR", "H0_CORR"}
METHODS = {
    "A": ("p_fc", "rej_fc", "fail_fc"),
    "B": ("p_fitted_f", "rej_fitted_f", "fail_fitted_f"),
    "C": ("p_C_svd", "rej_C_svd", "fail_C_svd"),
    "D": ("p_D_svd", "rej_D_svd", "fail_D_svd"),
}
ALPHA = 0.05
Z = NormalDist().inv_cdf(0.975)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def active_failure(column: pd.Series) -> pd.Series:
    return column.notna() & ~column.astype(str).str.lower().isin(
        ["", "false", "0", "none", "nan"])


def rate_record(rejected: int, total: int) -> dict:
    if not total:
        return {"rejected": rejected, "evaluated": 0, "rate": None,
                "wilson95": None, "display_4dp": "NA"}
    p = rejected / total
    d = 1 + Z * Z / total
    c = (p + Z * Z / (2 * total)) / d
    w = Z * math.sqrt(p * (1 - p) / total + Z * Z / (4 * total**2)) / d
    lo, hi = max(0., c - w), min(1., c + w)
    return {"rejected": rejected, "evaluated": total, "rate": p,
            "wilson95": [lo, hi],
            "display_4dp": f"{p:.4f} [{lo:.4f}-{hi:.4f}]"}


def load_pool(pool: str) -> tuple[pd.DataFrame, str]:
    raw = pd.read_parquet(FILES[f"{pool}_raw"])
    ident = "scenario" if pool == "s1" else "task"
    raw = raw[raw.scenario.isin(H0) if pool == "s1" else raw.klass.eq("H0")].copy()
    svd = pd.read_csv(FILES[f"{pool}_svd"], low_memory=False)
    keys = [ident, "n", "rel", "seed", "base_seed"]
    cols = keys + ["sentinel_ok", "p_C_svd", "rej_C_svd", "fail_C_svd",
                   "p_D_svd", "rej_D_svd", "fail_D_svd"]
    result = raw.merge(svd[cols], on=keys, how="left", validate="one_to_one",
                       indicator=True)
    assert result._merge.eq("both").all(), "Missing cached SVD matches"
    assert result.sentinel_ok.eq(True).all(), "H0 sentinel mismatch"
    return result, ident


def summarize(frame: pd.DataFrame, ident: str) -> tuple[dict, dict]:
    answer = {"h0_tasks": sorted(frame[ident].unique().tolist()), "by_n": {}}
    curves = {}
    for n, group in frame.groupby("n"):
        positive = group.rel.gt(0)
        data_ok = (~active_failure(group.data_fail) if "data_fail" in group
                   else pd.Series(True, index=group.index))
        cell = {
            "all_h0_rows": len(group),
            "zero_noise_rows_excluded": int(group.rel.eq(0).sum()),
            "positive_noise_rows": int(positive.sum()),
            "noise_levels_positive": sorted(group.loc[positive, "rel"].unique().tolist()),
            "grid_counts": [{"task": str(t), "rel": float(r), "rows": int(k)}
                            for (t, r), k in group.groupby([ident, "rel"]).size().items()],
            "methods": {},
        }
        for method, (pcol, rcol, fcol) in METHODS.items():
            p = pd.to_numeric(group[pcol], errors="coerce")
            finite = np.isfinite(p)
            ok = positive & data_ok & ~active_failure(group[fcol]) & finite & p.between(0, 1)
            observed = group.loc[ok, rcol].astype(bool)
            decision = p.loc[ok].lt(ALPHA)
            assert observed.eq(decision).all(), f"Stored rejection convention mismatch: {method}"
            record = rate_record(int(decision.sum()), int(ok.sum()))
            record.update({
                "positive_noise_excluded": int((positive & ~ok).sum()),
                "positive_noise_failure_flag": int((positive & active_failure(group[fcol])).sum()),
                "positive_noise_nonfinite_p": int((positive & ~finite).sum()),
                "zero_noise_failure_flag": int((group.rel.eq(0) & active_failure(group[fcol])).sum()),
                "rejected_if_failed_rows_counted_nonreject": int(group[rcol].astype(bool).sum()),
                "legacy_all_rows_denominator": len(group),
                "legacy_all_rows_rate": float(group[rcol].astype(bool).mean()),
                "p_equal_alpha": int(p.loc[ok].eq(ALPHA).sum()),
            })
            cell["methods"][method] = record
            curves[(int(n), method)] = np.sort(p.loc[ok].to_numpy())
        answer["by_n"][str(int(n))] = cell
    return answer, curves


def plot(curves: dict) -> Path:
    fig, axes = plt.subplots(2, 3, figsize=(12.0, 7.0), sharex=True, sharey=True)
    colors = ["#1565c0", "#d55e00", "#009e73", "#8e44ad"]
    for row, pool in enumerate(["s1", "s2"]):
        for col, n in enumerate([50, 200, 1000]):
            ax = axes[row, col]
            for method, color in zip(METHODS, colors):
                p = curves[pool][(n, method)]
                ax.step(p, np.arange(1, len(p) + 1) / len(p), where="post",
                        color=color, linewidth=1.25, label=method)
            ax.plot([0, 1], [0, 1], "k--", linewidth=.85, label="Uniform reference")
            ax.set_title(f"{pool.upper()} H0, n = {n}")
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.grid(alpha=.18)
            if row == 1:
                ax.set_xlabel("p-value")
            if col == 0:
                ax.set_ylabel("Empirical CDF")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5, bbox_to_anchor=(.5, .055))
    fig.suptitle("H0 p-values: positive noise, defined and failure-free evaluations", fontsize=13)
    fig.text(.5, .018, "A/B: original cached path; C/D: cached SVD path. "
             "ECDF uses p <= x; reported rejection rates use p < 0.05.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=[0, .105, 1, .95])
    output = BASE / "figures/review3_fig_D1.png"
    output.parent.mkdir(exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)
    return output


def main() -> None:
    evidence = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "operation": "Read-only aggregation of cached results; no simulation or fitting",
        "python": platform.python_version(), "pandas": pd.__version__,
        "numpy": np.__version__, "matplotlib": matplotlib.__version__,
        "alpha": ALPHA, "decision_rule": "p < alpha (asserted against stored reject flags)",
        "wilson_z": Z,
        "interval_scope": "Pointwise Wilson intervals for conditional simulation proportions; "
                          "not task-generalization or simultaneous coverage guarantees",
        "validity_filter": "H0, rel > 0, no data/method failure, finite p in [0,1], matching SVD sentinel",
        "inputs": {name: {"path": str(path), "sha256": digest(path), "bytes": path.stat().st_size}
                   for name, path in FILES.items()},
        "script_sha256": digest(Path(__file__)),
        "tables": {},
    }
    curves = {}
    for pool in ["s1", "s2"]:
        frame, ident = load_pool(pool)
        evidence["tables"][pool], curves[pool] = summarize(frame, ident)
    with FILES["s1_old_summary"].open(encoding="utf-8") as stream:
        legacy = json.load(stream)["h0_false_alarm_alpha05"]
    names = dict(A="fc", B="fitted_f", C="fitted_hc3", D="exog_hc3")
    for n, cell in evidence["tables"]["s1"]["by_n"].items():
        assert cell["positive_noise_rows"] == 6000
        assert cell["zero_noise_rows_excluded"] == 1500
        for method, row in cell["methods"].items():
            assert math.isclose(row["legacy_all_rows_rate"], legacy[n][names[method]]["rate"], abs_tol=1e-12)
    evidence["s1_legacy_check"] = "All 12 old Table 1 point rates reproduced only with zero-noise failures included in 7500 denominator; correct positive-noise denominator is 6000."
    figure = plot(curves)
    evidence["figure"] = {"path": str(figure), "sha256": digest(figure)}
    target = BASE / "evidence/size.json"
    target.parent.mkdir(exist_ok=True)
    target.write_text(json.dumps(evidence, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    for pool, table in evidence["tables"].items():
        print(pool.upper())
        for n, cell in table["by_n"].items():
            print(n, {m: f"{v['rejected']}/{v['evaluated']} {v['display_4dp']}"
                      for m, v in cell["methods"].items()})
    print(target)
    print(figure)


if __name__ == "__main__":
    main()
