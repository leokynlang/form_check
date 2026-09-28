"""Review-only numerical audit; does not alter benchmark/product or manuscript.

Same null/added subspaces, orthonormal coordinates, HC3 sandwich directly on Q.
This small independent-seed probe is not a replacement for the full benchmark.
"""
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import common
import srsd_data


def stable_test(x, y, variant):
    z = (x-x.mean(0))/x.std(0)
    base = np.column_stack([np.ones(len(y)), z])
    qb, _ = np.linalg.qr(base, mode="reduced")
    yp = qb @ (qb.T @ y)
    if variant == "fitted":
        u = (yp-yp.mean())/yp.std()
        extra = np.column_stack([u**2, u**3])
    else:
        extra = np.column_stack([z[:, j]**d for j in range(z.shape[1]) for d in [2, 3]])
    residual_extra = extra-qb @ (qb.T @ extra)
    if np.linalg.matrix_rank(residual_extra) != extra.shape[1]:
        raise ValueError("rank deficient added subspace")
    qe, _ = np.linalg.qr(residual_extra, mode="reduced")
    design = np.column_stack([qb, qe])
    assert np.allclose(design.T @ design, np.eye(design.shape[1]), atol=1e-9)
    coef = design.T @ y
    residual = y-design @ coef
    h = (design**2).sum(1)
    cov = design.T @ (design*(residual/(1-h))[:, None]**2)
    q = qe.shape[1]
    w = coef[-q:] @ np.linalg.solve(cov[-q:, -q:], coef[-q:])
    p = stats.f.sf(w/q, q, len(y)-design.shape[1])
    rss_u = residual @ residual
    rss_r = np.sum((y-yp)**2)
    f = ((rss_r-rss_u)/q)/(rss_u/(len(y)-design.shape[1]))
    return float(p), float(stats.f.sf(f, q, len(y)-design.shape[1]))


def main():
    rows = []
    tasks = [t for t in srsd_data.registry_with_meta() if t["klass"] == "H0"]
    for task in tasks:
        fn = srsd_data.make_dataset_fn(task["cls"])
        for seed in range(810000, 810050):
            df, features = fn(seed, 200, .03)
            x, y = np.log(df[features.split(",")].values), np.log(df.target.values)
            for variant in ["fitted", "exog"]:
                old = common.reset_test(x, y, variant=variant, vcov="hc3")
                stable, classical = stable_test(x, y, variant)
                rows.append({"task": task["eq_name"], "seed": seed, "variant": variant,
                             "original_p": old["p"], "original_failure": old["fail"],
                             "condition": old["cond"], "qr_hc3_p": stable,
                             "qr_classical_p": classical})
        print(task["eq_name"], "done", flush=True)
    out = ROOT / "runs" / "paper_v3_review"
    out.mkdir(exist_ok=True)
    data = pd.DataFrame(rows)
    data.to_csv(out/"hc3_coordinate_probe.csv", index=False)
    summary = {}
    for variant, g in data.groupby("variant"):
        summary[variant] = {"n": len(g), "original_failures": int(g.original_failure.notna().sum()),
                            "original_rejection_rate": float((g.original_p < .05).mean()),
                            "qr_hc3_rejection_rate": float((g.qr_hc3_p < .05).mean()),
                            "qr_classical_rejection_rate": float((g.qr_classical_p < .05).mean()),
                            "verdict_changes": int(((g.original_p < .05) != (g.qr_hc3_p < .05)).sum())}
    summary["scope"] = "13 H0 SRSD tasks, 50 new seeds, n=200, log noise=.03; numerical probe only"
    (out/"hc3_coordinate_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    # Existing S5: actual OR false-alarm cost, no refits or new data.
    data = pd.read_csv(ROOT/"runs/gate_cmp/results.csv")
    gate = json.loads((ROOT/"runs/gate_cmp/summary.json").read_text())
    rows = []
    for (bench, n), g in data.groupby(["bench", "n"]):
        h0 = g[g.klass == "H0"]
        threshold = gate["cells"][f"{bench}_n{n}"]["naive"]["d_star"]
        rows.append({"bench": bench, "n": n, "n_h0": len(h0),
                     "fc_fpr": float(h0.rej_fc.mean()),
                     "cv_fpr": float((h0.D_cv > threshold).mean()),
                     "or_fpr": float((h0.rej_fc | (h0.D_cv > threshold)).mean())})
    pd.DataFrame(rows).to_csv(out/"gate_actual_fpr.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    # S6 nested models must use the SAME countries to isolate added-feature effects.
    raw = ROOT/"runs/wb_real/raw"
    countries = json.loads((raw/"countries.json").read_text(encoding="utf-8"))[1]
    allowed = {r["id"] for r in countries if r["region"]["value"] != "Aggregates"}
    compare = []
    for year in [2018, 2019, 2022]:
        cols = {}
        for code, label in [("AG.YLD.CREL.KG", "yield"), ("AG.CON.FERT.ZS", "fert"),
                            ("EG.USE.PCAP.KG.OE", "energy")]:
            records = json.loads((raw/f"{code}_{year}.json").read_text(encoding="utf-8"))[1]
            cols[label] = {r["countryiso3code"]: float(r["value"]) for r in records
                           if r["value"] is not None and r["countryiso3code"] in allowed}
        df = pd.DataFrame(cols).dropna()
        df = df[(df > 0).all(axis=1)]
        for seed in [11, 22, 33, 44, 55]:
            for features in [["fert"], ["fert", "energy"]]:
                result = common.run_all_methods(df, "yield", features, seed)
                compare.append({"year": year, "n": len(df), "seed": seed,
                                "features": ",".join(features), **result})
    frame = pd.DataFrame(compare)
    frame.to_csv(out/"wb_same_countries_probe.csv", index=False)
    print(frame.groupby(["year", "features"])[["n", "p_fc", "p_fitted_f"]].median().to_string())


if __name__ == "__main__":
    main()
