"""Independent consistency checks on the saved real engineering evaluation."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from real_engineering_holdout import OUT, ROOT, digest, load_data, polynomial_fit
from dimreflex import fit_power_law


def main():
    manifest = json.loads((OUT / "manifest.json").read_text())
    assert manifest["complete"] and manifest["product_unchanged"]
    for name, sha in manifest["product_sha256"].items():
        assert digest(ROOT / name) == sha
    assert digest(Path(__file__).with_name("real_engineering_holdout.py")) == manifest["script_sha256"]
    # Source certificates are used ONLY here to verify parsing, after evaluation.
    # They were not used in any training/prediction or optimizer initialization.
    ch, feats = load_data("Chwirut1")
    x = ch.distance.to_numpy()
    yhat = np.exp(-.19027818370*x)/(.0061314004477+.010530908399*x)
    rss_ch = float(np.sum((ch.target-yhat)**2))
    assert np.isclose(rss_ch, 2.3844771393e3, rtol=1e-9)
    ne, _ = load_data("Nelson")
    yhat = 2.5906836021 - 5.6177717026e-9*ne.weeks*np.exp(.057701013174*ne.temperature_C)
    rss_ne = float(np.sum((np.log(ne.target)-yhat)**2))
    assert np.isclose(rss_ne, 3.7976833176, rtol=1e-9)
    # Independent product baseline equivalence (standardization must not change OLS).
    predict, _ = polynomial_fit(ch, feats, 1)
    fit = fit_power_law(ch, "target", feats)
    expected = np.log(fit["prefactor"]) + fit["exponents"]["distance"] * np.log(x)
    assert np.allclose(predict(ch), expected, atol=1e-11)
    splits = json.loads((OUT / "splits.json").read_text())
    assert len(splits) == 84
    for split in splits:
        assert not set(split["train_ids"]) & set(split["test_ids"])
        if split["split"] != "random_rows":
            assert split["overlap_groups"] == 0
        assert len(split["train_ids"]) + len(split["test_ids"]) == (214 if split["dataset"] == "Chwirut1" else 128)
    metrics = pd.read_csv(OUT / "metrics.csv")
    preds = pd.read_csv(OUT / "predictions.csv")
    assert len(metrics) == 546
    assert not metrics.duplicated(["dataset", "split", "seed", "model"]).any()
    checked = 0
    for key, rows in preds.groupby(["dataset", "split", "seed", "model", "partition"]):
        ds, split, seed, model, part = key
        match = metrics[(metrics.dataset == ds) & (metrics.split == split) &
                        (metrics.seed == seed) & (metrics.model == model)].iloc[0]
        assert len(rows) == match["n_"+part]
        err = rows.predicted_log - rows.observed_log
        rmse = float(np.sqrt(np.mean(err**2)))
        assert np.isclose(rmse, match[part+"_rmse_log"], rtol=1e-10, atol=1e-10)
        checked += 1
    checks = pd.read_csv(OUT / "checks.csv")
    assert len(checks) == 84
    for method in ["fc", "fitted_f", "fitted_hc3", "exog_hc3"]:
        good = checks[checks["fail_"+method].isna()]
        assert np.array_equal(good["p_"+method].to_numpy() < .05, good["rej_"+method])
    result = {"passed": True, "source_rss_chwirut1": rss_ch, "source_rss_nelson": rss_ne,
              "splits_verified": len(splits), "model_rows": len(metrics),
              "prediction_partitions_recomputed": checked,
              "failed_model_rows": int(metrics.failure.notna().sum()),
              "nonfinite_relative_error_rows": int((metrics.test_nonfinite_relative_errors > 0).sum()),
              "baseline_matches_product": True, "product_hashes_unchanged": True}
    (OUT / "audit.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
