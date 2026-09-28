"""Reproducible real-data validation; see real_data_protocol.md.

Run from any directory: python validation/real_engineering_holdout.py
Only writes runs/real_engineering_holdout; product modules are read-only.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import itertools
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import common
from dimreflex import CandidateModel, discover_auto, predict_candidate

OUT = ROOT / "runs" / "real_engineering_holdout"
URL_BASE = "https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/"
SEEDS = list(range(20260905, 20260925))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def safe_json(value):
    if isinstance(value, dict):
        return {str(k): safe_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe_json(v) for v in value]
    if isinstance(value, np.ndarray):
        return safe_json(value.tolist())
    if isinstance(value, (np.integer, np.bool_)):
        return value.item()
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if callable(value):
        return "<callable omitted; numeric model parameters retained>"
    return value


def load_data(name):
    path = OUT / "raw" / f"{name}.dat"
    content = path.read_text(encoding="ascii")
    lines = content.splitlines()
    start = next(i + 1 for i, line in enumerate(lines)
                 if line.startswith("Data:") and "y" in line
                 and "Response" not in line)
    values = np.loadtxt(io.StringIO("\n".join(lines[start:])))
    expected = {"Chwirut1": (214, 2), "Nelson": (128, 3)}[name]
    assert values.shape == expected, (name, values.shape)
    if name == "Chwirut1":
        df = pd.DataFrame(values, columns=["target", "distance"])
        features = ["distance"]
    else:
        df = pd.DataFrame(values, columns=["target", "weeks", "temperature_C"])
        df["temperature_K"] = df.temperature_C + 273.15
        features = ["weeks", "temperature_K"]
    assert np.isfinite(df.values).all() and (df.values > 0).all()
    df["row_id"] = np.arange(len(df))
    df["group"] = df.groupby(features, sort=True).ngroup()
    df.to_csv(OUT / f"{name}_clean.csv", index=False)
    return df, features


def split_data(df, name, seeds):
    ids = np.arange(len(df))
    groups = np.unique(df.group)
    for seed in seeds:
        rng = np.random.default_rng(seed)
        test = np.sort(rng.choice(ids, int(np.ceil(len(ids) * .25)), replace=False))
        yield "random_rows", seed, np.setdiff1d(ids, test), test
        rng = np.random.default_rng(seed)
        held = rng.choice(groups, int(np.ceil(len(groups) * .25)), replace=False)
        mask = df.group.isin(held).to_numpy()
        yield "heldout_groups", seed, ids[~mask], ids[mask]
    if name == "Chwirut1":
        levels = np.sort(df.distance.unique())
        upper = levels[int(np.floor(.75 * len(levels)))]
        mask = df.distance >= upper
        yield "high_distance", seeds[0], ids[~mask], ids[mask]
        middle = levels[int(np.floor(.25 * len(levels))):int(np.ceil(.75 * len(levels)))]
        mask = df.distance.isin(middle)
        yield "middle_to_tails", seeds[0], ids[mask], ids[~mask]
    else:
        levels = np.sort(df.weeks.unique())
        mask = df.weeks >= levels[int(np.floor(.75 * len(levels)))]
        yield "late_time", seeds[0], ids[~mask], ids[mask]
        mask = df.temperature_K == df.temperature_K.max()
        yield "high_temperature", seeds[0], ids[~mask], ids[mask]


def polynomial_fit(train, features, degree):
    lx = np.log(train[features].to_numpy())
    mu, sd = lx.mean(axis=0), lx.std(axis=0)
    if (sd == 0).any():
        raise ValueError("constant training feature")
    powers = [p for p in itertools.product(range(degree + 1), repeat=len(features))
              if sum(p) <= degree]
    powers.sort(key=lambda p: (sum(p), p))

    def design(data):
        z = (np.log(data[features].to_numpy()) - mu) / sd
        return np.column_stack([np.prod(z ** p, axis=1) for p in powers])

    A = design(train)
    coef, _, rank, _ = np.linalg.lstsq(A, np.log(train.target), rcond=None)
    if rank != A.shape[1]:
        raise ValueError(f"rank_deficient:{rank}/{A.shape[1]}")
    meta = {"degree": degree, "n_params": len(coef), "rank": int(rank),
            "condition": float(np.linalg.cond(A)), "coef": coef,
            "center": mu, "scale": sd, "powers": powers}
    return lambda data: design(data) @ coef, meta


def grouped_poly_select(train, features, degrees, seed):
    groups = np.unique(train.group)
    rng = np.random.default_rng(seed)
    folds = np.array_split(rng.permutation(groups), min(4, len(groups)))
    scores = {}
    for degree in degrees:
        squared_errors = []
        try:
            for valid_groups in folds:
                mask = train.group.isin(valid_groups)
                assert not set(train.loc[mask, "group"]) & set(train.loc[~mask, "group"])
                predict, _ = polynomial_fit(train.loc[~mask], features, degree)
                error = predict(train.loc[mask]) - np.log(train.loc[mask].target)
                if not np.isfinite(error).all():
                    raise ValueError("nonfinite inner prediction")
                squared_errors.extend(error ** 2)
            scores[degree] = {"mse": float(np.mean(squared_errors)), "failure": None}
        except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
            scores[degree] = {"mse": None, "failure": str(exc)}
    valid = [degree for degree in degrees if scores[degree]["mse"] is not None]
    if not valid:
        raise ValueError(f"all_grouped_cv_failed:{scores}")
    chosen = min(valid, key=lambda degree: (scores[degree]["mse"], degree))
    predict, meta = polynomial_fit(train, features, chosen)
    meta.update({"selected_degree": chosen, "group_cv": scores,
                 "inner_cv_rmse_log": np.sqrt(scores[chosen]["mse"])})
    return predict, meta


def reference_fit(train, name):
    """Fit source-specified form to training LOG responses; not NIST certification.

    Chwirut log form = a - b*x - log(1+c*x), b,c positive.
    Nelson log form = a - b*t*exp(c*(T_C-225)), b,c positive.
    Dimensionless parameters use known fixed units; no test-derived starts.
    """
    if name == "Chwirut1":
        def evaluate(theta, data):
            x = data.distance.to_numpy()
            return theta[0] - np.exp(theta[1])*x - np.log1p(np.exp(theta[2])*x)
        starts = [[float(np.log(train.target).max()), np.log(b), np.log(c)]
                  for b in [.01, .1, 1.] for c in [.1, 1., 10.]]
        upper_bounds = [20, 5, 5]
    else:
        def evaluate(theta, data):
            t = data.weeks.to_numpy()
            temperature = data.temperature_C.to_numpy()
            return theta[0] - np.exp(theta[1])*t*np.exp(np.exp(theta[2])*(temperature-225.))
        starts = [[float(np.log(train.target).max()), np.log(b), np.log(c)]
                  for b in [.0001, .01, .1] for c in [.005, .03, .1]]
        upper_bounds = [20, 5, 0]
    candidates = []
    for start in starts:
        fitted = least_squares(lambda theta: evaluate(theta, train)-np.log(train.target),
                               start, bounds=([-20, -25, -15], upper_bounds),
                               max_nfev=3000)
        if fitted.success and np.isfinite(fitted.fun).all():
            candidates.append(fitted)
    if not candidates:
        raise ValueError("reference_optimizer_failed")
    fitted = min(candidates, key=lambda item: np.sum(item.fun**2))
    return lambda data: evaluate(fitted.x, data), {
        "n_params": 3, "coef": fitted.x, "objective": "log-response SSE",
        "successful_starts": len(candidates), "optimizer_cost": float(fitted.cost)}


def auto_fit(train, features, seed):
    result = discover_auto(train, "target", features, random_state=seed)
    if not result["success"]:
        raise ValueError(result["reason"])
    model = CandidateModel.from_dict(result["final_model"])

    def predict(data):
        values = predict_candidate(model, data, features)
        if not np.isfinite(values).all() or (values <= 0).any():
            raise ValueError("nonfinite_or_nonpositive_auto_prediction")
        return np.log(values)

    rec = result["recommended"]
    return predict, {"n_params": model.n_params, "selected_family": rec["family"],
                     "final_family": model.family, "params": model.params,
                     "inner_cv_rmse_log": rec["validation"]["cv_rmse_log"],
                     "backend_adjustment": result["backend_adjustment"],
                     "assessments": result["assessments"]}


def wild_direction_check(train, features, seed, b=999):
    """Experimental HC3 wild sensitivity, NOT product form_check."""
    lx = np.log(train[features].to_numpy())
    ly = np.log(train.target.to_numpy())
    base = np.column_stack([np.ones(len(train)), lx])
    q = ((lx-lx.mean(axis=0))**2).mean(axis=1)
    A = np.column_stack([base, q])
    if np.linalg.matrix_rank(A) != A.shape[1]:
        return {"p_wild_sensitivity": None, "wild_failure": "rank_deficient"}
    pinv_base = np.linalg.pinv(base)
    pred = base @ pinv_base @ ly
    leverage = np.einsum("ij,ji->i", base, pinv_base)
    residual = (ly-pred) / np.maximum(1-leverage, 1e-8)
    direction = np.linalg.pinv(A)[-1]
    observed = float(direction @ ly)
    rng = np.random.default_rng(seed+919)
    null = (rng.choice([-1., 1.], size=(b, len(train))) * residual) @ direction
    return {"p_wild_sensitivity": float((1+np.sum(np.abs(null) >= abs(observed)))/(b+1)),
            "wild_failure": None}


def metrics(data, prediction):
    y = np.log(data.target.to_numpy())
    prediction = np.asarray(prediction)
    if not np.isfinite(prediction).all():
        raise ValueError("nonfinite log prediction")
    error = prediction-y
    with np.errstate(over="ignore", invalid="ignore"):
        relative = np.abs(np.expm1(error))
    sst = np.sum((y-y.mean())**2)
    return {"rmse_log": float(np.sqrt(np.mean(error**2))),
            "r2_log": float(1-np.sum(error**2)/sst) if sst > 0 else None,
            "median_ape": float(np.median(relative)),
            "p90_ape": float(np.quantile(relative, .9, method="inverted_cdf")),
            "nonfinite_relative_errors": int((~np.isfinite(relative)).sum())}


def write_summary(rows, checks):
    results = pd.DataFrame(rows)
    results.to_csv(OUT / "metrics.csv", index=False)
    pd.DataFrame(checks).to_csv(OUT / "checks.csv", index=False)
    summary = []
    for (dataset, split, model), g in results.groupby(["dataset", "split", "model"], sort=False):
        ok = g[g.failure.isna()]
        row = {"dataset": dataset, "split": split, "model": model,
               "runs": len(g), "successes": len(ok), "failures": len(g)-len(ok)}
        for key in ["train_rmse_log", "test_rmse_log", "train_r2_log", "test_r2_log",
                    "test_median_ape", "inner_cv_rmse_log"]:
            vals = pd.to_numeric(ok[key], errors="coerce").dropna() if key in ok else pd.Series(dtype=float)
            for q, label in [(.25, "q25"), (.5, "median"), (.75, "q75")]:
                row[key+"_"+label] = float(vals.quantile(q)) if len(vals) else None
        summary.append(row)
    pd.DataFrame(summary).to_csv(OUT / "summary.csv", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", type=int, default=20)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    product_files = sorted((ROOT / "dimreflex").glob("*.py")) + [ROOT / "common.py"]
    before = {str(p.relative_to(ROOT)): digest(p) for p in product_files}
    manifest = {"python": sys.version, "platform": platform.platform(),
                "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
                "seeds": SEEDS[:args.seeds], "script_sha256": digest(__file__),
                "protocol_sha256": digest(Path(__file__).with_name("real_data_protocol.md")),
                "product_sha256": before,
                "data": {name: {"url": URL_BASE+name+".dat", "sha256": digest(OUT/"raw"/(name+".dat"))}
                         for name in ["Chwirut1", "Nelson"]}}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    rows, checks, predictions, splits, models = [], [], [], [], []
    for name in ["Chwirut1", "Nelson"]:
        df, features = load_data(name)
        degrees = [1, 2, 3, 9] if name == "Chwirut1" else [1, 2, 3]
        for split, seed, itrain, itest in split_data(df, name, SEEDS[:args.seeds]):
            train, test = df.iloc[itrain].copy(), df.iloc[itest].copy()
            assert not set(train.row_id) & set(test.row_id)
            overlap = set(train.group) & set(test.group)
            if split != "random_rows":
                assert not overlap
            context = {"dataset": name, "split": split, "seed": seed,
                       "n_train": len(train), "n_test": len(test),
                       "train_groups": train.group.nunique(), "test_groups": test.group.nunique(),
                       "overlap_groups": len(overlap)}
            splits.append({**context, "train_ids": train.row_id.tolist(), "test_ids": test.row_id.tolist(),
                           "feature_ranges_train": {f: [float(train[f].min()), float(train[f].max())] for f in features},
                           "feature_ranges_test": {f: [float(test[f].min()), float(test[f].max())] for f in features}})
            check = {**context, **common.run_all_methods(train, "target", features, seed),
                     **wild_direction_check(train, features, seed)}
            checks.append(check)
            fitters = {("power_law" if d == 1 else f"log_poly_{d}"):
                       (lambda degree=d: polynomial_fit(train, features, degree)) for d in degrees}
            fitters.update({"group_cv_poly": lambda: grouped_poly_select(train, features, degrees, seed),
                            "discover_auto": lambda: auto_fit(train, features, seed),
                            "nist_form_refit": lambda: reference_fit(train, name)})
            for model_name, fit in fitters.items():
                row = {**context, "model": model_name, "failure": None}
                try:
                    predict, meta = fit()
                    models.append({**context, "model": model_name, "fit": meta})
                    row.update({k: meta[k] for k in ["n_params", "condition", "selected_degree",
                                                    "selected_family", "inner_cv_rmse_log"] if k in meta})
                    for label, data in [("train", train), ("test", test)]:
                        logpred = predict(data)
                        row.update({label+"_"+k: v for k, v in metrics(data, logpred).items()})
                        predictions.extend({"dataset": name, "split": split, "seed": seed,
                                            "model": model_name, "partition": label,
                                            "row_id": int(rid), "observed_log": float(obs),
                                            "predicted_log": float(pred)}
                                           for rid, obs, pred in zip(data.row_id, np.log(data.target), logpred, strict=True))
                except Exception as exc:
                    row["failure"] = f"{type(exc).__name__}: {exc}"
                rows.append(row)
            write_summary(rows, checks)
            print(f"{name} {split} seed={seed} n={len(train)}/{len(test)} "
                  f"p_fc={check['p_fc']:.3g} ({time.perf_counter()-start:.1f}s)", flush=True)
    pd.DataFrame(predictions).to_csv(OUT / "predictions.csv", index=False)
    (OUT / "splits.json").write_text(json.dumps(safe_json(splits), indent=2), encoding="utf-8")
    (OUT / "models.json").write_text(json.dumps(safe_json(models), indent=2), encoding="utf-8")
    assert before == {str(p.relative_to(ROOT)): digest(p) for p in product_files}, "product changed"
    manifest["wall_seconds"] = time.perf_counter()-start
    manifest["product_unchanged"] = True
    manifest["complete"] = True
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"COMPLETE: {len(rows)} model fits, {len(checks)} evaluations", flush=True)


if __name__ == "__main__":
    main()
