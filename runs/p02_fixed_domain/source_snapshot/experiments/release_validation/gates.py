"""Shared-dataset adapters; reuse the registered R2 definitions explicitly."""
from __future__ import annotations

import hashlib
import time

import numpy as np

from common import form_check_with_p, reset_test
from gate_compare import gate_model_compare
from paper.revision_loc.rerun_cd_svd import reset_test_svd

METHODS = tuple("ABCDEFG")


def dataset_digest(xlog, ly):
    h = hashlib.sha256()
    for arr in (xlog, ly):
        h.update(np.asarray(arr.shape, dtype="<i8").tobytes())
        h.update(np.asarray(arr, dtype="<f8", order="C").tobytes())
    return h.hexdigest()


def rbf_cv_score(xlog, ly, seed, landmarks=32, ridge=0.01):
    """Fixed RBF residual rival, training-only preprocessing/landmarks in each fold.

    Affine terms are unpenalized; only kernel coefficients get n_train*ridge.
    This is a calibrated finite rival score, not an omnibus-test guarantee.
    Folds match E exactly. No response-based parameter/landmark selection.
    """
    n, p = xlog.shape
    folds = np.array_split(np.random.default_rng(int(seed) * 7919 + 13).permutation(n), 5)
    ss_base = ss_rival = 0.0
    for fold_no, heldout in enumerate(folds):
        mask = np.ones(n, dtype=bool)
        mask[heldout] = False
        xtr, xte = xlog[mask], xlog[heldout]
        mean, sd = xtr.mean(axis=0), xtr.std(axis=0)
        if np.any(sd <= np.finfo(float).eps):
            raise ValueError("constant_feature")
        ztr, zte = (xtr - mean) / sd, (xte - mean) / sd
        basetr = np.column_stack([np.ones(len(ztr)), ztr])
        basete = np.column_stack([np.ones(len(zte)), zte])
        rng = np.random.default_rng(int(seed) * 104729 + fold_no)
        centers = ztr[rng.choice(len(ztr), min(landmarks, len(ztr)), replace=False)]

        def kernel(z, fold_centers=centers):
            sq = np.sum((z[:, None, :] - fold_centers[None, :, :]) ** 2, axis=2)
            return np.exp(-sq / (2 * p))

        ktr, kte = kernel(ztr), kernel(zte)
        # Remove affine overlap before ridge; equivalent to unpenalized affine
        # terms in a joint affine + kernel fit, without penalizing the base law.
        base_coef = np.linalg.lstsq(basetr, ly[mask], rcond=None)[0]
        projection = np.linalg.lstsq(basetr, ktr, rcond=None)[0]
        kres_tr, kres_te = ktr - basetr @ projection, kte - basete @ projection
        count = len(centers)
        aug = np.vstack([kres_tr, np.sqrt(len(ztr) * ridge) * np.eye(count)])
        rhs = np.concatenate([ly[mask] - basetr @ base_coef, np.zeros(count)])
        weights = np.linalg.lstsq(aug, rhs, rcond=None)[0]
        pred_base = basete @ base_coef
        pred_rival = pred_base + kres_te @ weights
        ss_base += float(np.sum((ly[heldout] - pred_base) ** 2))
        ss_rival += float(np.sum((ly[heldout] - pred_rival) ** 2))
    score = np.sqrt(ss_base / n) - np.sqrt(ss_rival / n)
    if not np.isfinite(score):
        raise ValueError("nonfinite_statistic")
    return float(score)


def run_gates(data, features, seed):
    xlog, ly = np.log(data[features].to_numpy()), np.log(data.target.to_numpy())
    result = {"dataset_sha256": dataset_digest(xlog, ly)}
    for method in METHODS:
        result.update({f"score_{method}": np.nan, f"fail_{method}": None})

    def record(method, score, failure, elapsed):
        result[f"score_{method}"] = float(score)
        result[f"fail_{method}"] = failure or (None if np.isfinite(score) else "nonfinite_statistic")
        result[f"seconds_{method}"] = elapsed

    started = time.perf_counter()
    try:
        a = form_check_with_p(data, "target", features, seed, n_boot=200, b_perm=999)
        record("A", -a["p_fc"], a["fail_fc"], time.perf_counter() - started)
    except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
        record("A", np.nan, f"adapter:{type(exc).__name__}:{exc}", time.perf_counter() - started)
    for method, variant, robust in [("B", "fitted", False), ("C", "fitted", True), ("D", "exog", True)]:
        started = time.perf_counter()
        try:
            answer = (reset_test_svd(xlog, ly, variant=variant) if robust
                      else reset_test(xlog, ly, variant=variant, vcov="classical"))
            record(method, -answer["p"], answer["fail"], time.perf_counter() - started)
        except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
            record(method, np.nan, f"adapter:{type(exc).__name__}:{exc}", time.perf_counter() - started)
    started = time.perf_counter()
    mc = gate_model_compare(xlog, ly, seed)
    elapsed = time.perf_counter() - started
    for method, key in [("E", "D_cv"), ("F", "delta_aicc")]:
        record(method, mc[key], mc["fail_mc"], elapsed)
    result["seconds_EF_shared"] = elapsed
    started = time.perf_counter()
    try:
        record("G", rbf_cv_score(xlog, ly, seed), None, time.perf_counter() - started)
    except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
        record("G", np.nan, f"adapter:{type(exc).__name__}:{exc}", time.perf_counter() - started)
    return result
