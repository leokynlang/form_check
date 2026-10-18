"""Frozen v2 H0-only robust affine scores; no empirical-CDF tail ceiling."""
from __future__ import annotations

import numpy as np

from experiments.release_validation.calibration import choose_threshold, valid_method

IQR_NORMALIZER = 1.3489795003921634
MAD_NORMALIZER = 1.482602218505602
P_FLOOR = np.finfo(np.float64).tiny


def transform_raw(raw, method):
    raw = np.asarray(raw, dtype=float)
    if method in "ABCD":
        p = -raw
        valid = np.isfinite(p) & (p >= 0) & (p <= 1)
        transformed = np.full(raw.shape, np.nan)
        transformed[valid] = -np.log10(np.maximum(p[valid], P_FLOOR))
        return transformed
    if method not in "EFG":
        raise ValueError(f"Unknown method {method}")
    return np.where(np.isfinite(raw), raw, np.nan)


def _guard(frame, phase):
    if not frame.klass.eq("H0").all() or not frame.phase.eq(phase).all():
        raise ValueError(f"{phase} requires H0 in the matching phase only")
    if "calibration_round" in frame and frame.calibration_round.nunique() > 1:
        raise ValueError("Do not pool independent calibration rounds")


def fit_reference(frame, methods):
    _guard(frame, "reference")
    reference = {}
    for method in methods:
        raw = frame[f"score_{method}"].to_numpy(dtype=float)
        transformed = transform_raw(raw, method)
        values = transformed[valid_method(frame, method) & np.isfinite(transformed)]
        record = {"usable": False, "center": None, "scale": None,
                  "scale_method": "unusable", "n": len(values),
                  "raw_transform": "negative_log10_p" if method in "ABCD" else "identity",
                  "reason": "empty_reference" if not len(values) else "degenerate_reference_scale",
                  "reference_p_zero_count": int(np.sum(raw == 0)) if method in "ABCD" else 0,
                  "reference_p_floor_count": int(np.sum(raw == -.001)) if method == "A" else 0,
                  "floating_floor_count": int(np.sum((-raw >= 0) & (-raw <= P_FLOOR))) if method in "ABCD" else 0}
        if len(values):
            center = float(np.median(values))
            quartiles = np.quantile(values, [.25, .75], method="linear")
            scales = [("iqr_normalized", float((quartiles[1] - quartiles[0]) / IQR_NORMALIZER)),
                      ("mad_normalized", float(np.median(np.abs(values - center)) * MAD_NORMALIZER)),
                      ("population_sd", float(np.std(values, ddof=0)))]
            record["center"] = center
            for label, scale in scales:
                if np.isfinite(scale) and scale > 0:
                    record.update({"usable": True, "scale": scale, "scale_method": label, "reason": None})
                    break
        reference[method] = record
    return reference


def policy_scores(frame, members, reference):
    if not members:
        raise ValueError("A policy must have a member")
    valid = np.ones(len(frame), dtype=bool)
    columns = []
    for method in members:
        ref = reference[method]
        transformed = transform_raw(frame[f"score_{method}"].to_numpy(dtype=float), method)
        ok = valid_method(frame, method) & np.isfinite(transformed)
        if ref["usable"]:
            with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                score = (transformed - ref["center"]) / ref["scale"]
            ok &= np.isfinite(score)
        else:
            ok[:] = False
            score = np.full(len(frame), np.nan)
        columns.append(score)
        valid &= ok
    score = np.max(np.column_stack(columns), axis=1)
    score[~valid] = np.nan
    return score, ~valid


def calibrate(frame, members, reference, budget):
    _guard(frame, "calibration")
    score, decline = policy_scores(frame, members, reference)
    return choose_threshold(score, decline, budget)


def apply_policy(frame, members, reference, fitted):
    if not fitted["feasible"]:
        raise ValueError("Infeasible policy cannot be treated as budget-compliant")
    score, decline = policy_scores(frame, members, reference)
    alarm = ~decline & (score > fitted["threshold"])
    return score, alarm, decline
