"""H0-only complete-policy calibration, with declines included in the budget."""
from __future__ import annotations

import numpy as np


def valid_method(frame, method):
    score = frame[f"score_{method}"].to_numpy(dtype=float)
    fail = frame[f"fail_{method}"].fillna("").astype(str).str.strip().ne("").to_numpy()
    return np.isfinite(score) & ~fail


def fit_reference(frame, methods):
    if not frame.klass.eq("H0").all() or not frame.phase.eq("reference").all():
        raise ValueError("CDF reference must contain reference-phase H0 only")
    return {method: np.sort(frame.loc[valid_method(frame, method), f"score_{method}"]
                            .to_numpy(dtype=float)).tolist() for method in methods}


def policy_scores(frame, members, reference):
    scores = []
    valid = np.ones(len(frame), dtype=bool)
    for method in members:
        ref = np.asarray(reference[method], dtype=float)
        ok = valid_method(frame, method)
        if not len(ref):
            ok[:] = False
        raw = frame[f"score_{method}"].to_numpy(dtype=float)
        # (# reference values strictly below observation)/(M+1).
        # Ties are not made artificially more anomalous; no random tie breaking.
        scores.append(np.searchsorted(ref, raw, side="left") / (len(ref) + 1))
        valid &= ok
    score = np.max(np.column_stack(scores), axis=1)
    score[~valid] = np.nan
    return score, ~valid


def choose_threshold(scores, decline, budget):
    scores, decline = np.asarray(scores, dtype=float), np.asarray(decline, dtype=bool)
    if not 0 <= budget < 1:
        raise ValueError("Budget must be in [0,1)")
    if scores.shape != decline.shape or np.any(~np.isfinite(scores[~decline])):
        raise ValueError("Non-declined scores must be finite")
    n = len(scores)
    cap = int(np.floor(budget * n + 1e-10))
    failures = int(decline.sum())
    if not n or failures > cap:
        return {"feasible": False, "threshold": None, "n": n,
                "declines": failures, "allowed_withheld": cap,
                "reason": "empty_screened_H0" if not n else "declines_exceed_budget"}
    values = np.sort(scores[~decline])
    allowed_alarms = cap - failures
    threshold = (-1.0 if allowed_alarms >= len(values)
                 else float(values[len(values) - allowed_alarms - 1]))
    alarms = int(np.sum(values > threshold))
    if alarms + failures > cap:
        raise AssertionError("Calibration exceeded its complete-policy budget")
    return {"feasible": True, "threshold": threshold, "n": n,
            "declines": failures, "alarms": alarms, "allowed_withheld": cap,
            "withheld_rate": (alarms + failures) / n}


def calibrate(frame, members, reference, budget):
    if not frame.klass.eq("H0").all() or not frame.phase.eq("calibration").all():
        raise ValueError("Threshold calibration must contain calibration-phase H0 only")
    scores, decline = policy_scores(frame, members, reference)
    return choose_threshold(scores, decline, budget)


def apply_policy(frame, members, reference, fitted):
    if not fitted["feasible"]:
        raise ValueError("Infeasible policy cannot be presented as budget-compliant")
    scores, decline = policy_scores(frame, members, reference)
    alarm = ~decline & (scores > fitted["threshold"])
    return scores, alarm, decline
