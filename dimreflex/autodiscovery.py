"""Validated autonomous discovery and auditable backend adjustment.

The existing function-family layer answers "what can fit these data?".  This
module adds the decision layer: deterministic k-fold validation, complexity
diagnostics, a one-standard-error preference for simpler models, and optional
backend adjustment for the selected power-law model.
"""


import numpy as np
import pandas as pd

from .anomaly import fit_with_guard
from .families import FAMILIES, CandidateModel


def predict_candidate(candidate: CandidateModel, data: pd.DataFrame,
                      features: list[str]) -> np.ndarray:
    """Evaluate a built-in candidate on new rows."""
    family = candidate.family
    params = candidate.params

    if family.startswith("power_law"):
        prediction = np.full(len(data), params["prefactor"], dtype=float)
        for feature, exponent in params["exponents"].items():
            prediction *= data[feature].to_numpy(dtype=float) ** exponent
        return prediction

    if family == "piecewise_power_law":
        feature = features[0]
        x = data[feature].to_numpy(dtype=float)
        segments = params["segments"]

        def segment_prediction(segment, mask):
            # 新 schema:exponents 字典(多特征);旧 schema:标量 exponent(单特征)
            exponents = segment.get("exponents")
            if exponents is not None:
                values = np.full(int(mask.sum()), segment["prefactor"],
                                 dtype=float)
                for feat, exponent in exponents.items():
                    values *= (data[feat].to_numpy(dtype=float)[mask]
                               ** exponent)
                return values
            return segment["prefactor"] * x[mask] ** segment["exponent"]

        if params["n_segments"] == 1:
            return segment_prediction(segments[0], np.ones(len(data), bool))
        breakpoint = params["breakpoints"][0]
        left, right = segments
        mask = x < breakpoint
        prediction = np.empty(len(data), dtype=float)
        prediction[mask] = segment_prediction(left, mask)
        prediction[~mask] = segment_prediction(right, ~mask)
        return prediction

    if family == "rational_log":
        feature = features[0]
        x = data[feature].to_numpy(dtype=float)
        L = np.log10(x) - params["log_center"]
        numerator = sum(
            coefficient * L ** power
            for power, coefficient in enumerate(params["num_coefs"])
        )
        denominator = 1.0 + sum(
            coefficient * L ** power
            for power, coefficient in enumerate(
                params["den_coefs"], start=1)
        )
        log_prediction = numerator / denominator
        # 多特征协变量项(旧 schema 无这些键,自然退化为单特征预测)
        for feat, coefficient, center in zip(params.get("extra_features") or [],
                                             params.get("extra_coefs") or [],
                                             params.get("extra_centers") or [],
                                             strict=True):
            log_prediction = log_prediction + coefficient * (
                np.log10(data[feat].to_numpy(dtype=float)) - center)
        return 10.0 ** log_prediction

    predictor = params.get("predict")
    if callable(predictor):
        return np.asarray(predictor(data), dtype=float)
    raise ValueError(f"函数族 '{family}' 没有跨样本预测接口,无法交叉验证")


def _validate_input(data: pd.DataFrame, target: str,
                    features: list[str]) -> str | None:
    columns = [target, *list(features)]
    missing = [column for column in columns if column not in data.columns]
    if missing:
        return f"数据缺少列 {missing}"
    if not features:
        return "至少需要一个特征列"
    try:
        values = data[columns].to_numpy(dtype=float)
    except (TypeError, ValueError):
        return "目标或特征包含非数值数据"
    if len(data) < 4:
        return "自主发现至少需要 4 个样本"
    if not np.isfinite(values).all():
        return "目标或特征包含 NaN/inf"
    if (values <= 0).any():
        return "当前自主函数族在 log 空间工作,要求目标和特征均为正"
    return None


def _folds(n_samples: int, n_splits: int,
           random_state: int) -> list[np.ndarray]:
    n_splits = min(max(2, int(n_splits)), n_samples)
    order = np.random.default_rng(random_state).permutation(n_samples)
    return [fold for fold in np.array_split(order, n_splits) if len(fold)]


def _information_criteria(y: np.ndarray, prediction: np.ndarray,
                          n_params: int) -> dict:
    log_y = np.log(y)
    log_prediction = np.log(prediction)
    residual = log_y - log_prediction
    n = len(y)
    sse = max(float(np.sum(residual ** 2)), np.finfo(float).tiny)
    k = max(1, int(n_params))
    aic = n * np.log(sse / n) + 2 * k
    correction = (2 * k * (k + 1) / (n - k - 1)
                  if n > k + 1 else float("inf"))
    return {
        "train_rmse_log": float(np.sqrt(sse / n)),
        "aic": float(aic),
        "aicc": float(aic + correction),
        "bic": float(n * np.log(sse / n) + k * np.log(n)),
    }


def assess_family(data: pd.DataFrame, target: str, features: list[str],
                  family: str, n_splits: int = 5,
                  random_state: int = 42) -> dict:
    """Fit and out-of-fold validate one registered family."""
    if family not in FAMILIES:
        return {"family": family, "success": False,
                "reason": f"未注册的函数族 '{family}'"}

    model = FAMILIES[family]().fit(data, target, features)
    assessment = {
        "family": family,
        "success": False,
        "model": model.to_dict(),
        "validation": {},
        "reason": model.reason,
    }
    if not model.success:
        return assessment

    n = len(data)
    folds = _folds(n, n_splits, random_state)
    predictions = np.full(n, np.nan)
    fold_r2 = []
    fold_failures = []
    all_indices = np.arange(n)
    for fold_index, valid in enumerate(folds):
        train = np.setdiff1d(all_indices, valid, assume_unique=True)
        fold_model = FAMILIES[family]().fit(
            data.iloc[train], target, features)
        if not fold_model.success:
            fold_failures.append(
                {"fold": fold_index, "reason": fold_model.reason})
            continue
        try:
            predicted = predict_candidate(
                fold_model, data.iloc[valid], features)
        except Exception as exc:  # noqa: BLE001 — 预测属函数族契约,任何崩溃都记为该折失败而非中断 CV
            fold_failures.append(
                {"fold": fold_index,
                 "reason": f"{type(exc).__name__}: {exc}"})
            continue
        if (not np.isfinite(predicted).all()) or (predicted <= 0).any():
            fold_failures.append(
                {"fold": fold_index, "reason": "预测含非有限值或非正值"})
            continue
        predictions[valid] = predicted
        y_fold = np.log(data.iloc[valid][target].to_numpy(dtype=float))
        p_fold = np.log(predicted)
        ss_res = float(np.sum((y_fold - p_fold) ** 2))
        ss_tot = float(np.sum((y_fold - y_fold.mean()) ** 2))
        fold_r2.append(1 - ss_res / ss_tot if ss_tot > 0 else 0.0)

    if fold_failures or not np.isfinite(predictions).all():
        assessment["reason"] = "交叉验证未完成"
        assessment["validation"] = {
            "n_splits": len(folds), "fold_failures": fold_failures}
        return assessment

    y = data[target].to_numpy(dtype=float)
    log_y = np.log(y)
    log_prediction = np.log(predictions)
    ss_res = float(np.sum((log_y - log_prediction) ** 2))
    ss_tot = float(np.sum((log_y - log_y.mean()) ** 2))
    cv_r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    cv_rmse = float(np.sqrt(ss_res / n))
    cv_se = (float(np.std(fold_r2, ddof=1) / np.sqrt(len(fold_r2)))
             if len(fold_r2) > 1 else 0.0)

    full_prediction = predict_candidate(model, data, features)
    if (not np.isfinite(full_prediction).all()
            or (full_prediction <= 0).any()):
        assessment["reason"] = "完整模型预测含非有限值或非正值"
        return assessment

    assessment.update(_information_criteria(
        y, full_prediction, model.n_params))
    assessment["validation"] = {
        "n_splits": len(folds),
        "cv_r_squared": float(cv_r2),
        "cv_rmse_log": cv_rmse,
        "cv_r2_standard_error": cv_se,
        "fold_r_squared": [float(value) for value in fold_r2],
        "fold_failures": [],
    }
    assessment["success"] = True
    assessment["reason"] = ""
    return assessment


def discover_auto(data: pd.DataFrame, target: str, features: list[str],
                  families: list[str] | None = None,
                  n_splits: int = 5, random_state: int = 42,
                  adjust_backend: bool = True) -> dict:
    """Autonomously validate, rank, recommend, and optionally adjust a model."""
    reason = _validate_input(data, target, features)
    if reason:
        return {"success": False, "reason": reason, "assessments": [],
                "recommended": None, "final_model": None}

    names = families if families is not None else list(FAMILIES)
    assessments = [
        assess_family(data, target, features, name,
                      n_splits=n_splits, random_state=random_state)
        for name in names
    ]
    successful = [item for item in assessments if item["success"]]
    if not successful:
        return {
            "success": False,
            "reason": "没有函数族通过完整拟合与交叉验证",
            "assessments": assessments,
            "recommended": None,
            "final_model": None,
        }

    best_validation = max(
        successful,
        key=lambda item: item["validation"]["cv_r_squared"],
    )
    best_cv = best_validation["validation"]["cv_r_squared"]
    tolerance = max(
        0.005,
        best_validation["validation"]["cv_r2_standard_error"],
    )
    competitive = [
        item for item in successful
        if item["validation"]["cv_r_squared"] >= best_cv - tolerance
    ]
    family_priority = {
        "power_law": 0,
        "piecewise_power_law": 1,
        "rational_log": 2,
    }
    recommended = min(
        competitive,
        key=lambda item: (
            round(item["aicc"], 10),
            round(item["bic"], 10),
            item["model"]["n_params"],
            family_priority.get(item["family"], 100),
            -item["validation"]["cv_r_squared"],
        ),
    )

    successful.sort(
        key=lambda item: (
            item is not recommended,
            -item["validation"]["cv_r_squared"],
            item["aicc"],
        ))
    failed = [item for item in assessments if not item["success"]]
    assessments = successful + failed
    for rank, item in enumerate(assessments, start=1):
        item["rank"] = rank if item["success"] else None
        item["recommended"] = item is recommended

    backend = None
    final_model = recommended["model"]
    if adjust_backend and recommended["family"] == "power_law":
        backend = fit_with_guard(
            data, target, features, adjustment="auto")
        if backend.get("success") and backend.get("after") is not None:
            final_model = backend["after"]
    elif adjust_backend:
        backend = {
            "success": True,
            "applied_adjustment": "not_applicable",
            "reason": (
                f"推荐族 {recommended['family']} 不是全局幂律;"
                "当前后端稳健调整仅作用于 power_law"),
        }

    return {
        "success": True,
        "target": target,
        "features": list(features),
        "assessments": assessments,
        "recommended": recommended,
        "final_model": final_model,
        "backend_adjustment": backend,
        "selection": {
            "rule": (
                "k-fold log-space validation; models within one standard "
                "error (minimum tolerance 0.005) are ordered by AICc/BIC"
            ),
            "best_cv_r_squared": float(best_cv),
            "one_se_tolerance": float(tolerance),
            "n_splits": int(n_splits),
            "random_state": int(random_state),
        },
    }
