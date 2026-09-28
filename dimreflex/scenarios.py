"""Multi-scenario power-law fitting with validated structure selection.

Engineering data often come from several machines, materials, sites, or
operating modes.  This module compares three explicit assumptions:

* pooled: one correlation for every scenario;
* shared_exponents: common mechanism exponents, scenario-specific prefactors;
* independent: a complete correlation per scenario.

Selection uses scenario-stratified k-fold validation and AICc/BIC, so scenarios
are split only when the data justify the extra freedom.
"""


import numpy as np
import pandas as pd

from ._ols import ols as _ols

STRATEGIES = ("pooled", "shared_exponents", "independent")


def _scenario_key(value) -> str:
    return str(value)


def _validate(data: pd.DataFrame, target: str, features: list[str],
              scenario: str) -> str | None:
    columns = [target, scenario, *list(features)]
    missing = [column for column in columns if column not in data.columns]
    if missing:
        return f"数据缺少列 {missing}"
    if not features:
        return "多场景拟合至少需要一个特征"
    if data[scenario].isna().any():
        return f"场景列 '{scenario}' 包含缺失值"
    labels = list(pd.unique(data[scenario]))
    label_keys = [_scenario_key(label) for label in labels]
    if len(set(label_keys)) != len(label_keys):
        return "不同场景标签在字符串化后发生碰撞，请统一场景列类型"
    try:
        values = data[[target, *features]].to_numpy(dtype=float)
    except (TypeError, ValueError):
        return "目标或特征包含非数值数据"
    if not np.isfinite(values).all():
        return "目标或特征包含 NaN/inf"
    if (values <= 0).any():
        return "多场景幂律拟合要求目标和特征均为正"
    counts = data.groupby(scenario, dropna=False, sort=False).size()
    minimum = len(features) + 3
    too_small = {
        _scenario_key(label): int(count)
        for label, count in counts.items() if count < minimum
    }
    if too_small:
        return f"部分场景样本不足(每场景至少 {minimum} 条): {too_small}"
    return None


def _fit_linear(design: np.ndarray, log_y: np.ndarray) -> dict:
    fit = _ols(design, log_y)
    return {
        "coef": fit["coef"],
        "prediction_log": fit["prediction_log"],
        "sse": fit["ss_res"],
        "rank": fit["rank"],
        "condition_number": fit["condition_number"],
    }


def _fit_strategy(data: pd.DataFrame, target: str, features: list[str],
                  scenario: str, strategy: str) -> dict:
    if strategy not in STRATEGIES:
        return {"success": False, "strategy": strategy,
                "reason": f"未知多场景策略 '{strategy}'"}

    log_y = np.log(data[target].to_numpy(dtype=float))
    log_x = np.log(data[features].to_numpy(dtype=float))
    labels = sorted(data[scenario].unique(), key=_scenario_key)
    label_keys = [_scenario_key(label) for label in labels]
    scenario_values = data[scenario].map(_scenario_key).to_numpy()
    n = len(data)

    if strategy == "pooled":
        design = np.column_stack([np.ones(n), log_x])
        fitted = _fit_linear(design, log_y)
        coef = fitted["coef"]
        params = {
            "prefactor": float(np.exp(coef[0])),
            "exponents": {
                feature: float(value)
                for feature, value in zip(features, coef[1:], strict=True)
            },
        }
        n_params = len(features) + 1

    elif strategy == "shared_exponents":
        indicators = np.column_stack(
            [(scenario_values == key).astype(float) for key in label_keys])
        design = np.column_stack([indicators, log_x])
        fitted = _fit_linear(design, log_y)
        coef = fitted["coef"]
        params = {
            "scenario_prefactors": {
                key: float(np.exp(coef[index]))
                for index, key in enumerate(label_keys)
            },
            "exponents": {
                feature: float(value)
                for feature, value in zip(
                        features, coef[len(labels):], strict=True)
            },
        }
        n_params = len(labels) + len(features)

    else:
        scenario_models = {}
        predictions = np.empty(n)
        total_rank = 0
        conditions = []
        for key in label_keys:
            mask = scenario_values == key
            design_local = np.column_stack(
                [np.ones(int(mask.sum())), log_x[mask]])
            local = _fit_linear(design_local, log_y[mask])
            coef = local["coef"]
            predictions[mask] = local["prediction_log"]
            total_rank += local["rank"]
            conditions.append(local["condition_number"])
            scenario_models[key] = {
                "prefactor": float(np.exp(coef[0])),
                "exponents": {
                    feature: float(value)
                    for feature, value in zip(features, coef[1:], strict=True)
                },
            }
        fitted = {
            "prediction_log": predictions,
            "sse": float(np.sum((log_y - predictions) ** 2)),
            "rank": total_rank,
            "condition_number": float(max(conditions)),
        }
        params = {"scenario_models": scenario_models}
        n_params = len(labels) * (len(features) + 1)

    if fitted["rank"] < n_params:
        return {
            "success": False,
            "strategy": strategy,
            "reason": (
                f"设计矩阵秩不足({fitted['rank']} < {n_params})，"
                "无法辨识该多场景结构"),
        }

    sse = max(fitted["sse"], np.finfo(float).tiny)
    ss_tot = float(np.sum((log_y - log_y.mean()) ** 2))
    r2 = 1 - sse / ss_tot if ss_tot > 0 else 0.0
    aic = n * np.log(sse / n) + 2 * n_params
    correction = (
        2 * n_params * (n_params + 1) / (n - n_params - 1)
        if n > n_params + 1 else float("inf")
    )
    return {
        "success": True,
        "strategy": strategy,
        "target": target,
        "features": list(features),
        "scenario_column": scenario,
        "scenario_labels": label_keys,
        "params": params,
        "n_params": int(n_params),
        "n_samples": int(n),
        "r_squared": float(r2),
        "train_rmse_log": float(np.sqrt(sse / n)),
        "aicc": float(aic + correction),
        "bic": float(n * np.log(sse / n) + n_params * np.log(n)),
        "design_rank": int(fitted["rank"]),
        "condition_number": fitted["condition_number"],
    }


def predict_scenarios(model: dict, data: pd.DataFrame) -> np.ndarray:
    """Predict rows using a model returned by fit_scenarios."""
    if not model.get("success"):
        raise ValueError("不能使用失败的多场景模型预测")
    features = model["features"]
    scenario = model["scenario_column"]
    missing = [column for column in [*features, scenario]
               if column not in data.columns]
    if missing:
        raise KeyError(f"预测数据缺少列 {missing}")
    x = data[features].to_numpy(dtype=float)
    if (not np.isfinite(x).all()) or (x <= 0).any():
        raise ValueError("预测特征包含非有限值或非正值")
    scenario_values = data[scenario].map(_scenario_key).to_numpy()
    unknown = sorted(set(scenario_values) - set(model["scenario_labels"]))
    if unknown:
        raise ValueError(f"出现训练时未见场景 {unknown}")

    strategy = model["strategy"]
    params = model["params"]
    if strategy == "pooled":
        prediction = np.full(len(data), params["prefactor"], dtype=float)
        for index, feature in enumerate(features):
            prediction *= x[:, index] ** params["exponents"][feature]
        return prediction

    prediction = np.empty(len(data))
    if strategy == "shared_exponents":
        for key in model["scenario_labels"]:
            mask = scenario_values == key
            values = np.full(
                int(mask.sum()), params["scenario_prefactors"][key],
                dtype=float)
            for index, feature in enumerate(features):
                values *= x[mask, index] ** params["exponents"][feature]
            prediction[mask] = values
        return prediction

    for key in model["scenario_labels"]:
        mask = scenario_values == key
        local = params["scenario_models"][key]
        values = np.full(int(mask.sum()), local["prefactor"], dtype=float)
        for index, feature in enumerate(features):
            values *= x[mask, index] ** local["exponents"][feature]
        prediction[mask] = values
    return prediction


def _stratified_folds(data: pd.DataFrame, scenario: str, n_splits: int,
                      random_state: int) -> list[np.ndarray]:
    groups = list(data.groupby(scenario, sort=False).indices.values())
    splits = min(max(2, int(n_splits)),
                 *(len(indices) for indices in groups))
    rng = np.random.default_rng(random_state)
    parts = [np.array_split(rng.permutation(indices), splits)
             for indices in groups]
    return [
        np.concatenate([scenario_parts[fold] for scenario_parts in parts])
        for fold in range(splits)
    ]


def _assess_strategy(data: pd.DataFrame, target: str, features: list[str],
                     scenario: str, strategy: str, n_splits: int,
                     random_state: int) -> dict:
    model = _fit_strategy(data, target, features, scenario, strategy)
    assessment = {"strategy": strategy, "success": False, "model": model}
    if not model["success"]:
        assessment["reason"] = model["reason"]
        return assessment

    folds = _stratified_folds(data, scenario, n_splits, random_state)
    predictions = np.full(len(data), np.nan)
    all_indices = np.arange(len(data))
    fold_r2 = []
    failures = []
    for fold_index, valid in enumerate(folds):
        train = np.setdiff1d(all_indices, valid, assume_unique=True)
        fold_model = _fit_strategy(
            data.iloc[train], target, features, scenario, strategy)
        try:
            predicted = predict_scenarios(fold_model, data.iloc[valid])
        except Exception as exc:  # noqa: BLE001 — 单折预测失败只记入 failures,不中断整轮评估
            failures.append({
                "fold": fold_index,
                "reason": f"{type(exc).__name__}: {exc}",
            })
            continue
        predictions[valid] = predicted
        actual_log = np.log(
            data.iloc[valid][target].to_numpy(dtype=float))
        predicted_log = np.log(predicted)
        ss_res = float(np.sum((actual_log - predicted_log) ** 2))
        ss_tot = float(np.sum((actual_log - actual_log.mean()) ** 2))
        fold_r2.append(1 - ss_res / ss_tot if ss_tot > 0 else 0.0)

    if failures or not np.isfinite(predictions).all():
        assessment["reason"] = "场景分层交叉验证未完成"
        assessment["validation"] = {"fold_failures": failures}
        return assessment

    actual_log = np.log(data[target].to_numpy(dtype=float))
    predicted_log = np.log(predictions)
    ss_res = float(np.sum((actual_log - predicted_log) ** 2))
    ss_tot = float(np.sum((actual_log - actual_log.mean()) ** 2))
    assessment.update({
        "success": True,
        "reason": "",
        "aicc": model["aicc"],
        "bic": model["bic"],
        "validation": {
            "n_splits": len(folds),
            "cv_r_squared": float(
                1 - ss_res / ss_tot if ss_tot > 0 else 0.0),
            "cv_rmse_log": float(np.sqrt(ss_res / len(data))),
            "cv_r2_standard_error": (
                float(np.std(fold_r2, ddof=1) / np.sqrt(len(fold_r2)))
                if len(fold_r2) > 1 else 0.0
            ),
            "fold_r_squared": [float(value) for value in fold_r2],
            "fold_failures": [],
        },
    })
    return assessment


def fit_scenarios(data: pd.DataFrame, target: str, features: list[str],
                  scenario: str, strategies: list[str] | None = None,
                  n_splits: int = 5, random_state: int = 42) -> dict:
    """Fit and automatically select a multi-scenario correlation structure."""
    reason = _validate(data, target, features, scenario)
    if reason:
        return {"success": False, "reason": reason, "assessments": [],
                "recommended": None, "final_model": None}
    selected = strategies if strategies is not None else list(STRATEGIES)
    assessments = [
        _assess_strategy(
            data, target, features, scenario, strategy,
            n_splits, random_state)
        for strategy in selected
    ]
    successful = [item for item in assessments if item["success"]]
    if not successful:
        return {"success": False, "reason": "没有多场景结构通过验证",
                "assessments": assessments, "recommended": None,
                "final_model": None}

    best = max(successful,
               key=lambda item: item["validation"]["cv_r_squared"])
    best_cv = best["validation"]["cv_r_squared"]
    tolerance = max(
        0.005, best["validation"]["cv_r2_standard_error"])
    competitive = [
        item for item in successful
        if item["validation"]["cv_r_squared"] >= best_cv - tolerance
    ]
    priority = {"pooled": 0, "shared_exponents": 1, "independent": 2}
    recommended = min(
        competitive,
        key=lambda item: (
            item["model"]["n_params"],
            round(item["aicc"], 10),
            round(item["bic"], 10),
            priority.get(item["strategy"], 100),
        ),
    )
    successful.sort(
        key=lambda item: (
            item is not recommended,
            -item["validation"]["cv_r_squared"],
            item["aicc"],
        ))
    failures = [item for item in assessments if not item["success"]]
    assessments = successful + failures
    for rank, assessment in enumerate(assessments, start=1):
        assessment["rank"] = rank if assessment["success"] else None
        assessment["recommended"] = assessment is recommended

    counts = data.groupby(scenario, sort=False).size()
    return {
        "success": True,
        "target": target,
        "features": list(features),
        "scenario_column": scenario,
        "scenario_counts": {
            _scenario_key(label): int(count)
            for label, count in counts.items()
        },
        "assessments": assessments,
        "recommended": recommended,
        "final_model": recommended["model"],
        "selection": {
            "rule": (
                "scenario-stratified k-fold validation; within one standard "
                "error prefer fewer parameters, then lower AICc/BIC"
            ),
            "best_cv_r_squared": float(best_cv),
            "one_se_tolerance": float(tolerance),
            "n_splits": int(n_splits),
            "random_state": int(random_state),
        },
    }
