"""Stateful monitoring for engineering correlations in batch operation.

The monitor combines four complementary signals instead of treating every
unusual row as a bad measurement:

* schema and numerical validity (hard data failures),
* operating-envelope excursions (extrapolation risk),
* feature and scenario drift (the process population changed),
* residual drift (the fitted engineering relationship changed).

It also keeps a consecutive-abnormal-batch counter so a weak but persistent
signal can be escalated without overreacting to one noisy batch.
"""

from copy import deepcopy
from threading import RLock

import numpy as np
import pandas as pd

from .autodiscovery import discover_auto, predict_candidate
from .families import CandidateModel
from .scenarios import fit_scenarios, predict_scenarios

DEFAULT_THRESHOLDS = {
    "feature_psi_warning": 0.20,
    "feature_psi_critical": 0.35,
    "outside_rate_warning": 0.10,
    "outside_rate_critical": 0.25,
    "residual_outlier_z": 3.5,
    "residual_outlier_rate_warning": 0.05,
    "residual_outlier_rate_critical": 0.15,
    "residual_shift_z_warning": 2.5,
    "residual_shift_z_critical": 4.0,
    "rmse_ratio_warning": 1.5,
    "rmse_ratio_critical": 2.5,
    "scenario_mix_warning": 0.20,
    "scenario_mix_critical": 0.35,
    "persistent_batches": 3,
}


def _robust_scale(values: np.ndarray) -> float:
    median = float(np.median(values))
    mad_scale = 1.4826 * float(np.median(np.abs(values - median)))
    standard = float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
    return max(mad_scale, standard, 1e-8)


def _psi_reference(values: np.ndarray) -> dict:
    quantiles = np.quantile(values, [0.2, 0.4, 0.6, 0.8])
    internal = np.unique(quantiles)
    edges = np.concatenate(([-np.inf], internal, [np.inf]))
    counts, _ = np.histogram(values, bins=edges)
    probabilities = counts.astype(float) / max(1, counts.sum())
    return {"edges": edges, "probabilities": probabilities}


def population_stability_index(reference: dict,
                               current: np.ndarray) -> float:
    """Return PSI for a current numeric population and a stored reference."""
    values = np.asarray(current, dtype=float)
    counts, _ = np.histogram(values, bins=reference["edges"])
    current_probabilities = counts.astype(float) / max(1, counts.sum())
    expected = np.asarray(reference["probabilities"], dtype=float)
    epsilon = 1e-6
    actual = np.clip(current_probabilities, epsilon, None)
    expected = np.clip(expected, epsilon, None)
    return float(np.sum((actual - expected) * np.log(actual / expected)))


def _as_candidate(model: CandidateModel | dict) -> CandidateModel:
    if isinstance(model, CandidateModel):
        return model
    if isinstance(model, dict) and "family" in model:
        return CandidateModel.from_dict(model)
    raise TypeError("普通监视模型必须是 CandidateModel 或其字典形式")


def _severity(value: float, warning: float, critical: float) -> str | None:
    if value >= critical:
        return "critical"
    if value >= warning:
        return "warning"
    return None


class EngineeringMonitor:
    """Monitor new engineering batches against a fitted baseline.

    Parameters are deliberately explicit: the fitted model remains immutable,
    while only the alert streak and bounded history change between batches.
    """

    def __init__(self, baseline: pd.DataFrame, target: str,
                 features: list[str], model: CandidateModel | dict,
                 scenario_column: str | None = None,
                 thresholds: dict | None = None,
                 history_limit: int = 100,
                 monitor_scenario_mix: bool = True):
        self.target = target
        self.features = list(features)
        self.scenario_column = scenario_column
        self.model = deepcopy(model)
        self._is_scenario_model = (
            isinstance(model, dict)
            and model.get("strategy") in {
                "pooled", "shared_exponents", "independent"})
        if self._is_scenario_model:
            if not model.get("success"):
                raise ValueError("不能使用失败的多场景模型建立监视器")
            model_scenario = model.get("scenario_column")
            if scenario_column is None:
                self.scenario_column = model_scenario
            elif model_scenario != scenario_column:
                raise ValueError("模型与监视器的场景列不一致")
            if list(model.get("features", [])) != self.features:
                raise ValueError("模型与监视器的特征列表不一致")
            if model.get("target") != target:
                raise ValueError("模型与监视器的目标列不一致")
        else:
            candidate = _as_candidate(model)
            if not candidate.success:
                raise ValueError("不能使用失败的候选模型建立监视器")
            if scenario_column is not None:
                raise ValueError("提供了场景列，但模型不是多场景模型")
        self.thresholds = dict(DEFAULT_THRESHOLDS)
        if thresholds:
            unknown = set(thresholds) - set(self.thresholds)
            if unknown:
                raise ValueError(f"未知监视阈值 {sorted(unknown)}")
            self.thresholds.update(thresholds)
        self.history_limit = max(1, int(history_limit))
        self.monitor_scenario_mix = bool(monitor_scenario_mix)
        self.history: list[dict] = []
        self.consecutive_abnormal_batches = 0
        self._lock = RLock()

        reason = self._validate_frame(baseline, require_target=True)
        if reason:
            raise ValueError(f"基线数据无效: {reason}")

        prediction = self._predict(baseline)
        if (not np.isfinite(prediction).all()) or (prediction <= 0).any():
            raise ValueError("基线模型预测包含非有限值或非正值")
        residual = np.log(baseline[target].to_numpy(dtype=float)) - np.log(
            prediction)
        self._residual_center = float(np.median(residual))
        self._residual_scale = _robust_scale(residual)
        self._baseline_rmse_log = max(
            float(np.sqrt(np.mean(residual ** 2))), 1e-8)

        self._feature_reference = {}
        for feature in self.features:
            log_values = np.log(baseline[feature].to_numpy(dtype=float))
            q01, q25, q75, q99 = np.quantile(
                log_values, [0.01, 0.25, 0.75, 0.99])
            spread = max(float(q75 - q25), 1e-12)
            self._feature_reference[feature] = {
                "lower_log": float(q01 - 0.5 * spread),
                "upper_log": float(q99 + 0.5 * spread),
                "psi": _psi_reference(log_values),
            }

        self._scenario_probabilities = None
        if self.scenario_column is not None:
            keys = baseline[self.scenario_column].map(str)
            counts = keys.value_counts(normalize=True)
            self._scenario_probabilities = {
                key: float(counts.get(key, 0.0))
                for key in self.model["scenario_labels"]
            }

    @classmethod
    def fit_baseline(cls, baseline: pd.DataFrame, target: str,
                     features: list[str],
                     scenario_column: str | None = None,
                     n_splits: int = 5, random_state: int = 42,
                     thresholds: dict | None = None,
                     history_limit: int = 100,
                     monitor_scenario_mix: bool = True) -> dict:
        """Discover a model and build a monitor in one guarded operation."""
        if scenario_column is None:
            discovery = discover_auto(
                baseline, target, features, n_splits=n_splits,
                random_state=random_state)
        else:
            discovery = fit_scenarios(
                baseline, target, features, scenario_column,
                n_splits=n_splits, random_state=random_state)
        if not discovery.get("success"):
            return {
                "success": False,
                "reason": discovery.get("reason", "基线模型发现失败"),
                "discovery": discovery,
                "monitor": None,
            }
        try:
            monitor = cls(
                baseline, target, features, discovery["final_model"],
                scenario_column=scenario_column, thresholds=thresholds,
                history_limit=history_limit,
                monitor_scenario_mix=monitor_scenario_mix)
        except Exception as exc:  # noqa: BLE001 — 监视器构造失败必须折叠为结构化失败,不抛给调用方
            return {
                "success": False,
                "reason": f"{type(exc).__name__}: {exc}",
                "discovery": discovery,
                "monitor": None,
            }
        return {"success": True, "reason": "", "discovery": discovery,
                "monitor": monitor}

    def _validate_frame(self, data: pd.DataFrame,
                        require_target: bool) -> str | None:
        if not isinstance(data, pd.DataFrame):
            return "输入必须是 pandas.DataFrame"
        required = list(self.features)
        if require_target:
            required.append(self.target)
        if self.scenario_column is not None:
            required.append(self.scenario_column)
        missing = [column for column in required if column not in data.columns]
        if missing:
            return f"缺少列 {missing}"
        if len(data) == 0:
            return "批次为空"
        numeric = list(self.features)
        if require_target:
            numeric.append(self.target)
        try:
            values = data[numeric].to_numpy(dtype=float)
        except (TypeError, ValueError):
            return "目标或特征包含非数值数据"
        if not np.isfinite(values).all():
            return "目标或特征包含 NaN/inf"
        if (values <= 0).any():
            return "幂律监视要求目标和特征均为正"
        if (self.scenario_column is not None
                and data[self.scenario_column].isna().any()):
            return f"场景列 '{self.scenario_column}' 包含缺失值"
        return None

    def _predict(self, data: pd.DataFrame) -> np.ndarray:
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            if self._is_scenario_model:
                return predict_scenarios(self.model, data)
            candidate = _as_candidate(self.model)
            return predict_candidate(candidate, data, self.features)

    @staticmethod
    def _alert(code: str, severity: str, message: str,
               metrics: dict | None = None) -> dict:
        return {"code": code, "severity": severity, "message": message,
                "metrics": metrics or {}}

    def _record(self, report: dict) -> None:
        compact = {
            "batch_id": report["batch_id"],
            "status": report["status"],
            "n_rows": report["metrics"].get("n_rows", 0),
            "alert_codes": [alert["code"] for alert in report["alerts"]],
            "consecutive_abnormal_batches": (
                report["consecutive_abnormal_batches"]),
        }
        self.history.append(compact)
        if len(self.history) > self.history_limit:
            del self.history[:-self.history_limit]

    def _failure_report(self, batch_id, reason: str,
                        n_rows: int = 0) -> dict:
        self.consecutive_abnormal_batches += 1
        report = {
            "success": False,
            "batch_id": batch_id,
            "status": "critical",
            "alerts": [self._alert(
                "data_validation_failure", "critical", reason)],
            "metrics": {"n_rows": int(n_rows)},
            "row_anomaly_positions": [],
            "consecutive_abnormal_batches": (
                self.consecutive_abnormal_batches),
        }
        self._record(report)
        return report

    def inspect(self, batch: pd.DataFrame, batch_id=None,
                require_target: bool = True) -> dict:
        """Inspect one batch and update the persistent alert streak/history."""
        with self._lock:
            reason = self._validate_frame(batch, require_target=require_target)
            if reason:
                n_rows = len(batch) if hasattr(batch, "__len__") else 0
                return self._failure_report(batch_id, reason, n_rows)

            alerts = []
            metrics: dict = {
                "n_rows": len(batch),
                "target_available": bool(require_target),
                "features": {},
            }

            if self.scenario_column is not None:
                scenario_values = batch[self.scenario_column].map(str)
                known = set(self.model["scenario_labels"])
                unknown = sorted(set(scenario_values) - known)
                if unknown:
                    return self._failure_report(
                        batch_id, f"出现训练时未见场景 {unknown}", len(batch))
                proportions = scenario_values.value_counts(normalize=True)
                current = {
                    key: float(proportions.get(key, 0.0)) for key in known}
                total_variation = 0.5 * sum(
                    abs(current[key] - self._scenario_probabilities[key])
                    for key in known)
                metrics["scenario_mix"] = {
                    "baseline": dict(self._scenario_probabilities),
                    "current": current,
                    "total_variation": float(total_variation),
                }
                if self.monitor_scenario_mix:
                    level = _severity(
                        total_variation,
                        self.thresholds["scenario_mix_warning"],
                        self.thresholds["scenario_mix_critical"])
                    if level:
                        alerts.append(self._alert(
                            "scenario_mix_drift", level,
                            "场景构成相对基线发生漂移",
                            {"total_variation": float(total_variation)}))

            for feature in self.features:
                values_log = np.log(batch[feature].to_numpy(dtype=float))
                reference = self._feature_reference[feature]
                outside = ((values_log < reference["lower_log"])
                           | (values_log > reference["upper_log"]))
                outside_rate = float(np.mean(outside))
                psi = population_stability_index(reference["psi"], values_log)
                metrics["features"][feature] = {
                    "outside_rate": outside_rate,
                    "psi": psi,
                    "observed_min": float(batch[feature].min()),
                    "observed_max": float(batch[feature].max()),
                    "envelope_min": float(np.exp(reference["lower_log"])),
                    "envelope_max": float(np.exp(reference["upper_log"])),
                }
                outside_level = _severity(
                    outside_rate,
                    self.thresholds["outside_rate_warning"],
                    self.thresholds["outside_rate_critical"])
                if outside_level:
                    alerts.append(self._alert(
                        "operating_envelope_excursion", outside_level,
                        f"特征 '{feature}' 超出基线运行包络",
                        {"feature": feature, "outside_rate": outside_rate}))
                psi_level = _severity(
                    psi, self.thresholds["feature_psi_warning"],
                    self.thresholds["feature_psi_critical"])
                if psi_level:
                    alerts.append(self._alert(
                        "feature_distribution_drift", psi_level,
                        f"特征 '{feature}' 的分布相对基线漂移",
                        {"feature": feature, "psi": psi}))

            row_anomalies: list[int] = []
            if require_target:
                try:
                    prediction = self._predict(batch)
                except Exception as exc:  # noqa: BLE001 — 批次预测失败转为失败报告,监视循环不得中断
                    return self._failure_report(
                        batch_id, f"模型预测失败: {type(exc).__name__}: {exc}",
                        len(batch))
                if ((not np.isfinite(prediction).all())
                        or (prediction <= 0).any()):
                    return self._failure_report(
                        batch_id, "模型预测包含非有限值或非正值", len(batch))
                residual = np.log(batch[self.target].to_numpy(dtype=float)) \
                    - np.log(prediction)
                z = (residual - self._residual_center) / self._residual_scale
                mask = np.abs(z) >= self.thresholds["residual_outlier_z"]
                row_anomalies = np.flatnonzero(mask).astype(int).tolist()
                outlier_rate = float(np.mean(mask))
                residual_center = float(np.median(residual))
                shift_z = abs(
                    residual_center - self._residual_center) \
                    / self._residual_scale
                rmse = float(np.sqrt(np.mean(residual ** 2)))
                rmse_ratio = rmse / self._baseline_rmse_log
                metrics["residuals"] = {
                    "median_log_residual": residual_center,
                    "median_shift_z": float(shift_z),
                    "rmse_log": rmse,
                    "baseline_rmse_log": self._baseline_rmse_log,
                    "rmse_ratio": float(rmse_ratio),
                    "outlier_rate": outlier_rate,
                    "max_abs_robust_z": float(np.max(np.abs(z))),
                }
                outlier_level = _severity(
                    outlier_rate,
                    self.thresholds["residual_outlier_rate_warning"],
                    self.thresholds["residual_outlier_rate_critical"])
                if row_anomalies:
                    alerts.append(self._alert(
                        "residual_outliers", outlier_level or "warning",
                        "批次中检测到模型残差异常点",
                        {"outlier_rate": outlier_rate,
                         "positions": row_anomalies}))
                shift_level = _severity(
                    shift_z,
                    self.thresholds["residual_shift_z_warning"],
                    self.thresholds["residual_shift_z_critical"])
                if shift_level:
                    alerts.append(self._alert(
                        "residual_center_drift", shift_level,
                        "模型残差中心发生系统性偏移",
                        {"median_shift_z": float(shift_z)}))
                rmse_level = _severity(
                    rmse_ratio,
                    self.thresholds["rmse_ratio_warning"],
                    self.thresholds["rmse_ratio_critical"])
                if rmse_level:
                    alerts.append(self._alert(
                        "model_error_degradation", rmse_level,
                        "模型批次误差相对基线恶化",
                        {"rmse_ratio": float(rmse_ratio)}))

            if alerts:
                self.consecutive_abnormal_batches += 1
            else:
                self.consecutive_abnormal_batches = 0
            persistent = int(self.thresholds["persistent_batches"])
            if alerts and self.consecutive_abnormal_batches >= persistent:
                alerts.append(self._alert(
                    "persistent_abnormal_batches", "critical",
                    "工程异常已连续出现，建议升级处置而非继续观察",
                    {"consecutive_batches":
                     self.consecutive_abnormal_batches}))

            status = "normal"
            if any(alert["severity"] == "critical" for alert in alerts):
                status = "critical"
            elif alerts:
                status = "warning"
            report = {
                "success": True,
                "batch_id": batch_id,
                "status": status,
                "alerts": alerts,
                "metrics": metrics,
                "row_anomaly_positions": row_anomalies,
                "consecutive_abnormal_batches": (
                    self.consecutive_abnormal_batches),
            }
            self._record(report)
            return report

    def snapshot(self) -> dict:
        """Return a JSON-friendly summary of monitor state and baseline."""
        with self._lock:
            return {
                "target": self.target,
                "features": list(self.features),
                "scenario_column": self.scenario_column,
                "baseline_rmse_log": self._baseline_rmse_log,
                "residual_center": self._residual_center,
                "residual_scale": self._residual_scale,
                "thresholds": dict(self.thresholds),
                "consecutive_abnormal_batches": (
                    self.consecutive_abnormal_batches),
                "history": deepcopy(self.history),
            }

    def reset_history(self) -> None:
        """Clear operational state without changing the calibrated baseline."""
        with self._lock:
            self.history.clear()
            self.consecutive_abnormal_batches = 0
