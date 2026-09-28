"""Engineering batch monitoring, drift, and persistent alert tests."""

import dimreflex as dr
import numpy as np
import pandas as pd


def power_batch(n=400, seed=11, factor=1.0, x_factor=1.0):
    rng = np.random.default_rng(seed)
    x = np.exp(rng.uniform(np.log(1.0), np.log(100.0), n)) * x_factor
    y = factor * 2.0 * x ** 0.8 * np.exp(rng.normal(0.0, 0.025, n))
    return pd.DataFrame({"x": x, "y": y})


def fitted_monitor():
    baseline = power_batch()
    model = dr.PowerLawFamily().fit(baseline, "y", ["x"])
    return dr.EngineeringMonitor(baseline, "y", ["x"], model)


class TestBatchMonitoring:
    def test_clean_batch_stays_normal(self):
        report = fitted_monitor().inspect(power_batch(seed=19), batch_id="ok")
        assert report["success"]
        assert report["status"] == "normal"
        assert report["alerts"] == []

    def test_operating_envelope_and_feature_drift_are_detected(self):
        report = fitted_monitor().inspect(
            power_batch(seed=19, x_factor=100.0), batch_id="ood")
        codes = {alert["code"] for alert in report["alerts"]}
        assert report["status"] == "critical"
        assert "operating_envelope_excursion" in codes
        assert "feature_distribution_drift" in codes

    def test_relationship_drift_is_detected_without_feature_drift(self):
        report = fitted_monitor().inspect(
            power_batch(seed=19, factor=1.8), batch_id="degraded")
        codes = {alert["code"] for alert in report["alerts"]}
        assert "residual_center_drift" in codes
        assert "model_error_degradation" in codes
        assert report["metrics"]["features"]["x"]["psi"] < 0.20

    def test_row_anomaly_positions_are_reported(self):
        batch = power_batch(seed=19)
        batch.loc[7, "y"] *= 20.0
        report = fitted_monitor().inspect(batch)
        assert 7 in report["row_anomaly_positions"]
        assert any(alert["code"] == "residual_outliers"
                   for alert in report["alerts"])

    def test_feature_only_monitoring_does_not_require_target(self):
        report = fitted_monitor().inspect(
            power_batch(seed=19).drop(columns="y"),
            require_target=False)
        assert report["success"]
        assert report["metrics"]["target_available"] is False
        assert "residuals" not in report["metrics"]

    def test_invalid_batch_is_structured_critical_failure(self):
        batch = power_batch(seed=19)
        batch.loc[0, "x"] = 0.0
        report = fitted_monitor().inspect(batch, batch_id="bad")
        assert not report["success"]
        assert report["status"] == "critical"
        assert report["alerts"][0]["code"] == "data_validation_failure"

    def test_non_dataframe_is_structured_failure(self):
        report = fitted_monitor().inspect(None)
        assert not report["success"]
        assert "DataFrame" in report["alerts"][0]["message"]


class TestStatefulEscalation:
    def test_three_abnormal_batches_escalate_persistence(self):
        monitor = fitted_monitor()
        for index in range(3):
            report = monitor.inspect(
                power_batch(seed=20 + index, factor=1.8),
                batch_id=f"shift-{index}")
        assert report["consecutive_abnormal_batches"] == 3
        assert any(alert["code"] == "persistent_abnormal_batches"
                   for alert in report["alerts"])
        assert len(monitor.snapshot()["history"]) == 3

    def test_normal_batch_resets_streak_and_reset_clears_history(self):
        monitor = fitted_monitor()
        monitor.inspect(power_batch(seed=21, factor=1.8))
        normal = monitor.inspect(power_batch(seed=22))
        assert normal["consecutive_abnormal_batches"] == 0
        monitor.reset_history()
        assert monitor.snapshot()["history"] == []


class TestAutonomousMonitorAndScenarios:
    def test_fit_baseline_discovers_and_builds_monitor(self):
        built = dr.EngineeringMonitor.fit_baseline(
            power_batch(), "y", ["x"])
        assert built["success"]
        assert built["discovery"]["success"]
        assert isinstance(built["monitor"], dr.EngineeringMonitor)

    def test_scenario_monitor_rejects_unknown_scenario(self):
        rng = np.random.default_rng(5)
        scenario = np.repeat(["A", "B"], 100)
        x = np.exp(rng.uniform(0.0, 3.0, len(scenario)))
        prefactor = np.where(scenario == "A", 1.0, 3.0)
        y = prefactor * x ** 0.6 * np.exp(rng.normal(0.0, 0.02, len(x)))
        baseline = pd.DataFrame({"scenario": scenario, "x": x, "y": y})
        built = dr.EngineeringMonitor.fit_baseline(
            baseline, "y", ["x"], scenario_column="scenario",
            monitor_scenario_mix=False)
        assert built["success"]
        unknown = pd.DataFrame({"scenario": ["C"], "x": [2.0], "y": [2.0]})
        report = built["monitor"].inspect(unknown)
        assert not report["success"]
        assert "未见场景" in report["alerts"][0]["message"]
