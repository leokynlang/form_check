"""Validated autonomous discovery and automatic backend adjustment tests."""

from pathlib import Path

import dimreflex as dr
import numpy as np
import pandas as pd

DATA_PATH = Path(__file__).resolve().parent.parent / "heat_exchanger_data.csv"


def power_data(n=90, seed=13, outlier=None):
    rng = np.random.default_rng(seed)
    x = np.exp(rng.uniform(np.log(1.0), np.log(100.0), n))
    y = 2.0 * x ** 0.8 * rng.normal(1.0, 0.02, n)
    if outlier is not None:
        y[outlier] *= 30.0
    return pd.DataFrame({"x": x, "y": y})


def regime_data(seed=17):
    rng = np.random.default_rng(seed)
    x = np.concatenate([
        rng.uniform(500, 2200, 50),
        rng.uniform(2400, 50000, 50),
    ])
    y = np.where(x < 2300, 0.5 * x ** (1 / 3), 0.023 * x ** 0.8)
    y *= rng.normal(1.0, 0.01, len(y))
    return pd.DataFrame({"Re": x, "Nu": y})


def heat_system():
    system = dr.System("heat")
    system.add_variable("rho", "kg/m**3", "material_property")
    system.add_variable("mu", "Pa*s", "material_property")
    system.add_variable("cp", "J/(kg*K)", "material_property")
    system.add_variable("k", "W/(m*K)", "material_property")
    system.add_variable("D", "m", "geometry")
    system.add_variable("v", "m/s", "state")
    system.add_variable("h", "W/(m**2*K)", "target")
    return system


class TestAutonomousSelection:
    def test_pure_power_law_recommends_power_law(self):
        result = dr.discover_auto(
            power_data(), "y", ["x"],
            families=["power_law", "piecewise_power_law", "rational_log"],
        )
        assert result["success"]
        assert result["recommended"]["family"] == "power_law"
        assert result["recommended"]["validation"]["cv_r_squared"] > 0.99
        assert result["selection"]["n_splits"] == 5

    def test_two_regimes_recommend_piecewise(self):
        result = dr.discover_auto(
            regime_data(), "Nu", ["Re"],
            families=["power_law", "piecewise_power_law", "rational_log"],
        )
        assert result["success"]
        assert result["recommended"]["family"] == "piecewise_power_law"
        assert result["recommended"]["model"]["params"]["n_segments"] == 2

    def test_invalid_data_degrades_explicitly(self):
        data = power_data()
        data.loc[3, "x"] = np.nan
        result = dr.discover_auto(data, "y", ["x"])
        assert not result["success"]
        assert result["recommended"] is None
        assert "NaN" in result["reason"]

    def test_assessments_contain_validation_and_complexity(self):
        result = dr.discover_auto(
            power_data(), "y", ["x"], families=["power_law"])
        assessment = result["assessments"][0]
        assert assessment["success"]
        assert np.isfinite(assessment["aicc"])
        assert np.isfinite(assessment["bic"])
        assert assessment["validation"]["fold_failures"] == []


class TestAutomaticBackendAdjustment:
    def test_outlier_automatically_uses_robust_fit(self):
        result = dr.discover_auto(
            power_data(outlier=10), "y", ["x"],
            families=["power_law"])
        backend = result["backend_adjustment"]
        assert backend["success"]
        assert backend["applied_adjustment"] == "robust"
        assert backend["after"]["family"] == "power_law_robust"
        exponent = backend["after"]["params"]["exponents"]["x"]
        assert abs(exponent - 0.8) < 0.05

    def test_collinearity_automatically_uses_ridge(self):
        data = power_data()
        data["x2"] = data["x"] ** 1.00001
        guarded = dr.fit_with_guard(
            data, "y", ["x", "x2"], adjustment="auto")
        assert guarded["success"]
        assert guarded["anomalies"]["collinearity"] is not None
        assert guarded["applied_adjustment"] == "ridge"
        assert guarded["after"]["family"] == "power_law_ridge"
        assert guarded["after"]["params"]["ridge_alpha"] > 0

    def test_clean_data_is_left_unchanged(self):
        guarded = dr.fit_with_guard(
            power_data(), "y", ["x"], adjustment="auto",
            leverage_factor=3.0)
        assert guarded["success"]
        assert guarded["applied_adjustment"] == "report_only"
        assert guarded["after"] is None

    def test_non_numeric_data_is_structured_failure(self):
        data = pd.DataFrame({"x": [1, 2, 3, 4],
                             "y": ["bad", "2", "3", "4"]})
        guarded = dr.fit_with_guard(data, "y", ["x"], adjustment="auto")
        assert not guarded["success"]
        assert any(issue["type"] == "non_numeric"
                   for issue in guarded["health"])


class TestPiSpaceAutonomy:
    def test_auto_pi_discovers_validated_heat_transfer_model(self):
        data = pd.read_csv(DATA_PATH)
        result = dr.discover_auto_pi(
            data, heat_system(), "h",
            families=["power_law", "piecewise_power_law", "rational_log"],
        )
        assert result["autonomous"]["success"]
        assert result["autonomous"]["recommended"]["family"] == "power_law"
        assert result["autonomous"]["recommended"]["validation"][
            "cv_r_squared"] > 0.99
        assert result["data_report"]["level"] == 3
        target_group = next(group for group in result["groups"] if "h" in group)
        assert target_group == {"k": -1.0, "D": 1.0, "h": 1.0}

    def test_discover_method_auto_uses_same_pipeline(self):
        data = pd.read_csv(DATA_PATH)
        result = dr.discover(data, heat_system(), "h", method="auto")
        assert result["autonomous"]["success"]
        assert result["fit"]["success"]
        assert result["assessments"]
