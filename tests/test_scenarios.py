"""Multi-scenario fitting and validated structure selection tests."""

import dimreflex as dr
import numpy as np
import pandas as pd


def scenario_data(mode="shared_exponents", n_per_scenario=90, seed=7):
    rng = np.random.default_rng(seed)
    labels = np.repeat(["A", "B", "C"], n_per_scenario)
    n = len(labels)
    x1 = np.exp(rng.uniform(-1.0, 1.0, n))
    x2 = np.exp(rng.uniform(-1.0, 1.0, n))
    noise = np.exp(rng.normal(0.0, 0.04, n))
    scenario_index = np.repeat(np.arange(3), n_per_scenario)

    if mode == "pooled":
        y = 2.0 * x1 ** 0.8 * x2 ** -0.3
    elif mode == "shared_exponents":
        y = np.array([1.0, 2.0, 4.0])[scenario_index] \
            * x1 ** 0.8 * x2 ** -0.3
    elif mode == "independent":
        exponents = [(0.4, -0.1), (0.8, -0.3), (1.3, 0.2)]
        prefactors = [1.0, 2.0, 4.0]
        y = np.empty(n)
        for index, (b1, b2) in enumerate(exponents):
            mask = scenario_index == index
            y[mask] = prefactors[index] * x1[mask] ** b1 * x2[mask] ** b2
    else:
        raise ValueError(mode)
    return pd.DataFrame(
        {"scenario": labels, "x1": x1, "x2": x2, "y": y * noise})


class TestScenarioSelection:
    def test_no_scenario_effect_prefers_pooled_model(self):
        result = dr.fit_scenarios(
            scenario_data("pooled"), "y", ["x1", "x2"], "scenario")
        assert result["success"]
        assert result["final_model"]["strategy"] == "pooled"

    def test_prefactor_effect_prefers_shared_exponents(self):
        result = dr.fit_scenarios(
            scenario_data("shared_exponents"), "y", ["x1", "x2"],
            "scenario")
        assert result["final_model"]["strategy"] == "shared_exponents"
        exponents = result["final_model"]["params"]["exponents"]
        assert abs(exponents["x1"] - 0.8) < 0.03
        assert abs(exponents["x2"] + 0.3) < 0.03

    def test_mechanism_effect_prefers_independent_models(self):
        result = dr.fit_scenarios(
            scenario_data("independent"), "y", ["x1", "x2"],
            "scenario")
        assert result["final_model"]["strategy"] == "independent"
        assert result["recommended"]["validation"]["cv_r_squared"] > 0.99

    def test_every_assessment_has_stratified_validation(self):
        result = dr.fit_scenarios(
            scenario_data(), "y", ["x1", "x2"], "scenario")
        assert all(item["success"] for item in result["assessments"])
        assert all(item["validation"]["fold_failures"] == []
                   for item in result["assessments"])
        assert "fewer parameters" in result["selection"]["rule"]


class TestScenarioPredictionAndGuards:
    def test_prediction_uses_correct_scenario_prefactor(self):
        data = scenario_data("shared_exponents")
        model = dr.fit_scenarios(
            data, "y", ["x1", "x2"], "scenario")["final_model"]
        probe = pd.DataFrame({
            "scenario": ["A", "B", "C"],
            "x1": [2.0, 2.0, 2.0],
            "x2": [1.0, 1.0, 1.0],
        })
        prediction = dr.predict_scenarios(model, probe)
        ratios = prediction / prediction[0]
        assert np.allclose(ratios, [1.0, 2.0, 4.0], rtol=0.04)

    def test_unseen_scenario_is_rejected(self):
        data = scenario_data()
        model = dr.fit_scenarios(
            data, "y", ["x1", "x2"], "scenario")["final_model"]
        probe = pd.DataFrame(
            {"scenario": ["new"], "x1": [2.0], "x2": [1.0]})
        try:
            dr.predict_scenarios(model, probe)
        except ValueError as exc:
            assert "未见场景" in str(exc)
        else:
            raise AssertionError("unseen scenario should fail")

    def test_small_scenario_degrades_explicitly(self):
        data = scenario_data(n_per_scenario=5)
        data = pd.concat([data, data[data["scenario"] != "C"]],
                         ignore_index=True)
        data = data.drop(data[data["scenario"] == "C"].index[3:])
        result = dr.fit_scenarios(
            data, "y", ["x1", "x2"], "scenario")
        assert not result["success"]
        assert "样本不足" in result["reason"]

    def test_nonpositive_data_degrades_explicitly(self):
        data = scenario_data()
        data.loc[0, "x1"] = 0.0
        result = dr.fit_scenarios(
            data, "y", ["x1", "x2"], "scenario")
        assert not result["success"]
        assert "均为正" in result["reason"]

    def test_stringified_scenario_collision_is_rejected(self):
        data = scenario_data("pooled")
        data["scenario"] = pd.Series(
            [1 if index % 2 == 0 else "1" for index in range(len(data))],
            dtype=object,
        )
        result = dr.fit_scenarios(
            data, "y", ["x1", "x2"], "scenario")
        assert not result["success"]
        assert "碰撞" in result["reason"]
