"""form_check_battery(形式合格性电池)补丁验证测试。

钉住五件事:
1. 真幂律:电池 checked 且无报警,可行族含 power
2. 对数律(非幂律非电池族):报警=无可行族升级路径
3. rational 真值:wild 零分布下 viable 含 rational(wild 修复的回归钉)
4. discover() 全链路:form_check_battery 键存在、passthrough 三入口一致
5. 确定性:同输入两次调用 viable 逐位一致
"""
import importlib
import sys

import numpy as np
import pandas as pd

import dimreflex as dr
from dimreflex.discover import _battery_check, FAM_NAMES


def make_power_law_data(n=60, noise=0.005, seed=7):
    rng = np.random.default_rng(seed)
    x = np.logspace(1, 5, n)
    y = 2.5 * x ** 0.73 * np.exp(noise * rng.standard_normal(n))
    return pd.DataFrame({"x": x, "y": y})


def make_rational_data(n=60, noise=0.005, seed=7):
    """饱和曲线真值 y = x/(0.2+0.1x),乘性噪声(rational 族的最好情形)。"""
    rng = np.random.default_rng(seed)
    x = np.logspace(0, 2, n)
    y = x / (0.2 + 0.1 * x) * np.exp(noise * rng.standard_normal(n))
    return pd.DataFrame({"x": x, "y": y})


def make_log_law_data(n=60, noise=0.005, seed=7):
    """对数型:y = 0.8 + 1.1*ln(x),域平移保证正值。不属于电池任何族。"""
    rng = np.random.default_rng(seed)
    x = np.linspace(1.5, 100, n)
    y = (0.8 + 1.1 * np.log(x)) * np.exp(noise * rng.standard_normal(n))
    return pd.DataFrame({"x": x, "y": y})


class TestBatteryUnit:
    def test_power_truth_power_viable_no_alarm(self):
        """真幂律:power 族留在 viable。多种子扫描断言(≥8/10)——
        单一固定种子会抽进 wild 尺寸尾巴(n=60 实测 6.5%),那是
        校准行为不是缺陷;测试只钉机制,不钉彩票。"""
        kept = 0
        for s in range(10):
            res = _battery_check(make_power_law_data(seed=s), "y", ["x"],
                                 power_payload={"checked": True,
                                                "alarm": False})
            assert res["checked"]
            kept += "power" in res["viable"]
        assert kept >= 8, f"power viable 率 {kept}/10,低于校准预期"

    def test_logarithmic_truth_escalates(self):
        """对数型真值不在电池族中:报警 + 无可行族(升级语义)。"""
        fc = {"checked": True, "alarm": True}
        res = _battery_check(make_log_law_data(), "y", ["x"],
                             power_payload=fc)
        assert res["checked"]
        assert res["alarm"], res["reason"]
        assert res["viable"] == []
        assert "无可行族" in res["reason"]

    def test_rational_truth_viable_wild(self):
        """rational 真值必须留在 viable(wild 修复回归钉;排列版会假拒绝)。"""
        res = _battery_check(make_rational_data(), "y", ["x"],
                             power_payload={"checked": True, "alarm": True})
        assert res["checked"]
        assert "rational" in res["viable"], res["reason"]

    def test_power_self_endorsement_blocked(self):
        """幂律被 form_check 拒绝时,power 不得为自己背书。"""
        df = make_power_law_data()
        res = _battery_check(df, "y", ["x"],
                             power_payload={"checked": True, "alarm": True})
        assert "power" not in res["viable"]

    def test_unchecked_power_payload_skips(self):
        res = _battery_check(make_power_law_data(), "y", ["x"],
                             power_payload={"checked": False, "reason": "x"})
        assert not res["checked"]
        assert "跳过" in res["reason"]

    def test_nonpositive_structured_reject(self):
        df = pd.DataFrame({"x": [1.0, -2.0, 3.0] * 8,
                           "y": [1.0, 2.0, 3.0] * 8})
        res = _battery_check(df, "y", ["x"],
                             power_payload={"checked": True, "alarm": False})
        assert not res["checked"]
        assert "非正值" in res["reason"]

    def test_deterministic(self):
        a = _battery_check(make_rational_data(), "y", ["x"],
                           power_payload={"checked": True, "alarm": True})
        b = _battery_check(make_rational_data(), "y", ["x"],
                           power_payload={"checked": True, "alarm": True})
        assert a["viable"] == b["viable"]
        assert a["p_rational"] == b["p_rational"]


class TestBatteryIntegration:
    def _system(self):
        s = dr.System("battery_test")
        s.add_variable("dp", "Pa", "target")
        s.add_variable("v", "m/s", "state")
        s.add_variable("D", "m", "geometry")
        s.add_variable("rho", "kg/m**3", "material_property")
        s.add_variable("mu", "Pa*s", "material_property")
        return s

    def _frame(self, df):
        v = df["x"].to_numpy()
        rng = np.random.default_rng(11)
        return pd.DataFrame({
            "dp": df["y"].to_numpy() * 100.0,
            "v": v,
            "D": np.full(len(df), 0.1),
            "rho": np.full(len(df), 998.0) * (1 + rng.normal(0, 0.001, len(df))),
            "mu": np.full(len(df), 1.0e-3) * (1 + rng.normal(0, 0.001, len(df))),
        })

    def test_discover_passthrough_and_no_noise(self):
        """真幂律全链路:键存在、无电池级警告、可行族含 power。
        3 个数据种子扫描(≥2/3 含 power),全部 0 电池级警告。"""
        ok = 0
        for s in range(3):
            df = make_power_law_data(seed=100 + s)
            result = dr.discover(self._frame(df), self._system(),
                                 target="dp")
            assert result["fit"] is not None and result["fit"]["success"]
            bat = result["form_check_battery"]
            assert bat is not None and bat["checked"]
            assert not bat["alarm"]
            assert not any("形式合格性电池" in w
                           for w in result["data_report"]["warnings"])
            ok += "power" in bat["viable"]
        assert ok >= 2, f"power viable {ok}/3,低于校准预期"

    def test_auto_and_multi_passthrough(self):
        df = make_rational_data()
        r1 = dr.discover_auto_pi(self._frame(df), self._system(), target="dp")
        assert r1["form_check_battery"] is not None
        r2 = dr.discover_multi_pi(self._frame(df), self._system(), target="dp")
        assert r2["form_check_battery"] is not None
        assert r1["form_check_battery"]["viable"] == \
               r2["form_check_battery"]["viable"]
