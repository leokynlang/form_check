"""form_check(幂律形式合适性检查)单元测试。

背景:负对照 validate_friction.py(Prandtl-Karman 对数律=非幂律真值)
实锤"幂律路径对形式失配静默"。这里钉住三件事:
1. 真幂律数据不误报(特异性)
2. 对数律真值报警(灵敏度)
3. discover() 全链路把报警留痕到 warnings(form_check 键存在)
"""

import dimreflex as dr
import numpy as np
import pandas as pd
from dimreflex import form_check


def make_power_law_data(n=60, noise=0.005, seed=7):
    """纯幂律真值 y = 2.5 * x^0.73,log 空间加性噪声(乘性原始空间)。"""
    rng = np.random.default_rng(seed)
    x = np.logspace(1, 5, n)
    y = 2.5 * x ** 0.73 * np.exp(noise * rng.standard_normal(n))
    return pd.DataFrame({"x": x, "y": y})


def make_log_law_data(n=60, noise=0.005, seed=7):
    """非幂律真值:1/sqrt(f) = 2*log10(Re*sqrt(f)) - 0.8 的光滑管律,
    y = 1/f(对数律的化身)。任何有限幂律都覆盖不了它。"""
    from scipy.optimize import brentq
    rng = np.random.default_rng(seed)

    def f_of(re):
        g = lambda f: 1.0 / np.sqrt(f) - 2 * np.log10(re * np.sqrt(f)) + 0.8
        return brentq(g, 1e-6, 0.1)

    x = np.logspace(np.log10(4e3), 7, n)
    f = np.array([f_of(r) for r in x])
    # 噪声加在 log10(f) 上,保持对数律隐式结构不被破坏太多
    log_f = np.log10(f) + noise * rng.standard_normal(n)
    y = 10.0 ** log_f
    return pd.DataFrame({"x": x, "y": y})


class TestFormCheck:
    def test_power_law_data_no_alarm(self):
        """真幂律 + 同量级噪声 → 不得报警(判据 C5 的单元化)。"""
        res = form_check(make_power_law_data(), "y", ["x"])
        assert res["checked"]
        assert not res["alarm"], res["reason"]

    def test_log_law_data_alarms(self):
        """Prandtl-Karman 对数律 → 必须报警(判据 C1 的单元化)。"""
        res = form_check(make_log_law_data(), "y", ["x"])
        assert res["checked"]
        assert res["alarm"], res["reason"]
        assert abs(res["t_stat"]) > 3

    def test_structured_rejection_nonpositive(self):
        """非正值 → 结构化拒绝,不抛异常(同构 schema 纪律)。"""
        df = pd.DataFrame({"x": [1.0, -2.0, 3.0] * 8, "y": [1.0, 2.0, 3.0] * 8})
        res = form_check(df, "y", ["x"])
        assert not res["checked"]
        assert "非正值" in res["reason"]

    def test_too_few_samples_skipped(self):
        df = make_power_law_data(n=8)
        res = form_check(df, "y", ["x"])
        assert not res["checked"]
        assert "过少" in res["reason"]

    def test_deterministic(self):
        """同 seed 两次调用结果逐位一致(框架可复现性纪律)。"""
        a = form_check(make_log_law_data(), "y", ["x"])
        b = form_check(make_log_law_data(), "y", ["x"])
        assert a["t_stat"] == b["t_stat"]
        assert a["curvature"] == b["curvature"]


class TestDiscoverIntegration:
    """discover() 全链路:报警必须留痕到 data_report.warnings。"""

    def _system(self):
        s = dr.System("form_check_test")
        s.add_variable("dp", "Pa", "target")
        s.add_variable("v", "m/s", "state")
        s.add_variable("D", "m", "geometry")
        s.add_variable("rho", "kg/m**3", "material_property")
        s.add_variable("mu", "Pa*s", "material_property")
        return s

    def _friction_frame(self, kind):
        df = make_log_law_data() if kind == "log_law" else make_power_law_data()
        v = df["x"].to_numpy()  # 借用 x 列当流速网格
        rng = np.random.default_rng(11)
        return pd.DataFrame({
            "dp": df["y"].to_numpy() * 100.0,
            "v": v,
            "D": np.full(len(df), 0.1),
            "rho": np.full(len(df), 998.0) * (1 + 0.0 * v + rng.normal(0, 0.001, len(df))),
            "mu": np.full(len(df), 1.0e-3) * (1 + 0.0 * v + rng.normal(0, 0.001, len(df))),
        })

    def test_power_law_path_flags_log_law(self):
        """v0.4 的静默缺口(负对照 C2 预期 FAIL)应在 v0.5 翻转为留痕。"""
        result = dr.discover(self._friction_frame("log_law"),
                             self._system(), target="dp")
        assert result["fit"] is not None and result["fit"]["success"]
        assert result["fit"]["r_squared"] > 0.98  # 拟合本身依旧"漂亮"
        assert result["form_check"]["checked"]
        assert result["form_check"]["alarm"]
        assert any("形式失配" in w
                   for w in result["data_report"]["warnings"])

    def test_power_law_data_no_false_alarm(self):
        """真幂律链路 → 无形式警告(不破坏既有干净输出)。"""
        result = dr.discover(self._friction_frame("power_law"),
                             self._system(), target="dp")
        assert result["fit"] is not None and result["fit"]["success"]
        assert not result["form_check"]["alarm"]
        assert not any("形式失配" in w
                       for w in result["data_report"]["warnings"])
