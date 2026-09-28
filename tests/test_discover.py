"""discover.py / templates.py 集成测试:用合成数据验证完整链路。"""

from pathlib import Path

import dimreflex as dr
import numpy as np
import pandas as pd
import pytest

DATA_PATH = Path(__file__).resolve().parent.parent / "heat_exchanger_data.csv"


def make_system():
    sys = dr.System("tube_side_heat_transfer")
    sys.add_variable("rho", "kg/m**3", "material_property")
    sys.add_variable("mu", "Pa*s", "material_property")
    sys.add_variable("cp", "J/(kg*K)", "material_property")
    sys.add_variable("k", "W/(m*K)", "material_property")
    sys.add_variable("D", "m", "geometry")
    sys.add_variable("v", "m/s", "state")
    sys.add_variable("h", "W/(m**2*K)", "target")
    return sys


@pytest.fixture(scope="module")
def data():
    return pd.read_csv(DATA_PATH)


def test_template_recovers_dittus_boelter(data):
    """模板路径应恢复 Nu = 0.023 * Re^0.8 * Pr^0.4(合成数据由该公式生成)。"""
    sys = make_system()
    Re = dr.Group({"rho": 1, "v": 1, "D": 1, "mu": -1}, "Re")
    Pr = dr.Group({"cp": 1, "mu": 1, "k": -1}, "Pr")
    Nu = dr.Group({"h": 1, "D": 1, "k": -1}, "Nu")
    result = dr.fit_with_templates(data, sys, target_group=Nu, feature_groups=[Re, Pr])
    assert result["success"]
    assert result["r_squared"] > 0.99
    assert result["prefactor"] == pytest.approx(0.023, rel=0.1)
    assert result["exponents"]["Re"] == pytest.approx(0.8, abs=0.05)
    assert result["exponents"]["Pr"] == pytest.approx(0.4, abs=0.05)


def test_discover_succeeds_with_target_isolation(data):
    """自动路径应成功且拟合质量高;含 target 的群恰好一个。"""
    sys = make_system()
    result = dr.discover(data, sys, target="h")
    assert result["fit"] is not None
    assert result["fit"]["success"]
    assert result["fit"]["r_squared"] > 0.99

    groups = result["groups"]
    target_groups = [g for g in groups if "h" in g]
    assert len(target_groups) == 1


def test_fit_power_law_rejects_nonpositive():
    """含非正值的输入应在 log 空间拟合前被拒绝。"""
    df = pd.DataFrame({"a": [1.0, 2.0, 3.0], "b": [1.0, 2.0, 0.0]})
    result = dr.fit_power_law(df, target="a", features=["b"])
    assert not result["success"]


def test_check_data_levels(data):
    """数据就绪分级:完整数据(80 行 ≥ 10×7)应达到 Level 2。"""
    sys = make_system()
    report = sys.check_data(data)
    assert report["level"] >= 2
    assert report["missing_vars"] == []


def test_check_data_missing_variable(data):
    """缺列时应报告缺失且等级为 0,绝不硬输出。"""
    sys = make_system()
    bad = data.drop(columns=["cp"])
    report = sys.check_data(bad)
    assert report["level"] == 0
    assert "cp" in report["missing_vars"]


# ---------------------------------------------------------------------------
# NaN / 非有限值防线(fit_power_law 曾因 NaN 直接抛 LinAlgError 崩溃)
# ---------------------------------------------------------------------------

class TestNonFiniteGuard:
    """框架哲学"数据不足时显式降级,绝不硬输出"对 NaN/inf 同样适用。"""

    def test_fit_power_law_nan_feature_returns_graceful_failure(self):
        """特征列含 NaN 时应返回 success=False,而不是让 SVD 崩溃。"""
        df = pd.DataFrame({"a": [1.0, 2.0, 3.0, 4.0, 5.0],
                           "b": [1.0, 2.0, np.nan, 4.0, 5.0]})
        result = dr.fit_power_law(df, target="a", features=["b"])
        assert result["success"] is False
        assert "non-finite" in result["reason"] or "非有限" in result["reason"]

    def test_fit_power_law_inf_target_returns_graceful_failure(self):
        """目标列含 inf 时应返回 success=False。"""
        df = pd.DataFrame({"a": [1.0, 2.0, np.inf, 4.0, 5.0],
                           "b": [1.0, 2.0, 3.0, 4.0, 5.0]})
        result = dr.fit_power_law(df, target="a", features=["b"])
        assert result["success"] is False

    def test_fit_with_templates_nan_returns_graceful_failure(self, data):
        """模板路径同样要拦截 NaN。"""
        sys = make_system()
        Re = dr.Group({"rho": 1, "v": 1, "D": 1, "mu": -1}, "Re")
        Pr = dr.Group({"cp": 1, "mu": 1, "k": -1}, "Pr")
        Nu = dr.Group({"h": 1, "D": 1, "k": -1}, "Nu")
        bad = data.copy()
        bad.loc[3, "mu"] = np.nan
        result = dr.fit_with_templates(bad, sys, target_group=Nu,
                                       feature_groups=[Re, Pr])
        assert result["success"] is False

    def test_discover_with_nan_never_raises(self, data):
        """端到端:discover() 遇到 NaN 应显式降级,fit 为 None/失败,绝不抛异常。"""
        sys = make_system()
        bad = data.copy()
        bad.loc[10, "v"] = np.nan
        result = dr.discover(bad, sys, target="h")
        assert result["fit"] is None or result["fit"].get("success") is False


# ---------------------------------------------------------------------------
# P0 合流:函数族在 π 空间运行(discover method='multi' 与 discover_multi_pi)
# ---------------------------------------------------------------------------

class TestPiSpaceMultiFamily:
    def test_discover_method_multi_runs_families_in_pi_space(self, data):
        """method='multi' 应在 π 空间运行全部函数族,目标隔离仍然成立。"""
        sys = make_system()
        result = dr.discover(data, sys, target="h", method="multi")
        assert result["candidates"]
        assert {"power_law", "rational_log"} <= {c.family
                                                 for c in result["candidates"]}
        assert result["fit"] is not None and result["fit"]["success"]
        assert len([g for g in result["groups"] if "h" in g]) == 1

    def test_discover_multi_pi_returns_candidates_and_best(self, data):
        """discover_multi_pi 应返回同构候选列表与最优候选,并升级到 Level 3。"""
        sys = make_system()
        result = dr.discover_multi_pi(data, sys, target="h")
        assert result["candidates"]
        assert result["best"] is not None
        assert result["data_report"]["level"] >= 3
        assert len([g for g in result["groups"] if "h" in g]) == 1

    def test_discover_multi_pi_power_law_recovers_db(self, data):
        """π 空间中 PowerLawFamily 仍应恢复 Dittus-Boelter(R² 极高)。"""
        sys = make_system()
        result = dr.discover_multi_pi(data, sys, target="h")
        pl = next(c for c in result["candidates"]
                  if c.family == "power_law")
        assert pl.success
        assert pl.r_squared > 0.99

    def test_discover_multi_pi_degrades_when_missing_variable(self, data):
        """缺列时 discover_multi_pi 同样显式降级,绝不硬输出。"""
        sys = make_system()
        bad = data.drop(columns=["cp"])
        result = dr.discover_multi_pi(bad, sys, target="h")
        assert result["candidates"] == []
        assert result["note"]


# ---------------------------------------------------------------------------
# 逐特征降级与共享变量检测(空气动力学真实数据验证驱动的架构改进)
# ---------------------------------------------------------------------------

def make_tiny_system():
    """最小可控系统:a[m], b[m], c[s], u[m/s] → 2 个 π 群。

    零空间基(经目标隔离)为 {b/a} 与 {u·c/a}:特征群与目标群共享 a,
    是构造共享变量警告的最小案例。
    """
    sys = dr.System("tiny")
    sys.add_variable("a", "m", "geometry")
    sys.add_variable("b", "m", "geometry")
    sys.add_variable("c", "s", "state")
    sys.add_variable("u", "m/s", "target")
    return sys


def make_tiny_data():
    n = 40  # ≥ 10 × 4 个变量,满足 check_data Level 2
    a = np.linspace(1.0, 2.0, n)
    b = np.linspace(1.0, 3.0, n)
    c = np.linspace(0.5, 2.0, n)
    u = 3.0 * a / c
    return pd.DataFrame({"a": a, "b": b, "c": c, "u": u})


class TestPerFeatureDegradation:
    """零点(合法物理值)只剔除受影响的特征群;NaN 仍维持整体失败。"""

    def test_zero_in_feature_variable_drops_only_affected_groups(self, data):
        """含零变量所在的所有特征群被显式剔除并留痕,其余群继续拟合。"""
        sys = make_system()
        groups = sys.dimensionless_groups(target="h")
        y_idx = next(i for i, g in enumerate(groups) if "h" in g)
        feat_vars = set().union(
            *[set(g) for i, g in enumerate(groups) if i != y_idx])
        candidates = sorted(feat_vars - set(groups[y_idx]) - {"h"})
        assert candidates, "需要一个只出现在特征群中的变量来构造测试"
        zeroed = candidates[0]

        bad = data.copy()
        bad[zeroed] = 0.0
        result = dr.discover_multi_pi(bad, sys, target="h")

        assert result["dropped_features"], "含零的特征群应被显式剔除"
        assert any("剔除" in w for w in result["data_report"]["warnings"])
        pl = next((c for c in result["candidates"]
                   if c.family == "power_law"), None)
        assert pl is not None and pl.success, "剩余特征群应继续完成拟合"

    def test_all_features_dropped_degrades_explicitly(self):
        """全部特征群含非正值时显式降级,绝不硬输出。"""
        sys = make_tiny_system()
        bad = make_tiny_data()
        bad["b"] = 0.0  # 唯一特征群 b/a 变为零
        result = dr.discover(bad, sys, target="u")
        assert result["fit"] is None
        assert result["note"]
        assert result["dropped_features"]

    def test_nan_feature_still_fails_globally(self, data):
        """不对称设计:NaN 是数据质量事故,维持整体失败,不做逐特征降级。"""
        sys = make_system()
        bad = data.copy()
        bad.loc[10, "v"] = np.nan
        result = dr.discover(bad, sys, target="h")
        assert result["fit"] is None or result["fit"].get("success") is False
        assert result["dropped_features"] == []


class TestLevel3UnifiedGating:
    """Level 3 只认 out-of-fold CV,不再用样本内 R² 升级。"""

    def test_promote_uses_cv_threshold(self):
        from dimreflex.discover import _promote_level_3
        report = {"level": 0, "warnings": []}
        _promote_level_3(report, 0.95, in_sample_r_squared=0.99)
        assert report["level"] == 3
        assert any("k-fold" in warning for warning in report["warnings"])

    def test_promote_refuses_when_cv_low(self):
        from dimreflex.discover import _promote_level_3
        report = {"level": 0, "warnings": []}
        _promote_level_3(report, 0.85, in_sample_r_squared=0.99)
        assert report["level"] == 0
        assert any("未升级" in warning for warning in report["warnings"])

    def test_promote_refuses_when_cv_unavailable(self):
        from dimreflex.discover import _promote_level_3
        report = {"level": 0, "warnings": []}
        _promote_level_3(report, None, in_sample_r_squared=0.99)
        assert report["level"] == 0
        assert report["warnings"] == []


class TestStdErrorsParity:
    """fit_with_templates 与 fit_power_law 的统计输出必须对等。

    背景:泵相似律正对照(validate_pump.py)时模板路径缺 std_errors,
    脚本被迫手算兜底(待办 2)。本类锁住补丁:同一(群定义等价的)设计
    矩阵下,两条路径的指数与标准误逐位一致(机器精度)。"""

    def test_templates_match_power_law(self, data):
        sys = make_system()
        Re = dr.Group({"rho": 1, "v": 1, "D": 1, "mu": -1}, "Re")
        Pr = dr.Group({"cp": 1, "mu": 1, "k": -1}, "Pr")
        Nu = dr.Group({"h": 1, "D": 1, "k": -1}, "Nu")
        tmpl = dr.fit_with_templates(data, sys, target_group=Nu,
                                     feature_groups=[Re, Pr])
        pl = dr.discover(data, sys, target="h", method="power_law")["fit"]

        assert tmpl["success"] and pl["success"]
        for key in ("std_errors", "exponents"):
            tm_vals = sorted(tmpl[key].values())
            pl_vals = sorted(pl[key].values())
            for a, b in zip(tm_vals, pl_vals, strict=True):
                assert np.isclose(a, b, rtol=1e-9)
        assert tmpl["identifiable"]
        assert tmpl["design_rank"] == tmpl["design_columns"] == 3
        assert np.isfinite(tmpl["condition_number"])


class TestSharedVariableDetection:
    """目标群与特征群共享变量 → 显式提示 R² 可能虚高(π 群伪相关)。"""

    def test_auto_path_flags_shared_variables(self):
        """tiny 系统的特征群 b/a 与目标群 u·c/a 共享 a,应被标注。"""
        sys = make_tiny_system()
        result = dr.discover(make_tiny_data(), sys, target="u")
        assert result["shared_vars"], "共享变量应被检出"
        shared_flat = {v for vals in result["shared_vars"].values()
                       for v in vals}
        assert "a" in shared_flat
        assert any("共享变量" in w
                   for w in result["data_report"]["warnings"])

    def test_templates_flag_shared_variables(self, data):
        """教科书情形:Nu 与 Re 共享 D、Nu 与 Pr 共享 k,应给出提示。"""
        sys = make_system()
        Re = dr.Group({"rho": 1, "v": 1, "D": 1, "mu": -1}, "Re")
        Pr = dr.Group({"cp": 1, "mu": 1, "k": -1}, "Pr")
        Nu = dr.Group({"h": 1, "D": 1, "k": -1}, "Nu")
        result = dr.fit_with_templates(data, sys, target_group=Nu,
                                       feature_groups=[Re, Pr])
        assert result["success"]
        assert any("共享变量" in w for w in result["warnings"])

    def test_templates_no_warning_when_disjoint(self):
        """变量集不相交时不应有共享变量警告。"""
        sys = dr.System("disjoint")
        sys.add_variable("a", "m", "geometry")
        sys.add_variable("b", "m", "geometry")
        sys.add_variable("c", "s", "state")
        sys.add_variable("e", "kg", "material_property")
        sys.add_variable("f", "kg/s", "target")
        n = 50
        df = pd.DataFrame({
            "a": np.linspace(1.0, 2.0, n), "b": np.linspace(2.0, 4.0, n),
            "c": np.linspace(0.5, 1.5, n), "e": np.linspace(0.8, 1.2, n),
            "f": np.linspace(1.0, 3.0, n)})
        target = dr.Group({"f": 1, "c": 1, "e": -1}, "flow_ratio")
        feature = dr.Group({"b": 1, "a": -1}, "aspect")
        result = dr.fit_with_templates(df, sys, target_group=target,
                                       feature_groups=[feature])
        assert result["success"]
        assert result["warnings"] == []
