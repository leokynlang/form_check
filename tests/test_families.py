"""families.py 测试:基底类契约、统一候选模型 schema 与幂律族迁移。

对应《项目评价与改进建议》2.1(函数族扩展)与第四节(统一候选模型 schema、
注册表模式):所有函数族实现同一基底类,产出同构的 CandidateModel。
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from dimreflex.families import (
    FAMILIES,
    CandidateModel,
    FunctionFamily,
    PiecewisePowerLawFamily,
    PowerLawFamily,
    RationalLogFamily,
    SymbolicRegressionFamily,
    discover_multi,
    register_family,
)

DATA_PATH = Path(__file__).resolve().parent.parent / "heat_exchanger_data.csv"


@pytest.fixture(scope="module")
def data():
    return pd.read_csv(DATA_PATH)


@pytest.fixture(scope="module")
def nu_data():
    """纯幂律合成数据:Nu = 0.023 * Re^0.8 * Pr^0.4 (+2% 噪声)。"""
    rng = np.random.default_rng(7)
    n = 60
    Re = rng.uniform(5000, 100000, n)
    Pr = rng.uniform(1.0, 10.0, n)
    Nu = 0.023 * Re**0.8 * Pr**0.4 * rng.normal(1.0, 0.02, n)
    return pd.DataFrame({"Re": Re, "Pr": Pr, "Nu": Nu})


@pytest.fixture(scope="module")
def two_regime_data():
    """经典双工况合成数据:层流(b=1/3)与湍流(b=0.8),变点在 Re=2300。"""
    rng = np.random.default_rng(11)
    Re = np.concatenate([
        rng.uniform(500, 2200, 40),      # 层流段
        rng.uniform(2400, 50000, 40),    # 湍流段
    ])
    Nu = np.where(
        Re < 2300,
        0.5 * Re**(1 / 3),
        0.023 * Re**0.8,
    )
    Nu *= rng.normal(1.0, 0.01, len(Re))
    return pd.DataFrame({"Re": Re, "Nu": Nu})


# ---------------------------------------------------------------------------
# 统一候选模型 schema(文档第四节:各层产出同构候选)
# ---------------------------------------------------------------------------

class TestCandidateModel:
    REQUIRED_KEYS = frozenset({"family", "formula", "params", "success",
                               "r_squared", "n_params", "n_samples"})

    def test_to_dict_contains_required_keys(self, nu_data):
        cm = PowerLawFamily().fit(nu_data, target="Nu", features=["Re", "Pr"])
        d = cm.to_dict()
        missing = self.REQUIRED_KEYS - set(d)
        assert not missing, f"schema 缺少字段: {missing}"

    def test_to_dict_roundtrip(self, nu_data):
        cm = PowerLawFamily().fit(nu_data, target="Nu", features=["Re", "Pr"])
        d = cm.to_dict()
        cm2 = CandidateModel.from_dict(d)
        assert cm2.family == cm.family
        assert cm2.formula == cm.formula
        assert cm2.success == cm.success
        assert cm2.params == cm.params


# ---------------------------------------------------------------------------
# 基底类契约:fit() 是模板方法,统一兜底,绝不向调用方抛数学异常
# ---------------------------------------------------------------------------

class ExplodingFamily(FunctionFamily):
    """_fit 抛数学异常的家族——验证基类兜底。"""
    name = "exploding"

    def _fit(self, data, target, features):
        raise np.linalg.LinAlgError("SVD did not converge")


class SilentFailureFamily(FunctionFamily):
    """_fit 返回失败但不给 reason——基类应补默认原因。"""
    name = "silent"

    def _fit(self, data, target, features):
        return CandidateModel(family=self.name, formula="", params={},
                              success=False)


class TestFunctionFamilyContract:
    def test_power_law_is_function_family(self):
        assert isinstance(PowerLawFamily(), FunctionFamily)

    def test_fit_catches_math_exception(self):
        """子类崩溃时 fit() 必须返回 success=False,绝不抛出。"""
        df = pd.DataFrame({"x": [1.0, 2.0, 3.0], "y": [1.0, 2.0, 3.0]})
        cm = ExplodingFamily().fit(df, target="y", features=["x"])
        assert cm.success is False
        assert cm.family == "exploding"
        assert cm.reason  # 有可读的原因

    def test_fit_fills_default_reason(self):
        df = pd.DataFrame({"x": [1.0, 2.0, 3.0], "y": [1.0, 2.0, 3.0]})
        cm = SilentFailureFamily().fit(df, target="y", features=["x"])
        assert cm.success is False
        assert cm.reason


# ---------------------------------------------------------------------------
# PowerLawFamily:现有幂律逻辑迁移为第一个插件(行为不回退)
# ---------------------------------------------------------------------------

class TestPowerLawFamily:
    def test_recovers_dittus_boelter(self, nu_data):
        fam = PowerLawFamily()
        cm = fam.fit(nu_data, target="Nu", features=["Re", "Pr"])
        assert cm.success is True
        assert cm.family == "power_law"
        assert cm.r_squared > 0.99
        assert np.isfinite(cm.r_squared)
        assert cm.params["prefactor"] == pytest.approx(0.023, rel=0.1)
        assert cm.params["exponents"]["Re"] == pytest.approx(0.8, abs=0.05)
        assert cm.params["exponents"]["Pr"] == pytest.approx(0.4, abs=0.05)
        assert "Nu" in cm.formula
        assert cm.n_params == 3  # prefactor + 2 指数
        assert cm.n_samples == len(nu_data)

    def test_nonpositive_data_is_graceful(self):
        df = pd.DataFrame({"a": [1.0, 2.0, 3.0], "b": [1.0, 2.0, 0.0]})
        cm = PowerLawFamily().fit(df, target="a", features=["b"])
        assert cm.success is False
        assert "非正" in cm.reason

    def test_nan_data_is_graceful(self):
        """继承既有 NaN 防线:非有限值显式失败(回归保护)。"""
        df = pd.DataFrame({"a": [1.0, 2.0, 3.0, 4.0, 5.0],
                           "b": [1.0, 2.0, np.nan, 4.0, 5.0]})
        cm = PowerLawFamily().fit(df, target="a", features=["b"])
        assert cm.success is False
        assert "非有限" in cm.reason


# ---------------------------------------------------------------------------
# PiecewisePowerLawFamily(2.1 第一步:分段幂律 + 变点检测)
# ---------------------------------------------------------------------------

class TestPiecewisePowerLawFamily:
    def test_pure_power_law_reports_single_segment(self, nu_data):
        """纯幂律数据应如实报告单一分段——变点检测不得虚报。"""
        cm = PiecewisePowerLawFamily(max_breaks=1).fit(
            nu_data, target="Nu", features=["Re"])
        assert cm.success is True
        assert cm.params["n_segments"] == 1

    def test_finds_breakpoint_between_regimes(self, two_regime_data):
        """双工况数据应检出变点:间隙包含真值 2300,两段指数接近真值。

        变点可识别性:过渡发生在数据的无约束间隙(2200,2400)内,数据只能
        把变点定位到该间隙,报告值取间隙几何中点;断言针对可识别的范围。
        """
        cm = PiecewisePowerLawFamily(max_breaks=1).fit(
            two_regime_data, target="Nu", features=["Re"])
        assert cm.success is True
        assert cm.params["n_segments"] == 2
        bp_lo, bp_hi = cm.params["breakpoint_interval"][0]
        assert bp_lo < 2300 < bp_hi, (
            f"过渡区间 [{bp_lo}, {bp_hi}] 应包含真实转变点 2300")
        bp = cm.params["breakpoints"][0]
        assert bp == pytest.approx(np.sqrt(bp_lo * bp_hi))
        seg = cm.params["segments"][0]
        assert seg["exponent"] == pytest.approx(1 / 3, abs=0.05)
        assert seg["prefactor"] == pytest.approx(0.5, rel=0.15)
        seg2 = cm.params["segments"][1]
        assert seg2["exponent"] == pytest.approx(0.8, abs=0.05)
        assert seg2["prefactor"] == pytest.approx(0.023, rel=0.2)

    def test_r_squared_beats_single_power_law(self, two_regime_data):
        """分段模型在真实分段数据上必须优于全局单幂律(否则不该分段)。

        宽动态范围下 log 空间 R² 天然偏高,增益阈值取 0.01 而非 0.05。
        """
        single = PowerLawFamily().fit(two_regime_data, target="Nu",
                                      features=["Re"])
        piecewise = PiecewisePowerLawFamily(max_breaks=1).fit(
            two_regime_data, target="Nu", features=["Re"])
        assert piecewise.r_squared > single.r_squared + 0.01

    def test_formula_is_readable(self, two_regime_data):
        cm = PiecewisePowerLawFamily(max_breaks=1).fit(
            two_regime_data, target="Nu", features=["Re"])
        assert cm.success
        assert "Re <" in cm.formula and "Re ≥" in cm.formula

    def test_insufficient_data_is_graceful(self):
        """数据太少(不足拟合两段)应显式失败,不硬输出。"""
        df = pd.DataFrame({"Re": [100.0, 200.0, 300.0],
                           "Nu": [1.0, 2.0, 3.0]})
        cm = PiecewisePowerLawFamily(max_breaks=1).fit(
            df, target="Nu", features=["Re"])
        assert cm.success is False
        assert cm.reason

    def test_nan_data_is_graceful(self):
        df = pd.DataFrame({"Re": [1.0, 2.0, np.nan, 4.0],
                           "Nu": [1.0, 2.0, 3.0, 4.0]})
        cm = PiecewisePowerLawFamily(max_breaks=1).fit(
            df, target="Nu", features=["Re"])
        assert cm.success is False

    def test_multi_feature_two_regime_recovers_exponents(self):
        """多特征双工况:变点仍沿第一特征检出,每段恢复完整指数向量。

        Pr 在两段取不同指数(0.35 vs 0.4),验证段间指数允许不同——
        这是多特征分段相对"协变量全局共享"的独立信息。
        """
        rng = np.random.default_rng(41)
        n = 120
        Re = np.concatenate([
            rng.uniform(500, 2200, n // 2),      # 层流段
            rng.uniform(2400, 50000, n // 2),    # 湍流段
        ])
        Pr = rng.uniform(1.0, 10.0, n)
        Nu = np.where(
            Re < 2300,
            0.9 * Re**(1 / 3) * Pr**0.35,
            0.023 * Re**0.8 * Pr**0.4,
        )
        Nu *= rng.normal(1.0, 0.01, n)
        df = pd.DataFrame({"Re": Re, "Pr": Pr, "Nu": Nu})
        cm = PiecewisePowerLawFamily(max_breaks=1).fit(
            df, target="Nu", features=["Re", "Pr"])
        assert cm.success is True, cm.reason
        assert cm.params["n_segments"] == 2
        bp_lo, bp_hi = cm.params["breakpoint_interval"][0]
        assert bp_lo < 2300 < bp_hi
        seg = cm.params["segments"][0]
        assert seg["exponents"]["Re"] == pytest.approx(1 / 3, abs=0.06)
        assert seg["exponents"]["Pr"] == pytest.approx(0.35, abs=0.06)
        seg2 = cm.params["segments"][1]
        assert seg2["exponents"]["Re"] == pytest.approx(0.8, abs=0.05)
        assert seg2["exponents"]["Pr"] == pytest.approx(0.4, abs=0.05)

    def test_multi_feature_cv_prediction_path(self):
        """多特征分段模型经 predict_candidate 的 CV 预测链路可用。

        这是 autodiscovery 交叉验证与 monitoring 批次预测的实际路径;
        数据打乱行序后拟合,断言预测与真值吻合(覆盖段边界的两侧)。
        """
        rng = np.random.default_rng(43)
        n = 120
        Re = np.concatenate([
            rng.uniform(500, 2200, n // 2),
            rng.uniform(2400, 50000, n // 2),
        ])
        Pr = rng.uniform(1.0, 10.0, n)
        Nu = np.where(Re < 2300, 0.9 * Re**(1 / 3) * Pr**0.35,
                      0.023 * Re**0.8 * Pr**0.4)
        df = pd.DataFrame({"Re": Re, "Pr": Pr, "Nu": Nu})
        df = df.sample(frac=1.0, random_state=9).reset_index(drop=True)
        cm = PiecewisePowerLawFamily(max_breaks=1).fit(
            df, target="Nu", features=["Re", "Pr"])
        assert cm.success is True, cm.reason
        from dimreflex.autodiscovery import predict_candidate
        pred = predict_candidate(cm, df, ["Re", "Pr"])
        assert np.allclose(pred, df["Nu"].values, rtol=0.05)


# ---------------------------------------------------------------------------
# RationalLogFamily(2.1 第二步:log 空间有理函数,覆盖 Gnielinski 类)
# ---------------------------------------------------------------------------

class TestRationalLogFamily:
    @pytest.fixture(scope="module")
    def rational_data(self):
        """有理形式合成数据:log10(Nu) = (2.0 + 0.5·L) / (1 + 0.15·L) + 噪声,
        其中 L = log10(Re) - 3.7(在数据中心处归零,条件数友好)。"""
        rng = np.random.default_rng(23)
        Re = np.exp(rng.uniform(np.log(4000), np.log(100000), 70))
        L = np.log10(Re) - 3.7
        Nu = 10 ** ((2.0 + 0.5 * L) / (1 + 0.15 * L))
        Nu *= rng.normal(1.0, 0.005, len(Re))
        return pd.DataFrame({"Re": Re, "Nu": Nu})

    def test_recovers_rational_form(self, rational_data):
        """应恢复分子 2 项 / 分母 1 项的有理结构,且精度极高。

        R² 断言按噪声标定:0.5% 乘性噪声在 log10 空间约 ±0.0022,
        R² 的可达天花板约 0.998,故断言 > 0.998;结构/数值精度由
        predict_log10 与真值函数的逐点比较保证(下一条断言)。
        线性化带来的系数衰减(q: 0.136 vs 真值 0.15)是已声明的
        errors-in-variables 限制。
        """
        cm = RationalLogFamily(num_degree=2, den_degree=1).fit(
            rational_data, target="Nu", features=["Re"])
        assert cm.success is True
        assert cm.r_squared > 0.998
        assert np.isfinite(cm.r_squared)
        assert cm.params["num_degree"] == 2
        assert cm.params["den_degree"] == 1
        assert len(cm.params["num_coefs"]) == 2
        assert len(cm.params["den_coefs"]) == 1
        # 恢复的有理函数在数据中位数处应接近真值函数。
        # 坐标系统一:predict_log10 接受的是拟合内部中心化的 L
        # (log_center = mean(log10 Re)),真值公式用 fixture 的中心 3.7,
        # 两边都换算到同一 Re 点再比较。
        Re_mid = rational_data["Re"].median()
        L0_fit = cm.params["log_center"]
        pred_val = cm.params["predict_log10"](np.log10(Re_mid) - L0_fit)
        L_true = np.log10(Re_mid) - 3.7
        true_val = (2.0 + 0.5 * L_true) / (1 + 0.15 * L_true)
        assert pred_val == pytest.approx(true_val, abs=0.02)
        # 中心化参数应等于数据 log10 均值
        assert L0_fit == pytest.approx(float(np.log10(rational_data["Re"]).mean()),
                                       abs=1e-9)

    def test_aic_prefers_correct_degree(self, rational_data):
        """过高的分母阶数应在 AIC 上更差(无必要就不该加复杂度)。"""
        fit_21 = RationalLogFamily(num_degree=2, den_degree=1).fit(
            rational_data, target="Nu", features=["Re"])
        fit_33 = RationalLogFamily(num_degree=3, den_degree=3).fit(
            rational_data, target="Nu", features=["Re"])
        assert fit_21.success and fit_33.success
        assert fit_21.params["aic"] <= fit_33.params["aic"]

    def test_pure_power_law_still_fits(self):
        """幂律数据(分子 2 项/分母 0 项的退化有理式)也应能拟合。

        用单特征数据:有理族按声明只用第一特征,多特征数据的其余
        方差无法解释(那是 PowerLawFamily 的职责)。
        """
        rng = np.random.default_rng(31)
        Re = np.exp(rng.uniform(np.log(5000), np.log(100000), 60))
        Nu = 0.023 * Re**0.8 * rng.normal(1.0, 0.01, 60)
        df = pd.DataFrame({"Re": Re, "Nu": Nu})
        cm = RationalLogFamily(num_degree=2, den_degree=0).fit(
            df, target="Nu", features=["Re"])
        assert cm.success is True
        assert cm.r_squared > 0.99

    def test_insufficient_data_is_graceful(self):
        df = pd.DataFrame({"x": [1.0, 2.0, 3.0], "y": [1.0, 2.0, 3.0]})
        cm = RationalLogFamily(num_degree=2, den_degree=1).fit(
            df, target="y", features=["x"])
        assert cm.success is False
        assert cm.reason

    def test_multi_feature_covariate_recovery(self):
        """多特征:有理主特征 + 幂律协变量,交替最小二乘应恢复两者。

        真值:log10(Nu) = (2.0+0.5L)/(1+0.15L) + 0.7·(log10(Pr)−c)。
        断言协变量系数恢复到 ±0.05,且整体 R² 达到噪声允许的天花板。
        """
        rng = np.random.default_rng(37)
        n = 70
        Re = np.exp(rng.uniform(np.log(4000), np.log(100000), n))
        Pr = np.exp(rng.uniform(np.log(0.5), np.log(20.0), n))
        L = np.log10(Re) - 3.7
        log10y = ((2.0 + 0.5 * L) / (1 + 0.15 * L)
                  + 0.7 * (np.log10(Pr) - np.log10(Pr).mean()))
        Nu = 10 ** log10y * rng.normal(1.0, 0.005, n)
        df = pd.DataFrame({"Re": Re, "Pr": Pr, "Nu": Nu})
        cm = RationalLogFamily(num_degree=2, den_degree=1).fit(
            df, target="Nu", features=["Re", "Pr"])
        assert cm.success is True, cm.reason
        assert cm.r_squared > 0.99
        assert cm.params["extra_features"] == ["Pr"]
        assert cm.params["extra_coefs"][0] == pytest.approx(0.7, abs=0.05)
        # predict_candidate 是 CV/监控的实际预测路径
        from dimreflex.autodiscovery import predict_candidate
        pred = predict_candidate(cm, df, ["Re", "Pr"])
        assert np.allclose(pred, Nu, rtol=0.02)
        # 单特征时代的 predict_log10 契约在多特征下仍可用(带 Zv)
        L0 = cm.params["log_center"]
        z_mid = np.log10(Pr[n // 2]) - cm.params["extra_centers"][0]
        manual = (cm.params["predict_log10"](np.log10(Re[n // 2]) - L0,
                                             np.array([z_mid])))
        assert manual == pytest.approx(log10y[n // 2], abs=0.02)

    def test_nan_is_graceful(self):
        df = pd.DataFrame({"x": [1.0, 2.0, np.nan], "y": [1.0, 2.0, 3.0]})
        cm = RationalLogFamily(num_degree=2, den_degree=1).fit(
            df, target="y", features=["x"])
        assert cm.success is False


# ---------------------------------------------------------------------------
# SymbolicRegressionFamily(2.1 第三步:PySR 占位——显式失败优于假装实现)
# ---------------------------------------------------------------------------

class TestSymbolicRegressionFamily:
    def test_unavailable_reports_honest_failure(self, nu_data):
        """PySR 未安装时应显式失败并说明原因,绝不假装成功。"""
        cm = SymbolicRegressionFamily().fit(nu_data, target="Nu",
                                            features=["Re", "Pr"])
        assert cm.success is False
        assert cm.reason  # 原因里应说明 PySR 不可用
        assert "PySR" in cm.reason or "pysr" in cm.reason.lower()

    def test_is_function_family(self):
        assert isinstance(SymbolicRegressionFamily(), FunctionFamily)


# ---------------------------------------------------------------------------
# 注册表 + discover_multi(文档第四节:注册表模式,各层产出同构候选)
# ---------------------------------------------------------------------------

class TestRegistryAndDiscoverMulti:
    def test_registry_has_builtins(self):
        assert "power_law" in FAMILIES
        assert "piecewise_power_law" in FAMILIES
        assert "rational_log" in FAMILIES
        assert "symbolic" in FAMILIES

    def test_register_custom_family(self):
        class AlwaysFails(FunctionFamily):
            name = "always_fails_test"
            def _fit(self, data, target, features):
                raise RuntimeError("boom")

        register_family(AlwaysFails)
        assert "always_fails_test" in FAMILIES

    def test_discover_multi_runs_all_families(self, nu_data):
        """所有注册家族都应产出同构候选,且崩溃家族不中断流程。"""
        candidates = discover_multi(nu_data, target="Nu",
                                    features=["Re", "Pr"])
        names = {c.family for c in candidates}
        assert {"power_law", "piecewise_power_law", "rational_log",
                "symbolic"} <= names
        for c in candidates:
            d = c.to_dict()  # 同构:每个候选都能走统一 schema 序列化
            assert CandidateModel.from_dict(d).family == c.family

    def test_discover_multi_survives_exploding_family(self, nu_data):
        """注册一个会崩溃的家族,discover_multi 仍应完整返回其他候选。"""
        class Exploding(FunctionFamily):
            name = "exploding_multi_test"
            def _fit(self, data, target, features):
                raise np.linalg.LinAlgError("SVD again")

        register_family(Exploding)
        candidates = discover_multi(nu_data, target="Nu",
                                    features=["Re", "Pr"])
        names = {c.family for c in candidates}
        assert "exploding_multi_test" in names
        failed = next(c for c in candidates
                      if c.family == "exploding_multi_test")
        assert failed.success is False
        assert failed.reason  # 崩溃被折叠为带原因的失败候选

    def test_discover_multi_family_filter(self, nu_data):
        candidates = discover_multi(nu_data, target="Nu",
                                    features=["Re", "Pr"],
                                    families=["power_law"])
        assert [c.family for c in candidates] == ["power_law"]

    def test_discover_multi_sorts_by_r2(self, nu_data):
        candidates = discover_multi(nu_data, target="Nu",
                                    features=["Re", "Pr"])
        successes = [c for c in candidates if c.success]
        r2s = [c.r_squared for c in successes]
        assert r2s == sorted(r2s, reverse=True)
