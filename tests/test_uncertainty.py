"""uncertainty.py(蒙特卡洛传播)单元测试。

背景:Variable.uncertainty 自 v0.1 声明、从未接入(待办 3)。
本模块把测量噪声传播成指数的置信区间,与回归 std_errors 分开报告。
"""

import dimreflex as dr
import numpy as np
import pandas as pd
from dimreflex import propagate_uncertainty


def make_system(ux=0.02, uy=0.01):
    s = dr.System("mc_test")
    # 量纲冗余是 π 分析的前提:三个长度量纲变量 → 零空间维 2 → 两个群。
    # (初版用 x[m]/y[s] 量纲独立,零空间维 0,连群都不存在——教训:
    # 玩具系统先手算零空间维数再写测试。)
    s.add_variable("x", "m", "state", uncertainty=ux)
    s.add_variable("L", "m", "geometry")
    s.add_variable("y", "m", "target", uncertainty=uy)
    return s


def make_data(n=80, seed=5):
    # seed 参数保留给将来加噪声的变体;当前数据是确定性的
    x = np.logspace(1, 4, n)
    L = 0.1
    y = 0.5 * L * (x / L) ** 0.73     # 真值: y/L = 0.5*(x/L)^0.73
    return pd.DataFrame({"x": x, "L": np.full(n, L), "y": y})


class TestPropagateUncertainty:
    def test_basic_propagation(self):
        """有声明 → 成功;MC 区间罩住基准拟合指数。

        注意指数的真值是框架表示的函数:此玩具系统的 π 基选出目标群
        y/x、特征群 x/L,而 y/x = 0.5*(x/L)^(0.73-1),故对 x/L 的
        真指数是 -0.27 而非 0.73——断言对照基准拟合值(框架自己的
        表示),不对照手写假设(与负对照"判据建在框架表示上"同款纪律)。"""
        res = propagate_uncertainty(make_data(), make_system(), target="y",
                                    n_sims=200)
        assert res["success"], res["reason"]
        assert set(res["variables_perturbed"]) == {"x", "y"}
        base_b = next(iter(res["fit"]["exponents"].values()))
        ex = next(iter(res["exponents_mc"].values()))
        assert ex["q025"] < base_b < ex["q975"]
        # 测量噪声小(u=2%/1%),MC 均值不应漂离基准拟合超过 3 个 MC 标准差
        assert abs(ex["mean"] - base_b) < 3 * ex["std"]
        assert 0 < ex["std"] < 0.5

    def test_no_declaration_refused(self):
        """没声明 uncertainty → 结构化拒绝,不虚构传播结果。"""
        res = propagate_uncertainty(make_data(), make_system(ux=None, uy=None),
                                    target="y")
        assert not res["success"]
        assert "uncertainty" in res["reason"]

    def test_invalid_uncertainty_refused(self):
        res = propagate_uncertainty(make_data(), make_system(ux=-0.1),
                                    target="y")
        assert not res["success"]
        assert "非法" in res["reason"]

    def test_missing_column_refused(self):
        """声明了不确定度的变量缺列 → 门卫拦截(不抛 KeyError)。"""
        s = make_system()
        s.add_variable("z", "K", "material_property", uncertainty=0.03)
        df = make_data()
        res = propagate_uncertainty(df, s, target="y")
        assert not res["success"]
        assert "缺列" in res["reason"]

    def test_deterministic(self):
        """同 seed → 区间逐位一致(可复现性纪律)。"""
        a = propagate_uncertainty(make_data(), make_system(), target="y")
        b = propagate_uncertainty(make_data(), make_system(), target="y")
        assert a["exponents_mc"] == b["exponents_mc"]

    def test_regression_se_reported_separately(self):
        """回归 std_errors 作为对照字段一并返回,语义区分写进 note。"""
        res = propagate_uncertainty(make_data(), make_system(), target="y")
        assert res["regression_std_errors"]
        assert "回归残差噪声" in res["note"]
