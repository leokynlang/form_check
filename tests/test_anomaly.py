"""anomaly.py 单元与集成测试:三层检测 + 三种调整策略 + 审计。"""

import dimreflex as dr
import numpy as np
import pandas as pd
from dimreflex import anomaly


def make_data(n=60, seed=7, outliers=None, mult=20.0, noise=0.02):
    """合成幂律数据 y = 2·x^0.8 (+2% 噪声), 可在指定索引注入离群点。"""
    rng = np.random.default_rng(seed)
    x = np.exp(rng.uniform(np.log(1.0), np.log(100.0), n))
    y = 2.0 * x ** 0.8 * rng.normal(1.0, noise, n)
    if outliers:
        for i in outliers:
            y[i] *= mult
    return pd.DataFrame({'x': x, 'y': y})


class TestDataHealth:
    def test_clean_data_no_issues(self):
        issues = anomaly.check_data_health(make_data(), 'y', ['x'])
        assert issues == []

    def test_non_positive_flagged_fatal(self):
        df = make_data()
        df.loc[3, 'y'] = -1.0
        issues = anomaly.check_data_health(df, 'y', ['x'])
        assert any(i['type'] == 'non_positive' and i['severity'] == 'fatal'
                   for i in issues)

    def test_zero_variance_flagged(self):
        df = make_data()
        df['c'] = 5.0
        issues = anomaly.check_data_health(df, 'y', ['x', 'c'])
        assert any(i['type'] == 'zero_variance' for i in issues)

    def test_missing_column(self):
        issues = anomaly.check_data_health(make_data(), 'y', ['nope'])
        assert any(i['type'] == 'missing_column' for i in issues)

    def test_non_finite_flagged(self):
        df = make_data()
        df.loc[2, 'x'] = np.nan
        issues = anomaly.check_data_health(df, 'y', ['x'])
        assert any(i['type'] == 'non_finite' for i in issues)


class TestDetection:
    def test_clean_residuals_no_outliers(self):
        df = make_data(seed=11)
        g = anomaly.fit_with_guard(df, 'y', ['x'], leverage_factor=3.0)
        assert g['success']
        assert g['anomalies']['residual_outliers']['indices'] == []
        assert g['after'] is None
        assert '无需调整' in g['audit']['note']

    def test_outlier_detected(self):
        df = make_data(seed=11, outliers=[5])
        g = anomaly.fit_with_guard(df, 'y', ['x'], leverage_factor=3.0)
        assert 5 in g['anomalies']['residual_outliers']['indices']
        assert g['audit']['flagged_sources']['5'] == 'residual'

    def test_leverage_detected_without_residual_outlier(self):
        """杠杆点在曲线上(残差小)也应被几何层抓住——与残差离群正交。"""
        df = make_data(seed=11)
        df.loc[7, 'x'] = 1.0e4
        df.loc[7, 'y'] = 2.0 * df.loc[7, 'x'] ** 0.8  # 恰好在曲线上
        g = anomaly.fit_with_guard(df, 'y', ['x'], leverage_factor=3.0)
        assert 7 in g['anomalies']['leverage']['indices']
        assert 7 not in g['anomalies']['residual_outliers']['indices']

    def test_collinearity_flagged_on_duplicate_feature(self):
        df = make_data(seed=11)
        df['x2'] = df['x']
        g = anomaly.fit_with_guard(df, 'y', ['x', 'x2'])
        assert g['anomalies']['collinearity'] is not None
        assert g['anomalies']['collinearity']['condition_number'] > 30

    def test_leverage_undetectable_when_n_too_small(self):
        df = make_data(n=3)
        X = np.column_stack([np.ones(3), np.log(df['x'].to_numpy())])
        lev = anomaly.detect_leverage(X)
        assert lev['indices'] == []
        assert '不可辨识' in lev['note']


class TestAdjustments:
    def test_robust_closer_to_truth_than_ols(self):
        df = make_data(seed=11)
        df.loc[df['x'].idxmax(), 'y'] *= 50.0  # 高影响位置的大离群
        g = anomaly.fit_with_guard(df, 'y', ['x'], adjustment='robust',
                                   leverage_factor=3.0)
        assert g['success'] and g['after'] is not None
        b_ols = g['before']['params']['exponents']['x']
        b_rob = g['after']['params']['exponents']['x']
        assert abs(b_rob - 0.8) < abs(b_ols - 0.8)
        assert g['after']['params']['n_downweighted'] >= 1
        # 样本内 R² 是 OLS 按构造最优的, 不能用来评价稳健修正;
        # 有意义的证据是: 稳健拟合相对 OLS 的系数位移方向正确
        assert g['comparison']['max_exponent_shift'] > 0

    def test_exclude_removes_injected_point_with_audit(self):
        df = make_data(seed=11, outliers=[5])
        g = anomaly.fit_with_guard(df, 'y', ['x'], adjustment='exclude',
                                   leverage_factor=3.0)
        assert g['after'] is not None
        assert g['after']['family'] == 'power_law_excluded'
        assert g['after']['n_samples'] == len(df) - 1
        exc = g['audit']['excluded']
        assert [e['index'] for e in exc] == [5]
        assert exc[0]['residual_z'] is not None
        assert g['comparison']['r2_after'] > g['comparison']['r2_before']

    def test_exclude_refuses_mass_flagging(self):
        """1/3 以上样本被标记 → 拒绝自动剔除(多半是模型形式错了)。"""
        df = make_data(n=60, seed=11)
        corrupt = list(range(25))
        g0 = anomaly.fit_with_guard(df, 'y', ['x'], adjustment='exclude',
                                    mad_threshold=2.0, leverage_factor=3.0)
        df_corrupt = df.copy()
        df_corrupt.loc[corrupt, 'y'] *= 20.0
        g = anomaly.fit_with_guard(df_corrupt, 'y', ['x'], adjustment='exclude',
                                   mad_threshold=2.0, leverage_factor=3.0)
        assert 'refused' in g['audit']
        assert g['after'] is None
        assert g0['success']  # 对照: 正常数据不触发拒绝

    def test_report_only_never_adjusts(self):
        df = make_data(seed=11, outliers=[5])
        g = anomaly.fit_with_guard(df, 'y', ['x'])  # 默认 report_only
        assert g['adjustment'] == 'report_only'
        assert g['after'] is None
        assert g['anomalies']['residual_outliers']['indices'] == [5]

    def test_exclude_leverage_semantics(self):
        """杠杆点(在曲线上)默认不剔除, 只报告; exclude_leverage=True 才剔除。"""
        df = make_data(seed=11)
        df.loc[7, 'x'] = 1.0e4
        df.loc[7, 'y'] = 2.0 * df.loc[7, 'x'] ** 0.8  # 杠杆点, 残差正常

        g = anomaly.fit_with_guard(df, 'y', ['x'], adjustment='exclude')
        assert g['after'] is None
        assert 7 in g['anomalies']['leverage']['indices']
        assert 'exclude_leverage=True' in g['audit']['note']

        g2 = anomaly.fit_with_guard(df, 'y', ['x'], adjustment='exclude',
                                    exclude_leverage=True)
        assert g2['after'] is not None
        assert [e['index'] for e in g2['audit']['excluded']] == [7]
        assert g2['audit']['excluded'][0]['source'] == '+leverage'

    def test_robust_not_blocked_by_leverage_only(self):
        """robust 是全局重拟合策略, 不被点级剔除语义早退拦截。"""
        df = make_data(seed=11)
        df.loc[7, 'x'] = 1.0e4
        df.loc[7, 'y'] = 2.0 * df.loc[7, 'x'] ** 0.8
        g = anomaly.fit_with_guard(df, 'y', ['x'], adjustment='robust')
        assert g['after'] is not None
        assert g['audit']['flagged_indices'] == []

    def test_unknown_adjustment_explicit_failure(self):
        df = make_data()
        g = anomaly.fit_with_guard(df, 'y', ['x'], adjustment='delete')
        assert not g['success']
        assert 'report_only' in g['reason']

    def test_fatal_health_blocks_fitting(self):
        df = make_data()
        df.loc[3, 'x'] = -2.0
        g = anomaly.fit_with_guard(df, 'y', ['x'], adjustment='robust')
        assert not g['success']
        assert any(i['severity'] == 'fatal' for i in g['health'])


class TestIRLS:
    def test_irls_converges_and_downweights(self):
        rng = np.random.default_rng(3)
        x = np.exp(rng.uniform(0, 3, 80))
        y = np.log(2.0 * x ** 0.8) + rng.normal(0, 0.02, 80)
        y[10] += 5.0
        X = np.column_stack([np.ones(80), np.log(x)])
        ir = anomaly.robust_irls(X, y)
        assert ir['converged']
        assert ir['weights'][10] < 0.5
        assert abs(ir['coef'][1] - 0.8) < 0.02


class TestPiSpaceIntegration:
    def test_guard_wraps_discover_pipeline(self):
        """守卫层与发现管线的集成: 对 π 空间表做守卫式拟合。"""
        from pathlib import Path
        data = pd.read_csv(Path(__file__).resolve().parent.parent
                           / 'heat_exchanger_data.csv')
        s = dr.System('tube')
        s.add_variable('rho', 'kg/m**3', 'material_property')
        s.add_variable('mu', 'Pa*s', 'material_property')
        s.add_variable('cp', 'J/(kg*K)', 'material_property')
        s.add_variable('k', 'W/(m*K)', 'material_property')
        s.add_variable('D', 'm', 'geometry')
        s.add_variable('v', 'm/s', 'state')
        s.add_variable('h', 'W/(m**2*K)', 'target')

        res = dr.discover(data, s, target='h')
        assert res['pi_data'] is not None
        y_i = next(i for i, g in enumerate(res['groups']) if 'h' in g)
        y_col = res['pi_data'].columns[y_i]
        x_cols = [c for c in res['pi_data'].columns if c != y_col]

        g = anomaly.fit_with_guard(res['pi_data'], y_col, x_cols,
                                   adjustment='robust')
        assert g['success']
        assert g['before']['success']
        assert g['after'] is not None
        assert g['after']['success']
        # 2% 噪声的干净合成数据: 系数位移应很小
        assert g['comparison']['max_exponent_shift'] < 0.05
