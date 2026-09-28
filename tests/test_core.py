"""core.py 单元测试:变量、量纲向量、π 群与目标隔离。"""

import dimreflex as dr
import numpy as np
import pytest


def make_system():
    """换热器管侧对流换热系统(与 data_analysis.py 一致)。"""
    sys = dr.System("tube_side_heat_transfer")
    sys.add_variable("rho", "kg/m**3", "material_property")
    sys.add_variable("mu", "Pa*s", "material_property")
    sys.add_variable("cp", "J/(kg*K)", "material_property")
    sys.add_variable("k", "W/(m*K)", "material_property")
    sys.add_variable("D", "m", "geometry")
    sys.add_variable("v", "m/s", "state")
    sys.add_variable("h", "W/(m**2*K)", "target")
    return sys


class TestVariable:
    def test_dimension_vector_of_velocity(self):
        v = dr.Variable("v", "m/s", "state")
        assert np.allclose(v.dimension_vector(), [0, 1, -1, 0, 0, 0, 0])

    def test_dimension_vector_of_density(self):
        rho = dr.Variable("rho", "kg/m**3", "material_property")
        assert np.allclose(rho.dimension_vector(), [1, -3, 0, 0, 0, 0, 0])

    def test_invalid_unit_raises(self):
        with pytest.raises(ValueError):
            dr.Variable("x", "not_a_unit_!!", "state")

    def test_invalid_role_raises(self):
        with pytest.raises(ValueError):
            dr.Variable("x", "m", "bad_role")


class TestSystem:
    def test_duplicate_variable_raises(self):
        sys = dr.System("s")
        sys.add_variable("rho", "kg/m**3", "material_property")
        with pytest.raises(ValueError):
            sys.add_variable("rho", "kg/m**3", "material_property")

    def test_dimension_matrix_shape(self):
        sys = make_system()
        D, names = sys.dimension_matrix()
        assert D.shape == (7, 7)
        assert names == list(sys.variables.keys())

    def test_pi_group_count(self):
        # 7 变量,量纲矩阵秩 4 → 3 个无量纲群 (Buckingham π: n - r)
        groups = make_system().dimensionless_groups()
        assert len(groups) == 3

    def test_groups_are_dimensionless(self):
        sys = make_system()
        D, names = sys.dimension_matrix()
        for g in sys.dimensionless_groups():
            vec = np.array([g.get(n, 0.0) for n in names], dtype=float)
            assert np.allclose(D @ vec, 0, atol=1e-6)

    def test_canonicalization_never_breaks_dimensionless_constraint(self):
        """回归:旧的逐分量取整会让这个矩阵产生约 0.049 的量纲残差。"""
        dimension_matrix = np.array([
            [-2, 0, 3, 2, 1, 1, -3],
            [2, 0, -3, -2, 0, 2, 0],
            [1, 3, 1, -1, 1, 1, -3],
            [-3, 1, -1, 1, -1, -2, -2],
        ])
        bases = ["kg", "m", "s", "K"]
        sys = dr.System("rounding_regression")
        for col in range(dimension_matrix.shape[1]):
            unit = "*".join(
                f"{bases[row]}**({int(dimension_matrix[row, col])})"
                for row in range(dimension_matrix.shape[0])
            )
            sys.add_variable(f"v{col}", unit, "state")

        D, names = sys.dimension_matrix()
        groups = sys.dimensionless_groups()

        assert len(groups) == len(names) - np.linalg.matrix_rank(D)
        for group in groups:
            vector = np.array([group.get(name, 0.0) for name in names])
            assert np.allclose(D @ vector, 0.0, rtol=0.0, atol=1e-9)
            assert all(float(exp).is_integer() for exp in group.values())

    def test_target_isolation_preserves_dimensionless_constraint(self):
        sys = make_system()
        D, names = sys.dimension_matrix()
        for group in sys.dimensionless_groups(target="h"):
            vector = np.array([group.get(name, 0.0) for name in names])
            assert np.allclose(D @ vector, 0.0, rtol=0.0, atol=1e-9)

    def test_target_isolation(self):
        """target 隔离后:恰好一个群含 target 且指数为 1,其余群不含 target。"""
        sys = make_system()
        groups = sys.dimensionless_groups(target="h")
        target_groups = [g for g in groups if "h" in g]
        assert len(target_groups) == 1
        assert target_groups[0]["h"] == pytest.approx(1.0)
        assert target_groups[0] == {"k": -1.0, "D": 1.0, "h": 1.0}

    def test_target_not_in_system_raises(self):
        sys = make_system()
        with pytest.raises(ValueError):
            sys.dimensionless_groups(target="does_not_exist")
