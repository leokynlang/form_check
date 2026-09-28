"""
dimreflex.core
==============

核心抽象:
- Variable: 一个带单位、量纲、角色、误差范围的变量对象
- System: 一组变量的集合,支持单位检查、量纲矩阵构建、无量纲群生成

设计哲学:
- 变量必须显式声明单位和角色,框架不猜
- 量纲检查在数据进入框架的最早阶段进行
- 数据不足时,框架显式报告而不是默默处理
"""

from dataclasses import dataclass
from fractions import Fraction
from math import gcd, lcm

import numpy as np
import pint
from scipy.linalg import null_space

# 全局 unit registry——所有变量共享一个,避免单位转换混乱
_ureg = pint.UnitRegistry()


# 七个 SI 基本量纲。变量的量纲会被表达为这七个的整数指数向量。
# 注意:这里我们处理的是"量纲"而不是"单位"。kg/m^3 和 g/cm^3 单位不同但量纲相同。
SI_BASE_DIMS = ['[mass]', '[length]', '[time]', '[temperature]',
                '[current]', '[substance]', '[luminosity]']


def _rational_null_space(
        matrix: np.ndarray,
        max_denominator: int = 10000,
) -> list[list[Fraction]] | None:
    """Return an exact rational basis for the matrix null space when possible.

    Pint dimensionalities are normally integers or small rational numbers.
    Floating-point SVD finds the correct subspace, but its orthonormal basis is
    an arbitrary rotation and therefore a poor source for readable exponents.
    Reconstruct the dimensional matrix as small fractions and use exact
    Gauss-Jordan elimination instead. None selects the numerical fallback when
    a custom dimensionality cannot be rationalized faithfully.
    """
    if not np.isfinite(matrix).all():
        return None

    rational = [
        [Fraction(float(value)).limit_denominator(max_denominator)
         for value in row]
        for row in matrix
    ]
    reconstructed = np.array(
        [[float(value) for value in row] for row in rational], dtype=float)
    if not np.allclose(reconstructed, matrix, rtol=0.0, atol=1e-12):
        return None

    n_rows, n_cols = matrix.shape
    pivot_cols: list[int] = []
    pivot_row = 0
    for col in range(n_cols):
        row = next((r for r in range(pivot_row, n_rows)
                    if rational[r][col] != 0), None)
        if row is None:
            continue
        rational[pivot_row], rational[row] = rational[row], rational[pivot_row]
        pivot = rational[pivot_row][col]
        rational[pivot_row] = [value / pivot
                               for value in rational[pivot_row]]
        for other in range(n_rows):
            if other == pivot_row or rational[other][col] == 0:
                continue
            factor = rational[other][col]
            rational[other] = [
                value - factor * pivot_value
                for value, pivot_value in zip(
                    rational[other], rational[pivot_row], strict=True)
            ]
        pivot_cols.append(col)
        pivot_row += 1
        if pivot_row == n_rows:
            break

    free_cols = [col for col in range(n_cols) if col not in pivot_cols]
    basis: list[list[Fraction]] = []
    for free_col in free_cols:
        vector = [Fraction(0) for _ in range(n_cols)]
        vector[free_col] = Fraction(1)
        for row, pivot_col in enumerate(pivot_cols):
            vector[pivot_col] = -rational[row][free_col]
        basis.append(vector)
    return basis


def _isolate_target_exact(
        basis: list[list[Fraction]],
        target_index: int,
) -> list[list[Fraction]]:
    """Isolate a target in a rational null-space basis without round-off."""
    containing = [index for index, vector in enumerate(basis)
                  if vector[target_index] != 0]
    if not containing:
        return basis
    pivot_index = max(
        containing, key=lambda index: abs(basis[index][target_index]))
    pivot_scale = basis[pivot_index][target_index]
    pivot = [value / pivot_scale for value in basis[pivot_index]]

    isolated: list[list[Fraction]] = []
    for index, vector in enumerate(basis):
        if index == pivot_index:
            isolated.append(pivot)
            continue
        factor = vector[target_index]
        isolated.append([value - factor * pivot_value
                         for value, pivot_value in zip(vector, pivot, strict=True)])
    return isolated


def _simplify_target_exact(
        basis: list[list[Fraction]],
        target_index: int,
        coefficient_limit: int = 6,
) -> list[list[Fraction]]:
    """Reduce target-group support by exact null-basis combinations.

    Target isolation is not unique: any target-free group may be added to the
    target group without changing its target exponent or dimensional validity.
    A small coordinate-descent search uses that freedom to prefer fewer
    variables and smaller exponents. All arithmetic remains exact.
    """
    target_group = next(
        (index for index, vector in enumerate(basis)
         if vector[target_index] != 0),
        None,
    )
    if target_group is None:
        return basis

    def score(vector: list[Fraction]) -> tuple:
        nonzero = [abs(float(value)) for value in vector if value != 0]
        return (len(nonzero), sum(nonzero), max(nonzero, default=0.0))

    current = basis[target_group]
    changed = True
    while changed:
        changed = False
        best_vector = current
        best_score = score(current)
        for index, feature_group in enumerate(basis):
            if index == target_group:
                continue
            for coefficient in range(-coefficient_limit,
                                     coefficient_limit + 1):
                if coefficient == 0:
                    continue
                candidate = [
                    value + coefficient * feature_value
                    for value, feature_value in zip(
                        current, feature_group, strict=True)
                ]
                candidate_score = score(candidate)
                if candidate_score < best_score:
                    best_vector = candidate
                    best_score = candidate_score
        if best_vector != current:
            current = best_vector
            changed = True

    simplified = [vector[:] for vector in basis]
    simplified[target_group] = current
    return simplified


def _canonical_fraction_vector(
        vector: list[Fraction],
        target_index: int | None = None,
) -> np.ndarray:
    """Convert an exact vector to canonical exponents by uniform scaling."""
    if target_index is not None and vector[target_index] != 0:
        scale = vector[target_index]
        return np.array(
            [float(value / scale) for value in vector], dtype=float)

    denominator_lcm = 1
    for value in vector:
        denominator_lcm = lcm(denominator_lcm, value.denominator)
    integers = [
        value.numerator * (denominator_lcm // value.denominator)
        for value in vector
    ]
    common = 0
    for value in integers:
        common = gcd(common, abs(value))
    if common:
        integers = [value // common for value in integers]
    first = next((value for value in integers if value != 0), 1)
    if first < 0:
        integers = [-value for value in integers]
    return np.array(integers, dtype=float)


@dataclass
class Variable:
    """
    一个工程变量的完整描述。

    这个对象是整个框架的"身份证"——任何进入框架的数据都必须先变成 Variable。
    没有单位的数字、没有角色的列、没有误差范围的测量,都不被接受。

    role 的取值约定:
        - 'state': 状态变量(温度、压力、速度等)
        - 'material_property': 物性参数(密度、粘度、导热系数等)
        - 'geometry': 几何参数(管长、直径、面积等)
        - 'target': 待预测的目标变量(Nu、效率、转化率等)
        - 'derived': 由其他变量计算得到的派生量
    """
    name: str
    unit: str
    role: str
    description: str = ""
    uncertainty: float | None = None  # 相对误差,如 0.05 表示 ±5%

    def __post_init__(self):
        # 立刻验证单位字符串能被 pint 解析——失败比延迟失败好
        try:
            self._quantity = _ureg.parse_expression(self.unit)
        except Exception as e:  # pint 对非法单位抛出多种异常类型,统一折叠为 ValueError
            raise ValueError(f"变量 '{self.name}' 的单位 '{self.unit}' 无法解析: {e}") from e

        if self.role not in {'state', 'material_property', 'geometry', 'target', 'derived'}:
            raise ValueError(f"变量 '{self.name}' 的角色 '{self.role}' 不在允许集合中")

    def dimension_vector(self) -> np.ndarray:
        """
        返回该变量的量纲指数向量,长度 7,对应 SI_BASE_DIMS。
        例:速度 m/s → [0, 1, -1, 0, 0, 0, 0]
        """
        dim_dict = self._quantity.dimensionality
        return np.array([dim_dict.get(d, 0) for d in SI_BASE_DIMS], dtype=float)

    def __repr__(self):
        return f"Variable({self.name}, [{self.unit}], role={self.role})"


class System:
    """
    一组变量的集合。System 是数据进入框架后的第一个组织层。

    典型用法:
        sys = System("heat_exchanger")
        sys.add(Variable("rho", "kg/m**3", "material_property"))
        sys.add(Variable("v", "m/s", "state"))
        ...
        groups = sys.dimensionless_groups()
    """

    def __init__(self, name: str):
        self.name = name
        self.variables: dict[str, Variable] = {}

    def add(self, var: Variable) -> 'System':
        """添加变量。返回 self 以支持链式调用。"""
        if var.name in self.variables:
            raise ValueError(f"变量 '{var.name}' 已存在于系统 '{self.name}' 中")
        self.variables[var.name] = var
        return self

    def add_variable(self, name: str, unit: str, role: str, **kwargs) -> 'System':
        """便利方法,直接传字段而不是先构造 Variable。"""
        return self.add(Variable(name=name, unit=unit, role=role, **kwargs))

    def dimension_matrix(self) -> tuple:
        """
        构建量纲矩阵 D。

        D 的行对应 7 个基本量纲,列对应系统中的变量。
        D[i, j] = 第 j 个变量在第 i 个基本量纲上的指数。

        返回 (D, var_names),var_names 是列对应的变量名顺序。
        """
        var_names = list(self.variables.keys())
        D = np.column_stack([self.variables[n].dimension_vector() for n in var_names])
        return D, var_names

    def dimensionless_groups(self, target: str | None = None) -> list[dict[str, float]]:
        """
        基于 Buckingham π 定理生成无量纲群。

        数学:无量纲群对应于量纲矩阵 D 的零空间(null space)中的整数向量。
        每一个零空间向量 x 都给出一个组合 ∏ var_j ^ x_j,该组合无量纲。

        参数:
            target: 可选的目标变量名。提供时执行"目标隔离":
                - 恰好一个群包含 target,且其指数归一化为 +1
                  (形如 Nu = h·D/k,而不是 h^11·… 这类不可解释的组合);
                - 其余群不含 target(避免"用 target 预测 target")。
                discover() 的默认路径应始终传入 target,保证输出物理可解释。

        返回:每个无量纲群是一个 {var_name: exponent} 字典。
        """
        D, var_names = self.dimension_matrix()

        target_index = None
        if target is not None:
            if target not in var_names:
                raise ValueError(f"目标变量 '{target}' 不在系统 '{self.name}' 中")
            target_index = var_names.index(target)

        # 首选精确有理数零空间。整个向量只做统一缩放，不执行会破坏
        # D @ x = 0 的逐分量取整。
        exact_basis = _rational_null_space(D)
        if exact_basis is not None:
            if not exact_basis:
                return []
            if target_index is not None:
                exact_basis = _isolate_target_exact(
                    exact_basis, target_index)
                exact_basis = _simplify_target_exact(
                    exact_basis, target_index)
            vectors = [
                _canonical_fraction_vector(vector, target_index)
                for vector in exact_basis
            ]
            if all(np.allclose(D @ vector, 0.0, rtol=0.0, atol=1e-9)
                   for vector in vectors):
                return [
                    {var_names[j]: float(vector[j])
                     for j in range(len(var_names))
                     if abs(vector[j]) > 1e-12}
                    for vector in vectors
                ]

        # 无法安全有理化的自定义量纲走数值回退。数值向量仅整体缩放，
        # 返回前仍强制验证无量纲约束。
        ns = null_space(D)
        if ns.size == 0:
            return []
        if target_index is not None:
            ns = self._isolate_target(ns, target_index)

        groups = []
        for i in range(ns.shape[1]):
            vector = ns[:, i].copy()
            if (target_index is not None
                    and abs(vector[target_index]) > 1e-12):
                vector = vector / vector[target_index]
            else:
                scale = np.max(np.abs(vector))
                if scale > 0:
                    vector = vector / scale
                first = next(
                    (value for value in vector if abs(value) > 1e-12), 1.0)
                if first < 0:
                    vector = -vector
            if not np.allclose(D @ vector, 0.0, rtol=0.0, atol=1e-9):
                raise RuntimeError(
                    "数值零空间验证失败:生成的 π 群不满足 D @ x = 0")
            groups.append({
                var_names[j]: float(vector[j])
                for j in range(len(var_names)) if abs(vector[j]) > 1e-12
            })
        return groups

    def _isolate_target(self, ns: np.ndarray, t: int) -> np.ndarray:
        """
        在零空间基上做"目标隔离"的初等列变换:

        1. 选 target 分量绝对值最大的列作为主向量,归一化使 target 指数为 +1;
        2. 从其余列中消去 target 分量(列线性组合,不改变张成的零空间)。

        返回新基——数学上仍是同一零空间的一组基,但只有主向量包含 target。
        """
        ns = ns.copy()
        n_groups = ns.shape[1]
        pivot = int(np.argmax(np.abs(ns[t, :])))
        if abs(ns[t, pivot]) < 1e-12:
            # 没有任何群包含 target(例如 target 变量本身无量纲)
            return ns
        ns[:, pivot] = ns[:, pivot] / ns[t, pivot]  # target 指数 -> 1
        for j in range(n_groups):
            if j == pivot:
                continue
            ns[:, j] = ns[:, j] - ns[t, j] * ns[:, pivot]
        return ns

    def check_data(self, data) -> dict:
        """
        Data Readiness Assessment (Level 1)。

        检查传入的 DataFrame 是否包含 system 中所有变量,并报告缺失情况。
        这是可信度分级的第一关——单位/量纲检查通过、变量齐全,才能升到 Level 2+。

        返回一个报告字典,包含:
            - level: 当前数据可达到的可信度等级(1-5)
            - missing_vars: 缺失的变量列表
            - extra_vars: 数据中存在但未在 schema 中声明的变量
            - n_samples: 数据点数
            - warnings: 文字警告
        """
        report = {
            'level': 0,
            'missing_vars': [],
            'extra_vars': [],
            'n_samples': len(data),
            'warnings': [],
        }

        declared = set(self.variables.keys())
        present = set(data.columns)

        report['missing_vars'] = sorted(declared - present)
        report['extra_vars'] = sorted(present - declared)

        if report['missing_vars']:
            report['warnings'].append(
                f"缺失变量 {report['missing_vars']},无法构建完整量纲矩阵")
            report['level'] = 0
            return report

        # 单位/量纲层面通过,达到 Level 1
        report['level'] = 1

        # 样本数足够吗?这是一个非常粗的启发式
        n_vars = len(declared)
        if report['n_samples'] < 10 * n_vars:
            report['warnings'].append(
                f"样本数 {report['n_samples']} 相对于变量数 {n_vars} 偏少,"
                f"建议至少 {10 * n_vars} 个样本以支持稳健拟合")
        else:
            # 数据量足够,可以构建无量纲群,达到 Level 2
            report['level'] = 2

        return report

    def __repr__(self):
        return f"System('{self.name}', {len(self.variables)} variables)"
