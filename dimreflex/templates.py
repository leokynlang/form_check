"""
dimreflex.templates
===================

无量纲群模板。

为什么需要这个:
    Buckingham π 定理告诉我们零空间的维数,但不告诉我们应该选哪组基。
    数学上等价的基,工程意义可能完全不同。Re、Pr、Nu 是几代工程师从无穷
    多种合法组合中筛选出的"好基",这种工程知识必须显式地告诉框架。

设计:
    用户用 Group({"rho": 1, "v": 1, "D": 1, "mu": -1}, name="Re") 这样
    的字典声明一个群,框架验证它确实无量纲,然后把数据投影到这组用户
    指定的群上。
"""


import numpy as np
import pandas as pd

from ._ols import ols as _ols


class Group:
    """一个用户声明的无量纲群,例如 Re = ρvD/μ。"""

    def __init__(self, exponents: dict[str, float], name: str):
        self.exponents = exponents
        self.name = name

    def evaluate(self, data: pd.DataFrame) -> np.ndarray:
        """从数据里计算该群的数值。"""
        out = np.ones(len(data))
        for var, exp in self.exponents.items():
            if var not in data.columns:
                raise KeyError(f"群 '{self.name}' 需要变量 '{var}',但数据中没有")
            out = out * (data[var].values ** exp)
        return out

    def check_dimensionless(self, system) -> bool:
        """验证该群在 system 中确实无量纲。"""
        total = np.zeros(7)
        for var, exp in self.exponents.items():
            if var not in system.variables:
                raise KeyError(f"群 '{self.name}' 引用未声明的变量 '{var}'")
            total += exp * system.variables[var].dimension_vector()
        return np.allclose(total, 0, atol=1e-9)

    def __repr__(self):
        expr = " * ".join(f"{v}^{e:g}" for v, e in self.exponents.items())
        return f"Group({self.name} = {expr})"


def fit_with_templates(data: pd.DataFrame,
                       system,
                       target_group: Group,
                       feature_groups: list[Group]) -> dict:
    """
    在用户声明的模板下拟合幂律。

    例:
        Re = Group({"rho": 1, "v": 1, "D": 1, "mu": -1}, "Re")
        Pr = Group({"cp": 1, "mu": 1, "k": -1}, "Pr")
        Nu = Group({"h": 1, "D": 1, "k": -1}, "Nu")
        fit_with_templates(data, sys, target_group=Nu, feature_groups=[Re, Pr])

    这里框架会拟合 Nu = a * Re^b * Pr^c。
    """
    # 第一步:验证所有群都无量纲——这是用户声明的工程约束的合法性检查
    for g in [target_group, *feature_groups]:
        if not g.check_dimensionless(system):
            return {
                'success': False,
                'reason': f"群 '{g.name}' 不是无量纲的,请检查指数",
            }

    # 共享变量提示:目标群与特征群共享变量(如 Nu 与 Re 共享 D)会在
    # 数学上抬高 R²(π 群伪相关)。不禁止,但必须显式标注。
    target_vars = set(target_group.exponents)
    overlap = {g.name: sorted(target_vars & set(g.exponents))
               for g in feature_groups}
    overlap = {name: vars_ for name, vars_ in overlap.items() if vars_}
    warnings = []
    if overlap:
        warnings.append(
            f"目标群 '{target_group.name}' 与特征群共享变量 {overlap}:"
            f"共享变量会在数学上抬高 R²(π 群伪相关),"
            f"该 R² 不宜与无共享变量的拟合直接比较。")

    # 第二步:计算群的数值
    y = target_group.evaluate(data)
    X = np.column_stack([g.evaluate(data) for g in feature_groups])

    if (not np.isfinite(y).all()) or (not np.isfinite(X).all()):
        return {
            'success': False,
            'reason': '数据中存在非有限值(NaN/inf),无法在 log 空间拟合幂律',
        }

    if (y <= 0).any() or (X <= 0).any():
        return {
            'success': False,
            'reason': '存在非正值,无法在 log 空间拟合幂律',
        }

    # 第三步:log-log 最小二乘(复用共享 OLS 引擎,避免与 fit_power_law
    # / PowerLawFamily 的统计定义漂移;std_errors 口径与二者逐位一致)。
    log_y = np.log(y)
    log_X = np.log(X)
    A = np.column_stack([np.ones(len(log_y)), log_X])
    fit = _ols(A, log_y)
    coef = fit["coef"]
    r2 = fit["r_squared"]
    p = fit["p"]

    a = float(np.exp(coef[0]))
    exponents = {g.name: float(c)
                 for g, c in zip(feature_groups, coef[1:], strict=True)}

    std_errors = dict(zip(['log_a'] + [g.name for g in feature_groups],
                          (float(v) for v in fit["std_errors"]),
                          strict=True))

    parts = [f"{a:.4g}"] + [f"{name}^{exp:.3g}" for name, exp in exponents.items()]
    formula = f"{target_group.name} = " + " * ".join(parts)

    return {
        'success': True,
        'formula': formula,
        'prefactor': a,
        'exponents': exponents,
        'r_squared': r2,
        'std_errors': std_errors,
        'design_rank': int(np.linalg.matrix_rank(A)),
        'design_columns': int(p),
        'identifiable': bool(np.linalg.matrix_rank(A) == p),
        'condition_number': float(np.linalg.cond(A)),
        'n_samples': len(y),
        'warnings': warnings,
    }
