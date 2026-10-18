"""
dimreflex._ols
==============

框架内唯一的 log 空间最小二乘(OLS)引擎。

历史背景(为何需要本模块):fit_power_law、fit_with_templates、
PowerLawFamily、anomaly 的基准拟合与 scenarios 的内层拟合曾各自维护
一份几乎相同的 lstsq + R² + 标准误代码。数值逻辑每复制一份,统计补丁
(VIF、岭、自由度修正等)就要同步改多处,并用对等性测试人工锁定漂移。
本模块把这些计算收拢到一处;各调用方只负责自己的 API 形态(公式字符串、
候选模型、键名),不再重复统计数学。

设计说明:
- design 的调用方负责提供"第一列为 1(截距)"的完整设计矩阵;
- response 是 log 空间响应(拟合对象为幂律时,即 log(target));
- ols() 一次返回系数、样本内 R²、标准误、条件数等全部诊断,
  调用方按需取用(计算开销对工程数据规模可忽略);
- 所有数值口径与迁移前逐位一致(见 tests/test_discover.py::TestStdErrorsParity)。
"""


import numpy as np


def ols(design: np.ndarray, response: np.ndarray) -> dict:
    """普通最小二乘 + 完整诊断(单一实现,所有 log 空间拟合共用)。

    参数
    ----
    design : (n, p) 数组,调用方应含截距列([1, log x ...])。
    response : (n,) log 空间响应。

    返回
    ----
    dict: coef / prediction_log / ss_res / ss_tot / r_squared /
          rank / condition_number / std_errors(截距+特征顺序)/
          sigma2 / n / p。rank < p 时 std_errors 仍按伪逆给出(与
          迁移前行为一致,identifiable 判定由调用方负责)。
    """
    design = np.asarray(design, dtype=float)
    response = np.asarray(response, dtype=float)
    n, p = design.shape
    if len(response) != n:
        raise ValueError("design 行数必须等于 response 长度")

    coef, _, rank, _ = np.linalg.lstsq(design, response, rcond=None)
    prediction = design @ coef
    ss_res = float(np.sum((response - prediction) ** 2))
    ss_tot = float(np.sum((response - response.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    # 系数标准误(假设残差正态):n > p 才有正的自由度,否则显式 NaN
    if n > p:
        sigma2 = ss_res / (n - p)
        cov = sigma2 * np.linalg.pinv(design.T @ design)
        std_errors = np.sqrt(np.diag(cov))
    else:
        sigma2 = float("nan")
        std_errors = np.full(p, np.nan)

    return {
        "coef": coef,
        "prediction_log": prediction,
        "prediction": prediction,
        "ss_res": ss_res,
        "ss_tot": ss_tot,
        "r_squared": float(r2),
        "rank": int(rank),
        "condition_number": float(np.linalg.cond(design)),
        "std_errors": std_errors,
        "sigma2": float(sigma2),
        "n": int(n),
        "p": int(p),
    }


def ols_r2(design: np.ndarray, response: np.ndarray,
           coef: np.ndarray) -> float:
    """对任意设计矩阵与给定系数(OLS/稳健/岭均可)计算样本内 R²。"""
    prediction = np.asarray(design, dtype=float) @ np.asarray(coef, dtype=float)
    response = np.asarray(response, dtype=float)
    ss_res = float(np.sum((response - prediction) ** 2))
    ss_tot = float(np.sum((response - response.mean()) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
