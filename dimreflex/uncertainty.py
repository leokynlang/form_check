"""
dimreflex.uncertainty
=====================

测量不确定度的蒙特卡洛传播。

为什么需要这个:
    Variable.uncertainty(相对误差)自 v0.1 起声明,但从未接入计算。
    std_errors 只回答"回归残差噪声有多大"(条件于测量无误差),回答不了
    "仪器精度会如何动摇指数"。两者是独立的噪声来源,必须分开报告:

    - std_errors(回归):拟合残差折算的系数标准误,n 越大越小
    - MC 传播(测量):按声明的相对误差扰动数据重新拟合,n 增大不消失
      (它反映仪器极限,不是样本量问题)

方法:
    1. 收集所有声明了 uncertainty 的变量(含 target),none 声明则显式拒绝
       ——框架不猜
    2. 固定 π 群定义(群定义由量纲决定,不随数值扰动变化),逐次模拟:
       对每个不确定变量乘 exp(N(0, u)) ——log 空间加性高斯,与框架的
       log-log 拟合假设自洽,且中位数无偏
    3. 重新投影 π 数值、重拟合幂律,收集指数分布
    4. 报告均值/标准差/95% 区间,并附回归 std_errors 作对照

诚实的附带效应:
    扰动特征变量会引入 EIV(误差-in-变量)衰减——MC 均值可能略低于
    无噪拟合值。这不是 bug:它正是测量噪声对指数的总体影响,
    包括 std_errors 原理上看不到的那部分。
"""


import numpy as np
import pandas as pd

from .discover import compute_groups, fit_power_law


def propagate_uncertainty(data: pd.DataFrame,
                          system,
                          target: str,
                          n_sims: int = 200,
                          random_state: int = 42,
                          progress: dict | None = None) -> dict:
    """按 Variable.uncertainty 声明做蒙特卡洛传播,给出指数的测量噪声区间。

    返回(失败也是结构化结果,绝不抛异常):
        {
            'success': bool,
            'reason': str,                    # 失败原因
            'variables_perturbed': {名: u},   # 实际扰动了谁、幅度多少
            'fit': 无扰动基准拟合(与 fit_power_law 同构),
            'exponents_mc': {群名: {'mean','std','q025','q975','summary'}},
            'log_a_mc': {...},
            'regression_std_errors': {...},   # 对照:回归残差噪声,勿与本结果混淆
            'note': 两种噪声的语义区分,
        }
    """
    out = {'success': False, 'reason': '', 'variables_perturbed': {},
           'fit': None, 'exponents_mc': {}, 'log_a_mc': {},
           'regression_std_errors': {}, 'n_sims': int(n_sims),
           'n_failed_sims': 0,
           'note': ('regression_std_errors 是回归残差噪声(条件于测量无误差);'
                    '本结果是测量噪声传播(条件于拟合形式正确)。两者独立,'
                    '不应相加或直接比较。')}

    # -- 收集不确定度声明:框架不猜,没声明就显式拒绝 -----------------------
    declared = {}
    for name, var in system.variables.items():
        if var.uncertainty is None:
            continue
        if name not in data.columns:
            out['reason'] = f"变量 '{name}' 声明了 uncertainty 但数据缺列"
            return out
        u = float(var.uncertainty)
        if not np.isfinite(u) or u < 0:
            out['reason'] = f"变量 '{name}' 的 uncertainty={u} 非法(需非负有限)"
            return out
        if u > 0:
            declared[name] = u
    if not declared:
        out['reason'] = ("没有任何变量声明 uncertainty(如 "
                         "sys.add_variable('v', 'm/s', 'state', uncertainty=0.02))"
                         "——框架不会在无测量噪声信息时虚构传播结果。")
        return out
    out['variables_perturbed'] = dict(declared)

    # -- 数据门卫:缺列等事故在 compute_groups 之前显式拒绝(绝不抛异常) ----
    report = system.check_data(data)
    if report['level'] < 2:
        out['reason'] = (f"数据等级 Level {report['level']},不足以进行"
                         f"不确定度传播: {report.get('missing', '')}")
        return out

    # -- π 群定义只算一次(由量纲决定,与数值无关) ---------------------------
    groups = system.dimensionless_groups(target=target)
    target_groups = [i for i, g in enumerate(groups) if target in g]
    if not groups or not target_groups or len(groups) < 2:
        out['reason'] = "无量纲群不足(需恰好隔离出目标群且另有特征群)"
        return out
    y_idx = target_groups[0]
    y_col = None
    for i, g in enumerate(groups):
        if i == y_idx:
            y_col = f"PI_{i}_" + "_".join(f"{v}^{e:g}" for v, e in g.items())
    x_cols = [f"PI_{i}_" + "_".join(f"{v}^{e:g}" for v, e in g.items())
              for i, g in enumerate(groups) if i != y_idx]

    base_fit = fit_power_law(compute_groups(data, groups), y_col, x_cols)
    if not base_fit.get('success'):
        out['reason'] = f"基准拟合失败: {base_fit.get('reason')}"
        return out
    out['fit'] = base_fit
    out['regression_std_errors'] = base_fit.get('std_errors', {})

    # -- 蒙特卡洛主循环 ----------------------------------------------------
    rng = np.random.default_rng(random_state)
    cols = {name: data[name].to_numpy(dtype=float) for name in declared}
    exp_samples, loga_samples = [], []
    n_failed = 0
    for _ in range(int(n_sims)):
        perturbed = data.copy()
        for name, u in declared.items():
            perturbed[name] = cols[name] * np.exp(rng.normal(0.0, u,
                                                             len(data)))
        pi = compute_groups(perturbed, groups)
        y = pi[y_col].to_numpy()
        ok = (np.isfinite(y).all() and (y > 0).all()
              and all(np.isfinite(pi[c]).all() and (pi[c] > 0).all()
                      for c in x_cols))
        if not ok:
            n_failed += 1          # 扰动把某群打进负值/溢出:该次作废留痕
            continue
        fit = fit_power_law(pi, y_col, x_cols)
        if not fit.get('success'):
            n_failed += 1
            continue
        exp_samples.append(list(fit['exponents'].values()))
        loga_samples.append(np.log(fit['prefactor']))

    out['n_failed_sims'] = n_failed
    if len(exp_samples) < max(30, n_sims // 10):
        out['reason'] = (f"有效模拟次数不足({len(exp_samples)}/{n_sims}):"
                         f"扰动导致过多次拟合失败,不确定度声明可能过大")
        return out
    if n_failed > n_sims * 0.1:
        out['note'] += (f" 警告: {n_failed}/{n_sims} 次模拟因扰动溢出被废弃,"
                        f"区间可能低估。")

    # -- 汇总:均值/标准差/95% 区间 + 人读摘要 ------------------------------
    arr = np.array(exp_samples)
    names = list(base_fit['exponents'])
    for j, name in enumerate(names):
        s = arr[:, j]
        out['exponents_mc'][name] = {
            'mean': float(s.mean()), 'std': float(s.std(ddof=1)),
            'q025': float(np.quantile(s, 0.025)),
            'q975': float(np.quantile(s, 0.975)),
            'summary': f"{name} = {s.mean():.3g} ± {s.std(ddof=1):.2g}",
        }
    la = np.array(loga_samples)
    out['log_a_mc'] = {
        'mean': float(la.mean()), 'std': float(la.std(ddof=1)),
        'q025': float(np.quantile(la, 0.025)),
        'q975': float(np.quantile(la, 0.975)),
    }
    out['success'] = True
    return out
