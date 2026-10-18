"""用 H0（严格单项式）任务建校准阈值，再重判 H1 结果。

论文协议：threshold = H0 池 p 值的第 5 百分位（matched 5% budget）。
这里的 H0 池与 H1 同构：同样的 PySR(+,- on log space) 候选生成流程，
只是数据来自真幂律方程。
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, r'D:\form_check')
import srsd_data  # noqa: E402
from run_pysr_gate import OUT, run_one  # noqa: E402

H0_POOL_PATH = os.path.join(OUT, 'h0_pool.csv')


def build_h0_pool(seeds=(0, 1), niterations=25):
    reg = srsd_data.registry_with_meta()
    h0 = [r for r in reg if r['klass'] == 'H0']
    print(f'H0 任务数: {len(h0)}')

    rows = []
    total = len(h0) * len(seeds)
    k = 0
    for r in h0:
        for seed in seeds:
            k += 1
            print(f'[{k}/{total}] {r["key"]} seed={seed} ...', flush=True)
            try:
                row = run_one(r['key'], r['cls'], 50, seed, 0.01, niterations,
                              eq_name=r.get('eq_name', ''))
            except Exception as exc:
                row = {'task': r['key'], 'eq_name': r.get('eq_name', ''),
                       'n': 50, 'seed': seed, 'rel': 0.01,
                       'fail': f'{type(exc).__name__}: {exc}'[:120]}
            rows.append(row)
            print(f'    p=({row.get("p_fc")},{row.get("p_fitted_f")},'
                  f'{row.get("p_fitted_hc3")},{row.get("p_exog_hc3")})',
                  flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(H0_POOL_PATH, index=False)
    print(f'\nH0 池已写入 {H0_POOL_PATH} ({len(df)} 行)')
    return df


def calibrate(h0df, alpha=0.05):
    """每个门卫取 H0 池 p 的第 alpha 分位（与论文口径一致）。"""
    pairs = {
        'A': 'p_fc', 'B': 'p_fitted_f',
        'C': 'p_fitted_hc3', 'D': 'p_exog_hc3',
    }
    tau = {}
    for nm, col in pairs.items():
        v = pd.to_numeric(h0df[col], errors='coerce').dropna()
        tau[nm] = float(np.quantile(v, alpha)) if len(v) else np.nan
        print(f'  tau[{nm}] = {tau[nm]:.4g}   (基于 {len(v)} 个 H0 p 值)')
    return tau


if __name__ == '__main__':
    h0 = build_h0_pool()
    print('\n=== 校准阈值 ===')
    tau = calibrate(h0)

    h1 = pd.read_csv(os.path.join(OUT, 'pysr_gate_runs.csv'))
    pairs = {'A': 'p_fc', 'B': 'p_fitted_f',
             'C': 'p_fitted_hc3', 'D': 'p_exog_hc3'}
    print('\n=== 用校准阈值重判 H1 ===')
    any_rej = pd.Series(False, index=h1.index)
    for nm, col in pairs.items():
        v = pd.to_numeric(h1[col], errors='coerce')
        rej = v < tau[nm]
        any_rej |= rej.fillna(True)   # NaN 计为不可评估 → 不算沉默
        print(f'  {nm}: 拒绝 {int(rej.sum())}/{len(h1)}')
    print(f'\n  **联合沉默（校准后）: '
          f'{int((~any_rej).sum())}/{len(h1)} = {(~any_rej).mean():.1%}**')
    if (~any_rej).any():
        for _, r in h1[~any_rej].iterrows():
            print(f'    {r["task"]} s{r["seed"]} '
                  f'R2={r["r2_log"]:.4f} far={r["extrap_far_relerr"]:.3f}')