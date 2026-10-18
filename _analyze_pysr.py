"""分析 PySR 端到端实验结果。"""
import json

import pandas as pd

df = pd.read_csv(r'D:\form_check\runs\pysr_gate\pysr_gate_runs.csv')
print(f'样本数: {len(df)}')
print('任务:', sorted(df.task.unique()))
print()

GATES = [('A', 'p_fc'), ('B', 'p_fitted_f'),
         ('C', 'p_fitted_hc3'), ('D', 'p_exog_hc3')]
ALPHA = 0.05

print('=' * 72)
print('一、逐样本明细')
print('=' * 72)
cols = ['task', 'seed', 'rmse_log', 'r2_log', 'extrap_far_relerr',
        'p_fc', 'p_fitted_f', 'p_fitted_hc3', 'p_exog_hc3']
sub = df[cols].copy()
for _, r in sub.iterrows():
    flags = ''
    for nm, col in GATES:
        v = r[col]
        if pd.notna(v) and v < ALPHA:
            flags += nm
    far = r['extrap_far_relerr']
    far_s = f'{far:.3f}' if pd.notna(far) else '  nan'
    print(f"  {r['task'][:22]:24s} s{r['seed']} "
          f"rmse={r['rmse_log']:.4f} R2={r['r2_log']:+.4f} "
          f"far={far_s:>7s}  拒={flags or '无'}")

print()
print('=' * 72)
print('二、各门卫拒绝率')
print('=' * 72)
for nm, col in GATES:
    v = df[col]
    valid = v.notna()
    rej = (v[valid] < ALPHA).sum()
    print(f'  {nm}: {rej}/{valid.sum()} = {rej/max(valid.sum(),1):.1%}')

print()
print('=' * 72)
print('三、联合沉默（四门卫全不拒）')
print('=' * 72)
any_rej = pd.Series(False, index=df.index)
any_decl = pd.Series(False, index=df.index)
for nm, col in GATES:
    v = df[col]
    any_decl |= v.isna()
    any_rej |= (v.fillna(1.0) < ALPHA)
silent = (~any_rej) & (~any_decl)
print(f'  联合沉默: {silent.sum()}/{len(df)} = {silent.mean():.1%}')
print(f'  有门卫拒绝: {any_rej.sum()}/{len(df)} = {any_rej.mean():.1%}')
print(f'  有门卫无法评估: {any_decl.sum()}')
if silent.any():
    print('  沉默的样本:')
    for _, r in df[silent].iterrows():
        print(f"    {r['task']} s{r['seed']} far={r['extrap_far_relerr']}")

print()
print('=' * 72)
print('四、分任务组：危险 vs 安全')
print('=' * 72)
DANG = ['FeynmanIICh38Eq14', 'FeynmanBonus15']
df['grp'] = df.task.apply(lambda t: 'dangerous' if t in DANG else 'safe')
for grp in ('dangerous', 'safe'):
    g = df[df.grp == grp]
    if not len(g):
        continue
    print(f'\n  [{grp}] n={len(g)}')
    print(f'    rmse 中位: {g.rmse_log.median():.4f}')
    print(f'    R2   中位: {g.r2_log.median():.4f}')
    print(f'    far  中位: {g.extrap_far_relerr.median():.3f}')
    for nm, col in GATES:
        print(f'    {nm} 拒绝率: {(g[col] < ALPHA).mean():.1%}')