"""S2 数据：H1 的 test_log_rmse 与外推误差分布。"""
import os

import pandas as pd

df = pd.read_parquet(r'D:\form_check\runs\srsd_full\raw_srsd.parquet')
h1 = df[(df.klass == 'H1') & (df.rel > 0)]

print('=== H1 test_log_rmse 分布 ===')
print(h1.groupby('n')['test_log_rmse'].describe()[
    ['count', 'mean', '50%', 'max']].to_string())
print()

print('=== 逐任务中位 test_log_rmse（n=50, rel=1%）===')
sub = h1[(h1.n == 50) & (h1.rel == 0.01)]
t = sub.groupby('task')['test_log_rmse'].median().sort_values()
for k, v in t.items():
    print(f'  {k:22s} {v:.4f}')
print()

print('=== 外推误差（far）分布，n=50 rel=1% ===')
print(sub['extrap_far_relerr'].describe()[['count', 'mean',
                                           '50%', 'max']].to_string())
print()
print('  far > 10% 占比:',
      f"{(sub['extrap_far_relerr'] > 0.10).mean():.1%}")
print()
print('=== 该层内 test_log_rmse 很小的比例（拟合好的）===')
for thr in (0.01, 0.02, 0.05):
    print(f'  rmse < {thr}: {(sub["test_log_rmse"] < thr).mean():.1%}')