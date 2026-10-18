"""S2 H1 单样本：拟合质量与外推危险的联合分布。"""
import pandas as pd

df = pd.read_parquet(r'D:\form_check\runs\srsd_full\raw_srsd.parquet')
h1 = df[(df.klass == 'H1') & (df.rel > 0) & (df.n == 50)]
h1 = h1[h1.extrap_far_relerr.notna()]

print(f'样本数: {len(h1)}')
print()
print('=== rmse 分位 ===')
for q in (0.1, 0.25, 0.5, 0.75, 0.9, 0.99):
    print(f'  q{q:<5} rmse = {h1.test_log_rmse.quantile(q):.4f}')
print()

print('=== 联合：rmse 门槛 × far>10% 比例 ===')
print(f'{"rmse<thr":>10} {"样本数":>8} {"其中far>10%":>12} {"占比":>8}')
for thr in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5):
    sub = h1[h1.test_log_rmse < thr]
    if len(sub):
        k = int((sub.extrap_far_relerr > 0.10).sum())
        print(f'{thr:>10.2f} {len(sub):>8d} {k:>12d} {k/len(sub):>7.1%}')

print()
print('=== 反查：far>10% 的样本，其 rmse 分布 ===')
bad = h1[h1.extrap_far_relerr > 0.10]
print(f'  far>10% 样本数: {len(bad)} / {len(h1)} ({len(bad)/len(h1):.1%})')
if len(bad):
    print(bad.test_log_rmse.describe()[['min', '25%', '50%', '75%', 'max']]
          .to_string(float_format=lambda x: f'{x:.4f}'))
    print()
    print('  其中 rmse<0.1 的（过筛后仍危险）:',
          int((bad.test_log_rmse < 0.1).sum()))
    t = bad[bad.test_log_rmse < 0.1].groupby('task').size()
    print(t.to_string() if len(t) else '  （无）')