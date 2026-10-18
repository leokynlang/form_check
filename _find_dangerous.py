"""在 S2 数据里找「拟合好但外推危险」的任务——论文关注的那一类。"""
import pandas as pd

df = pd.read_parquet(r'D:\form_check\runs\srsd_full\raw_srsd.parquet')
h1 = df[(df.klass == 'H1') & (df.rel > 0) & (df.n == 50)]

g = h1.groupby('task').agg(
    rmse_med=('test_log_rmse', 'median'),
    far_med=('extrap_far_relerr', 'median'),
    far_p75=('extrap_far_relerr', lambda s: s.quantile(0.75)),
    far_frac_gt10=('extrap_far_relerr', lambda s: (s > 0.10).mean()),
    n_obs=('extrap_far_relerr', 'count'),
).sort_values('far_med', ascending=False)

print('=== H1 任务：拟合质量 vs 外推危险（n=50, rel>0）===')
print(g.to_string(float_format=lambda x: f'{x:.4f}'))
print()

print('=== 论文关注的类别：rmse<0.05（拟合好）且 far_med>10%（危险）===')
sel = g[(g.rmse_med < 0.05) & (g.far_med > 0.10)]
print(sel.to_string(float_format=lambda x: f'{x:.4f}') if len(sel) else '  （无）')
print()

print('=== 放宽：rmse<0.05 且 far_frac_gt10 > 5% ===')
sel2 = g[(g.rmse_med < 0.05) & (g.far_frac_gt10 > 0.05)]
print(sel2.to_string(float_format=lambda x: f'{x:.4f}') if len(sel2) else '  （无）')