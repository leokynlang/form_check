"""查 S2 已有数据：H1 任务的 R2(log) 分布，确认幂律近似的有效性。"""
import os
import pandas as pd

BASE = r'D:\form_check\runs'
p = os.path.join(BASE, 'srsd_full', 'raw_srsd.parquet')
df = pd.read_parquet(p)
print('列:', list(df.columns))
print('形状:', df.shape)
print()

# 找 R2 相关列
r2cols = [c for c in df.columns if 'r2' in c.lower() or 'r_sq' in c.lower()]
print('R2 相关列:', r2cols)
print()

if r2cols:
    r2 = r2cols[0]
    h1 = df[df.klass == 'H1']
    print(f'=== H1 任务 {r2} 分布 ===')
    print(h1.groupby('n')[r2].describe()[['count', 'mean', 'min', '50%', 'max']])
    print()
    print('=== 有多少 H1 拟合 R2 >= 0.99（论文关注的门槛）===')
    for n in sorted(h1.n.unique()):
        sub = h1[h1.n == n]
        frac = (sub[r2] >= 0.99).mean()
        print(f'  n={n}: {frac:.1%}  ({int((sub[r2] >= 0.99).sum())}/{len(sub)})')
    print()
    print('=== 逐任务 R2 中位数（n=50）===')
    t = h1[h1.n == 50].groupby('task')[r2].median().sort_values()
    print(t.to_string())
