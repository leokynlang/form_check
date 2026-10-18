"""诊断：PySR 在 log 空间、不同算子集下的输出形态。"""
import numpy as np
import sys

sys.path.insert(0, r'D:\form_check')
import srsd_data  # noqa: E402
from pysr import PySRRegressor  # noqa: E402

reg = srsd_data.registry_with_meta()
h1 = [r for r in reg if r['klass'] == 'H1']
entry = next(r for r in h1 if 'ICh27Eq6' in r['key'])
print('任务:', entry['key'], '| 变量数:', entry['n_vars'])

fn = srsd_data.make_dataset_fn(entry['cls'])
df, feats = fn(0, 50, 0.01)
feature_names = feats.split(',')
print('特征:', feature_names)
print('数据形状:', df.shape)
print(df.head(3).to_string())

X = df[feature_names].to_numpy(float)
y = df['target'].to_numpy(float)
lx, ly = np.log(X), np.log(y)
print(f'\nlog(y) 范围: [{ly.min():.2f}, {ly.max():.2f}]  '
      f'std={ly.std():.3f}')

for label, ops in [
    ('只有乘除 (*,/)', ['*', '/']),
    ('只有加减 (+,-)', ['+', '-']),
    ('加减乘 (+,-,*)', ['+', '-', '*']),
    ('完整 (+,-,*,/)', ['+', '-', '*', '/']),
]:
    print(f'\n{"=" * 60}')
    print(f'算子集: {label}')
    print('=' * 60)
    m = PySRRegressor(
        niterations=15, populations=10,
        binary_operators=ops, unary_operators=[],
        maxsize=12, deterministic=True, parallelism='serial',
        random_state=0, temp_equation_file=False,
        verbosity=0, progress=False,
    )
    try:
        m.fit(lx, ly)
        print('  最佳:', m.sympy())
        print('  R2(log):', round(float(m.score(lx, ly)), 5))
        eq = m.equations_
        if eq is not None and len(eq):
            top = eq.nlargest(3, 'complexity')[['complexity', 'loss', 'equation']]
            print(top.to_string(index=False))
    except Exception as e:
        print('  失败:', type(e).__name__, str(e)[:100])