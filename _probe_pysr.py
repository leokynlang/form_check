"""最小探针：PySR 能否在本机跑通，输出是什么形式。"""
import time

import numpy as np

print('导入 PySR ...')
t0 = time.time()
from pysr import PySRRegressor
print(f'  导入耗时 {time.time() - t0:.1f}s')

# 造一个明确不是幂律的数据：y = 2*x1^2 + 3*x1（二次+线性）
rng = np.random.default_rng(0)
n = 80
x1 = rng.uniform(1.0, 5.0, n)
x2 = rng.uniform(1.0, 5.0, n)
y_true = 2.0 * x1 ** 2 + 3.0 * x1 + 0.5 * x1 * x2
y = y_true * np.exp(rng.normal(0, 0.01, n))

X = np.column_stack([x1, x2])
print(f'\n数据: n={n}, X.shape={X.shape}')

print('\n启动 PySR（首次可能编译 Julia，请稍候）...')
t0 = time.time()
model = PySRRegressor(
    niterations=10,
    populations=8,
    binary_operators=['+', '-', '*', '/'],
    unary_operators=['square'],
    maxsize=12,
    deterministic=True,
    parallelism='serial',
    random_state=0,
    temp_equation_file=False,
    verbosity=0,
    progress=False,
)
try:
    model.fit(X, y)
    print(f'  拟合完成，耗时 {time.time() - t0:.1f}s')
    print('\n=== 最佳表达式 ===')
    print('  sympy :', model.sympy())
    print('  latex :', model.latex())
    print('  R2    :', model.score(X, y))
    print('\n=== 候选表 ===')
    print(model.equations_[['complexity', 'loss', 'equation']].to_string())
except Exception as e:
    print(f'!! PySR 失败: {type(e).__name__}: {e}')