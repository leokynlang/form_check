"""在「幂律拟合良好」的 H1 任务上测 PySR 能否达标。"""
import sys

import numpy as np

sys.path.insert(0, r'D:\form_check')
import srsd_data  # noqa: E402
from pysr import PySRRegressor  # noqa: E402

# 取自 S2 实测：这些任务的中位 test_log_rmse ≈ 0.010–0.014
GOOD = ['FeynmanICh10Eq7', 'FeynmanIICh11Eq28', 'FeynmanIIICh4Eq32',
        'FeynmanIICh11Eq27', 'FeynmanICh41Eq16']

reg = srsd_data.registry_with_meta()
by_key = {r['key']: r for r in reg}

for key in GOOD:
    cand = [k for k in by_key if key.lower() in k.lower()]
    if not cand:
        cand = [k for k in by_key if key.replace('Feynman', '')[:8].lower()
                in k.lower()]
    if not cand:
        print(f'!! 未找到 {key}，可用示例: {list(by_key)[:5]}')
        continue
    entry = by_key[cand[0]]
    if entry['klass'] != 'H1':
        print(f'  跳过 {entry["key"]}（klass={entry["klass"]}）')
        continue

    fn = srsd_data.make_dataset_fn(entry['cls'])
    df, feats = fn(0, 50, 0.01)
    names = feats.split(',')
    X = df[names].to_numpy(float)
    y = df['target'].to_numpy(float)
    lx, ly = np.log(X), np.log(y)

    print(f'\n=== {entry["key"]}  (vars={len(names)}, '
          f'log(y) std={ly.std():.3f}) ===')
    for label, ops in [('+,-', ['+', '-']), ('+,-,*', ['+', '-', '*'])]:
        m = PySRRegressor(
            niterations=25, populations=12,
            binary_operators=ops, unary_operators=[],
            maxsize=14, deterministic=True, parallelism='serial',
            random_state=0, temp_equation_file=False,
            verbosity=0, progress=False,
        )
        m.fit(lx, ly)
        r2 = float(m.score(lx, ly))
        rmse = float(np.sqrt(np.mean((ly - m.predict(lx)) ** 2)))
        print(f'  [{label:6s}] rmse={rmse:.4f}  R2={r2:.4f}  '
              f'{str(m.sympy())[:70]}')