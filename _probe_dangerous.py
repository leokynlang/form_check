"""在「拟合好但外推危险」的任务上测 PySR + 门卫。"""
import json
import sys

import numpy as np

sys.path.insert(0, r'D:\form_check')
import srsd_data  # noqa: E402
from common import form_check_with_p, reset_test  # noqa: E402
from run_pysr_gate import (far_probe, parse_linear_powerlaw,  # noqa: E402
                           predict_powerlaw)
from pysr import PySRRegressor  # noqa: E402

DANGEROUS = ['FeynmanBonus15', 'FeynmanIICh38Eq14']
reg = {r['key']: r for r in srsd_data.registry_with_meta()}

for key in DANGEROUS:
    r = reg.get(key)
    if r is None:
        print(f'!! 缺 {key}')
        continue
    print(f'\n{"=" * 66}')
    print(f'{key}  ({r["eq_name"]})  vars={r["n_vars"]}')
    print('=' * 66)
    for seed in (0, 1):
        fn = srsd_data.make_dataset_fn(r['cls'])
        df, feats = fn(seed, 50, 0.01)
        names = feats.split(',')
        X = df[names].to_numpy(float)
        y = df['target'].to_numpy(float)
        lx, ly = np.log(X), np.log(y)

        m = PySRRegressor(
            niterations=25, populations=12,
            binary_operators=['+', '-'], unary_operators=[],
            maxsize=14, deterministic=True, parallelism='serial',
            random_state=seed, temp_equation_file=False,
            verbosity=0, progress=False,
        )
        m.fit(lx, ly)
        expr = m.sympy()
        pred = m.predict(lx)
        rmse = float(np.sqrt(np.mean((ly - pred) ** 2)))
        ss_tot = float(np.sum((ly - ly.mean()) ** 2))
        r2 = 1 - float(np.sum((ly - pred) ** 2)) / ss_tot if ss_tot else np.nan

        parsed = parse_linear_powerlaw(expr, names)
        if parsed is None:
            print(f'  seed={seed}: 无法解析为幂律: {expr}')
            continue
        pref, expo = parsed
        far = far_probe.__wrapped__(r['cls'], X, names, pref, expo) \
            if hasattr(far_probe, '__wrapped__') else far_probe(
                r['cls'], X, names, pref, expo)

        fc = form_check_with_p(df, 'target', names, seed)
        b = reset_test(lx, ly, variant='fitted', vcov='classical')
        c = reset_test(lx, ly, variant='fitted', vcov='hc3')
        d = reset_test(lx, ly, variant='exog', vcov='hc3')

        print(f'  seed={seed}: rmse={rmse:.4f} R2={r2:.5f} far={far:.3f}')
        print(f'    expr: {expr}')
        print(f'    p: fc={fc.get("p_fc")} f={b.get("p")} '
              f'C={c.get("p")} D={d.get("p")}')
        rej = []
        for nm, p in [('fc', fc.get('p_fc')), ('f', b.get('p')),
                      ('C', c.get('p')), ('D', d.get('p'))]:
            if p is not None and not np.isnan(p) and p < 0.05:
                rej.append(nm)
        print(f'    拒绝的门卫: {rej if rej else "无（联合沉默）"}')