"""查关键任务的真实 sympy 公式。"""
import sys

import sympy

sys.path.insert(0, r'D:\form_check')
import srsd_data  # noqa: E402

reg = {r['key']: r for r in srsd_data.registry_with_meta()}

for key in ['FeynmanBonus15', 'FeynmanIICh38Eq14', 'FeynmanIIICh4Eq32',
            'FeynmanIIICh4Eq33', 'FeynmanIICh11Eq28', 'FeynmanICh10Eq7']:
    r = reg.get(key)
    if r is None:
        print(f'!! 缺 {key}')
        continue
    inst = r['cls']()
    f = inst.sympy_eq
    print(f'{key:22s} ({r["eq_name"]:18s})')
    print(f'    sympy: {f}')
    try:
        print(f'    latex: {sympy.latex(f)}')
    except Exception:
        pass
    print(f'    变量数: {r["n_vars"]}')
    print()