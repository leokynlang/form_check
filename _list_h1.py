import sys

sys.path.insert(0, r'D:\form_check')
import srsd_data

reg = srsd_data.registry_with_meta()
h1 = [r for r in reg if r['klass'] == 'H1']
print('H1 总数:', len(h1))
print()
for r in h1:
    print(f"  key={r['key']:26s} eq_name={r['eq_name']:18s} "
          f"vars={r['n_vars']}")