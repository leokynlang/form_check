"""检查可用依赖 + 现有数据结构。"""
import importlib
import os

print('=' * 60)
print('SR / ML 库')
print('=' * 60)
for m in ('pysr', 'pysindy', 'gplearn', 'sympy', 'sklearn',
          'numpy', 'scipy', 'pandas', 'pyarrow'):
    try:
        mod = importlib.import_module(m)
        ver = getattr(mod, '__version__', '?')
        print(f'  {m:12s} OK   {ver}')
    except Exception as e:
        print(f'  {m:12s} --   {type(e).__name__}')

print()
print('=' * 60)
print('现有数据文件')
print('=' * 60)
paths = [
    'runs/full/raw.parquet',
    'runs/srsd_full/raw_srsd.parquet',
    'runs/cd_svd/s2_cd_svd.csv',
    'runs/cd_svd',
]
for p in paths:
    full = os.path.join(r'D:\form_check', p)
    if os.path.isdir(full):
        n = len([f for f in os.listdir(full)])
        print(f'  [DIR ] {p}  ({n} files)')
    elif os.path.exists(full):
        print(f'  [FILE] {p}  {os.path.getsize(full)/1024:.0f} KB')
    else:
        print(f'  [MISS] {p}')

print()
print('=' * 60)
print('SRSD 数据模块接口')
print('=' * 60)
import sys
sys.path.insert(0, r'D:\form_check')
try:
    import srsd_data
    fns = [n for n in dir(srsd_data) if not n.startswith('_')]
    print('  srsd_data 导出:', fns[:20])
except Exception as e:
    print('  导入失败:', e)