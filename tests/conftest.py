"""pytest 配置:保证 `import dimreflex` 可用。

2026-09-05 自 D:\\dimreflex 拆分:本文件夹自包含,vendored 包在 ../dimreflex/。
用 _bootstrap 的元路径查找器强制 vendored 版本(本机若 pip 安装过 dimreflex,
其 finder 会遮蔽本地副本);_bootstrap 不存在时回退可安装包/父目录搜索
(兼容重新组合的扁平仓库)。
"""

import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent          # tests/
_PKG = _HERE.parent / "dimreflex"                # vendored 包
if _PKG.is_dir():
    sys.path.insert(0, str(_HERE.parent))

try:
    import _bootstrap  # noqa: F401 — 强制 vendored dimreflex 优先
    import dimreflex  # noqa: F401
except ImportError:
    try:
        import dimreflex  # noqa: F401 — 已安装(editable)则直接使用
    except ImportError:
        sys.path.insert(0, str(_HERE.parent.parent))
