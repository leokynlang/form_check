"""强制本目录 vendored 的 dimreflex 包优先于任何 pip 安装版本。

2026-09-05 自 D:\\dimreflex 拆分(见 context.json)。本文件夹自包含:
dimreflex 包位于 ./dimreflex/。若机器上恰好 pip 安装过 dimreflex
(如 editable 指向旧仓库),其 meta_path finder 会遮蔽本目录副本,
故用元路径查找器显式拦截 dimreflex* 的导入。
"""
import importlib.abc
import importlib.util
import sys
from pathlib import Path

_PKG = Path(__file__).resolve().parent / "dimreflex"


def _spec_for(name):
    if name == "dimreflex":
        return importlib.util.spec_from_file_location(
            name, str(_PKG / "__init__.py"),
            submodule_search_locations=[str(_PKG)])
    rest = name[len("dimreflex."):]
    mod = _PKG / (rest.replace(".", "/") + ".py")
    if mod.is_file():
        return importlib.util.spec_from_file_location(name, str(mod))
    sub_init = _PKG / rest.replace(".", "/") / "__init__.py"
    if sub_init.is_file():
        return importlib.util.spec_from_file_location(
            name, str(sub_init), submodule_search_locations=[str(sub_init.parent)])
    return None


class _VendoredFinder(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "dimreflex" or fullname.startswith("dimreflex."):
            return _spec_for(fullname)
        return None


if _PKG.is_dir():
    sys.meta_path.insert(0, _VendoredFinder())
