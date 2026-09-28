"""
dimreflex.cli
=============

命令行门面:让 dimreflex 能作为独立程序运行(进而可被 PyInstaller 冻结
为 .exe)。它只做:解析参数、读数据、按模式调用库函数、打印报告。

两种模式
--------
1. 量纲感知模式(推荐):CSV 同目录存在 <数据文件名>.schema.json,或用
   --schema <path> 显式指定。schema 声明每列的单位/角色/目标,此时走
   System -> π 群 -> discover 的完整链路(目标隔离、逐特征降级、共享变量
   提示、形式检查、Level 门控全部生效)。
2. 非量纲感知模式(诚实兜底):没有任何 schema 时,仅做原始变量空间的
   log-log 幂律拟合,并在输出顶部明确声明"未校验量纲、不代表物理量纲
   结论"。该模式不再伪装成量纲分析。

退出码:0 成功(至少产生一个成功候选);1 文件/数据/降级错误;2 用法错误。
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__
from .core import System
from .discover import discover, form_check
from .families import discover_multi

_SCHEMA_SUFFIX = ".schema.json"
_ALLOWED_ROLES = {"state", "material_property", "geometry", "target", "derived"}


def _safe_print(text: str = "") -> None:
    """打印一行;控制台编码容不下某个字符时降级替换,绝不抛异常。"""
    try:
        print(text)
    except UnicodeEncodeError:
        enc = getattr(sys.stdout, "encoding", None) or "ascii"
        print(text.encode(enc, errors="replace").decode(enc))


def _strip_quotes(s: str) -> str:
    """剥掉拖拽/粘贴路径带上的引号与首尾空白。"""
    return (s or "").strip().strip('"').strip("'")


def _print_usage() -> None:
    _safe_print("用法: dimreflex <csv> [目标列] [特征列...] [--schema <schema.json>]")
    _safe_print("      不指定列时:最后一列作目标,其余作特征(非量纲模式)")
    _safe_print("      dimreflex               (内置双工况演示)")
    _safe_print("量纲模式: CSV 同目录放 <csv文件名>.schema.json,或传 --schema")


def _demo_data() -> pd.DataFrame:
    """内置演示数据:层流/湍流双工况(与测试 fixture 同一构造,seed=11)。"""
    rng = np.random.default_rng(11)
    Re = np.concatenate([
        rng.uniform(500, 2200, 40),
        rng.uniform(2400, 50000, 40),
    ])
    Nu = np.where(Re < 2300, 0.5 * Re ** (1 / 3), 0.023 * Re ** 0.8)
    Nu *= rng.normal(1.0, 0.01, len(Re))
    return pd.DataFrame({"Re": Re, "Nu": Nu})


def _load_csv(path: str, typed: bool) -> pd.DataFrame | None:
    """读 CSV;失败时打印错误并返回 None。typed=True 表示路径来自交互输入,
    可再给一次机会。"""
    try:
        return pd.read_csv(path)
    except FileNotFoundError:
        if typed:
            _safe_print(f"错误: 找不到 {path}")
        else:
            _safe_print(f"错误: 找不到文件 {path}")
    except Exception as exc:  # noqa: BLE001 — 门面负责给出可读错误
        _safe_print(f"错误: 无法读取 {path}: {exc}")
    return None


def _parse_schema_flag(argv):
    """从参数里抽出 --schema;未知选项视为用法错误。"""
    schema_path = None
    positional = []
    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg == "--schema":
            if i + 1 >= len(argv):
                return None, positional, "--schema 后缺少路径"
            schema_path = argv[i + 1]
            i += 2
        elif arg.startswith("--"):
            return None, positional, f"未知选项: {arg}"
        else:
            positional.append(arg)
            i += 1
    return schema_path, positional, None


def _sidecar_path(csv_path: str) -> Path:
    path = Path(csv_path)
    return path.with_name(path.stem + _SCHEMA_SUFFIX)


def _load_schema(schema_path: str):
    """读取并解析 schema JSON。失败返回 (None, error)。"""
    try:
        with open(schema_path, encoding="utf-8") as fh:
            schema = json.load(fh)
    except Exception as exc:  # noqa: BLE001 — 门面负责给出可读错误
        return None, f"无法读取 schema {schema_path}: {exc}"
    if not isinstance(schema, dict):
        return None, f"schema {schema_path} 顶层必须是 JSON 对象"
    return schema, None


def _build_system(schema: dict, csv_stem: str):
    """从 schema 建 System 并解析 target。失败返回 (None, None)。"""
    variables = schema.get("variables")
    if not isinstance(variables, list) or not variables:
        _safe_print("错误: schema 缺少非空 variables 列表")
        return None, None

    system = System(str(schema.get("system") or csv_stem))
    target = schema.get("target")
    for item in variables:
        if not isinstance(item, dict):
            _safe_print("错误: schema.variables 每项都必须是对象")
            return None, None
        name = item.get("name")
        unit = item.get("unit")
        role = item.get("role")
        if not name or not unit or not role:
            _safe_print(f"错误: schema 变量项缺少 name/unit/role: {item}")
            return None, None
        if role not in _ALLOWED_ROLES:
            _safe_print(f"错误: 变量 '{name}' 的角色 '{role}' 不合法;"
                        f"允许 {sorted(_ALLOWED_ROLES)}")
            return None, None
        try:
            system.add_variable(name, unit, role,
                                uncertainty=item.get("uncertainty"))
        except Exception as exc:  # noqa: BLE001 — pint/重名列等门面错误
            _safe_print(f"错误: 变量 '{name}' 声明无效: {exc}")
            return None, None

    if target is None:
        targets = [item.get("name") for item in variables
                   if item.get("role") == "target"]
        if len(targets) != 1:
            _safe_print("错误: schema 未指定 target,且 role='target' 的变量"
                        f"不是恰好一个({len(targets)} 个): {targets}")
            return None, None
        target = targets[0]
    if target not in system.variables:
        _safe_print(f"错误: target '{target}' 不在 schema 变量声明中")
        return None, None
    return system, target


def _print_report(data: pd.DataFrame, target: str, features,
                  candidates) -> None:
    _safe_print(f"数据: {len(data)} 行 | 目标: {target} | "
                f"特征: {', '.join(features)}")
    _safe_print("-" * 72)
    for c in candidates:
        mark = "✓" if c.success else "✗"
        if c.success:
            detail = f"R²={c.r_squared:.6f}  " + c.formula.replace("\n", " / ")
        else:
            detail = f"(失败) {c.reason}"
        _safe_print(f"  {mark} {c.family:<22} {detail}")


def _print_dimless_report(result: dict, system_name: str,
                          target: str) -> int:
    """打印量纲模式报告;返回退出码(0=有成功候选,1=显式降级/全失败)。"""
    report = result.get("data_report", {})
    _safe_print(f"System: {system_name} | 目标: {target}")
    _safe_print(f"数据等级: Level {report.get('level')}")
    for warning in report.get("warnings", []):
        _safe_print(f"  ⚠ {warning}")
    extra = report.get("extra_vars", [])
    if extra:
        _safe_print(f"  · 未声明列(框架忽略): {extra}")

    if result.get("note"):
        _safe_print(f"  ✗ 显式降级: {result['note']}")
        return 1

    candidates = result.get("candidates", [])
    if not candidates:
        _safe_print("  ✗ 没有可评估的函数族候选")
        return 1
    _safe_print("-" * 72)
    any_success = False
    for c in candidates:
        mark = "✓" if c.success else "✗"
        if c.success:
            any_success = True
            detail = f"R²={c.r_squared:.6f}  " + c.formula.replace("\n", " / ")
        else:
            detail = f"(失败) {c.reason}"
        _safe_print(f"  {mark} {c.family:<22} {detail}")
    return 0 if any_success else 1


def main(argv=None) -> int:
    raw = sys.argv[1:] if argv is None else list(argv)
    raw = [_strip_quotes(a) for a in raw]

    schema_path, argv, usage_error = _parse_schema_flag(raw)
    if usage_error:
        _safe_print(f"错误: {usage_error}")
        _print_usage()
        return 2

    # ---- 双击模式:只有当 stdout 真的是终端时才询问(管道/自动化不问) ----
    typed = False
    if not argv:
        try:
            if sys.stdin.isatty():
                _safe_print("dimreflex — 直接回车运行内置双工况演示,")
                _safe_print("或输入/拖入 CSV 文件路径后回车:")
                line = sys.stdin.readline() or ""
                path = _strip_quotes(line)
                if path:
                    argv = [path]
                    typed = True
        except Exception:  # noqa: BLE001,S110 — 询问失败就当没有交互
            pass

    # ---- 无参数 → 演示 ----
    if not argv:
        data = _demo_data()
        candidates = discover_multi(data, "Nu", ["Re"])
        _safe_print(f"dimreflex {__version__} — 内置双工况演示")
        _print_report(data, "Nu", ["Re"], candidates)
        return 0

    # ---- 有路径:加载数据(交互输错给一次重试机会) ----
    attempts = 0
    while True:
        data = _load_csv(argv[0], typed)
        if data is not None:
            break
        if typed and attempts < 1:
            attempts += 1
            _safe_print("请重新输入 CSV 路径(直接回车退出):")
            line = sys.stdin.readline() or ""
            path = _strip_quotes(line)
            if path:
                argv = [path]
                continue
        _print_usage()
        return 2 if typed else 1

    csv_path = argv[0]

    # ---- 量纲感知模式判定 ----
    if schema_path is None:
        sidecar = _sidecar_path(csv_path)
        if sidecar.exists():
            schema_path = str(sidecar)

    if schema_path is not None:
        if not Path(schema_path).exists():
            _safe_print(f"错误: 找不到 schema 文件 {schema_path}")
            return 1
        schema, err = _load_schema(schema_path)
        if err:
            _safe_print(f"错误: {err}")
            return 1
        system, target = _build_system(schema, Path(csv_path).stem)
        if system is None:
            return 1
        _safe_print(f"dimreflex {__version__} — 量纲感知模式(schema)")
        try:
            result = discover(data, system, target=target, method="multi")
        except Exception as exc:  # noqa: BLE001 — CLI 门面必须给可读错误
            _safe_print(f"错误: 分析失败: {exc}")
            return 1
        return _print_dimless_report(result, system.name, target)

    # ---- 非量纲模式:诚实标注,只做原始变量 log-log 幂律拟合 ----
    _safe_print(f"dimreflex {__version__} — 非量纲感知模式(无 schema)")
    _safe_print(f"  ⚠ 未找到 {Path(csv_path).stem}{_SCHEMA_SUFFIX},"
                "也未指定 --schema;以下为原始变量空间数值拟合,"
                "未校验单位/量纲,不代表物理量纲结论。")

    cols = list(data.columns)
    if len(argv) >= 2:
        target = argv[1]
        features = (argv[2:] if len(argv) >= 3
                    else [c for c in cols if c != target])
    else:
        target, features = cols[-1], cols[:-1]

    if not features:
        _safe_print("错误: 数据只有一列,无法区分目标与特征。")
        _print_usage()
        return 2
    missing = [c for c in [target, *features] if c not in cols]
    if missing:
        _safe_print(f"错误: 数据缺少列 {missing};实际列: {cols}")
        return 1

    try:
        fc = form_check(data, target, features)
        if fc.get("alarm"):
            _safe_print(f"  ⚠ 形式检查: {fc['reason']}")
    except Exception:  # noqa: BLE001,S110 — 形式检查失败不阻断裸模式报告
        pass

    candidates = discover_multi(data, target, features, families=["power_law"])
    _print_report(data, target, features, candidates)
    return 0 if any(c.success for c in candidates) else 1


if __name__ == "__main__":
    sys.exit(main())
