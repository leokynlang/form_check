# 第四部分 门面层与项目配置

> 本部分覆盖 `cli.py`(326 行)、`launcher.py`(22 行)、`__init__.py`(89 行)、`pyproject.toml`(66 行)。门面层不新增任何数学:它把库函数装配成"双击可用的 exe"和"pip install 后的命令",并且同样遵守"诚实兜底、显式降级"哲学。

---

## 4.1 `cli.py` — 命令行门面

**文件定位**:docstring:"它只做:解析参数、读数据、按模式调用库函数、打印报告。" 两种模式:量纲感知(schema 存在时走完整链路)与非量纲感知(无 schema 时的诚实兜底)。退出码契约:0=至少一个成功候选;1=文件/数据/降级错误;2=用法错误。

### 4.1.1 辅助函数(cli.py:37-119)

| 函数 | 计算 | 原因 | 导向 |
|---|---|---|---|
| `_safe_print(text)` | 打印时捕获 `UnicodeEncodeError`,降级为 replace 编码 | Windows GBK 控制台容不下 R²/π 等字符(CONTEXT.md 记录的教训),打印失败不能让报告死掉 | 所有输出路径 |
| `_strip_quotes(s)` | 剥掉拖拽路径带上的引号与空白 | Windows"拖拽 CSV 到 exe"会把路径带引号传入 | argv 预处理 |
| `_print_usage()` | 打印用法 | — | exit 2 路径 |
| `_demo_data()` | seed=11 生成层流/湍流双工况:`Re∈[500,2200]∪[2400,50000]`,`Nu = 0.5·Re^{1/3}(Re<2300) 或 0.023·Re^{0.8}`,乘 1% 噪声 | 双击 exe 无参数时的内置演示;双工况构造与测试 fixture 同源,演示即分段幂律族的活广告 | `main` 的无参数分支 |
| `_load_csv(path, typed)` | pd.read_csv,文件不存在/解析失败打印可读错误返回 None | 门面负责可读错误而非 traceback | main 的数据加载(交互输入可重试一次) |
| `_parse_schema_flag(argv)` | 线性扫描抽出 `--schema <path>`,未知 `--` 选项报用法错误 | 避免引入 argparse 的重量级行为(保持拖拽即用的语义) | main 的参数分离 |
| `_sidecar_path(csv_path)` | `<stem>.schema.json` 同目录副文件名 | "CSV 旁边放同名 schema 即自动进入量纲模式"的约定 | main 的模式判定 |
| `_load_schema(schema_path)` | json.load + 顶层必须是对象 | — | `_build_system` 输入 |
| `_build_system(schema, csv_stem)` | 遍历 schema.variables 构造 System(name/unit/role 必填、role 白名单、uncertainty 透传);target 缺省时取唯一 role='target' 的变量;target 不在声明中报错 | schema 是"自然语言之外的最薄声明层"——量纲模式的一切前提 | main 量纲模式分支 |

### 4.1.2 `_print_report` 与 `_print_dimless_report`(cli.py:166-210)

- **【计算】** 候选列表打印:成功者 `✓ family R²=... formula`,失败者 `✗ (失败) reason`;量纲模式额外打印 System 名、数据等级、warnings、extra_vars、显式降级 note;返回退出码(有成功候选 0,否则 1)。
- **【原因】** 失败候选与降级 note 必须出现在 stdout——用户能看到"为什么没有结果"与"框架的警告",这是显式降级哲学的展示面。
- **【导向】** main 的返回值 → PyInstaller exe 的进程退出码 → 拖拽/脚本调用方据此判断。

### 4.1.3 `main(argv=None)`(cli.py:213-322)

- **【计算】** 分支树:
  1. 参数解析与引号剥离;
  2. **双击模式**:无参数且 stdin 是 TTY → 提示输入路径(管道/自动化不问,不阻塞脚本);
  3. 仍无参数 → 内置双工况演示 `discover_multi(Nu~Re)`;
  4. 有路径 → 读 CSV(交互输错给一次重试);
  5. 量纲模式判定:`--schema` 显式 > 同目录副文件自动 > 无;
  6. 量纲模式:`discover(method='multi')` → 打印(若崩溃打印可读错误,exit 1);
  7. 非量纲模式:顶部打印"未校验单位/量纲、不代表物理量纲结论"声明 → 先跑 `form_check` 报警形式失配 → `discover_multi(families=['power_law'])` → 打印。
- **【原因】** 两模式的设计原因(docstring):"该模式不再伪装成量纲分析"——无 schema 时的诚实兜底是框架哲学在产品层的体现(宁可声明自己只是数值拟合,也不冒充物理结论)。非量纲模式也跑 form_check 是因为形式失配与量纲无关,纯数值拟合同样受益。
- **【导向】** `pyproject.toml` 的 `[project.scripts] dimreflex = "dimreflex.cli:main"` 与 PyInstaller 打包目标都指向这里;exit code 被 validate 类脚本与 CI 语义消费。

---

## 4.2 `launcher.py` — PyInstaller 薄启动器

```python
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dimreflex.cli import main
if __name__ == "__main__":
    rc = main()
    if len(sys.argv) <= 1:
        with contextlib.suppress(EOFError):
            input("\n按回车键退出...")
    sys.exit(rc)
```

- **【计算】** 把父目录插入 sys.path 后导入并执行 `cli.main`;双击运行(无参数)结束后等一次回车。
- **【原因】** 扁平布局下脚本位于包目录内部,直接运行时 Python 搜索路径不含父目录——不插 sys.path 就 import 不到 dimreflex;等待回车避免双击时结果窗口一闪而过(EOFError 抑制:管道下 stdin 关闭也不报错)。
- **【导向】** `pyinstaller --onefile --name dimreflex_app launcher.py` → `dist/dimreflex_app.exe`(CONTEXT.md:PyInstaller 6.22.2 重打包并冒烟通过,量纲模式 Level 3 + 四族候选)。

---

## 4.3 `__init__.py` — 公共 API 面

- **【计算】** 从 9 个子模块重导出 31 个符号(所有公开函数与类),`__version__ = "0.5.0"`,`__all__` 白名单。
- **【原因】** `import dimreflex as dr; dr.discover(...)` 的命名空间即项目地图——`__all__` 顺序(异常→autodiscovery→core→discover→families→monitoring→scenarios→templates→uncertainty)本身就是模块依赖关系的镜像。
- **【导向】** README 全部示例、validate_*.py、测试 conftest 的 fallback 导入都走这个命名空间。

---

## 4.4 `pyproject.toml` — 打包与质量门

- **【计算/内容】**
  - `[build-system]`:setuptools ≥68;
  - `[project]`:name/version 0.5.0/requires-python ≥3.10;依赖钉上下界(numpy≥1.21,<3;pandas≥1.3,<4;scipy≥1.7,<2;pint≥0.19,<1);可选依赖 test=[pytest≥7];脚本入口 `dimreflex = "dimreflex.cli:main"`;
  - `[tool.setuptools]`:**`package-dir = { "dimreflex" = "." }`**——扁平布局:仓库根本身即包,不移动任何文件即可 `pip install -e .`(注释说明 2026-05 的包布局重构曾被 git 回退,本方案绕开);
  - `[tool.ruff.lint]`:14 个规则族全开(F/W/I/UP/B/SIM/C4/BLE/RUF/S/PLC/PLW/PIE/ISC/FURB);显式豁免及**逐条理由**:RUF001-003(中文标点是项目语言)、PLC0415(函数内延迟导入是刻意设计:pysr 可选依赖、CLI 门面隔离重型导入);per-file-ignores:tests/validation 的 S101(断言是测试/预注册判据惯用法,与退出码联动)、experiments 的 S101/S603/S607(实验脚本从 PATH 调 git 审计);
  - `[tool.pytest.ini_options]`:`testpaths = ["tests"]`——避免把实验归档中的历史样例误当成正式测试。
- **【原因】** 依赖范围钉上下界 + `requirements-lock.txt` 锁实测版本(2026-09-04,Python 3.13)是可复现部署的两层;ruff 全绿(371 项修复)是 CONTEXT.md 记录的质量里程碑;豁免逐条注明理由对应"不用全局静默"的纪律。
- **【导向】** `pip install -e .` → editable 导入路径;`dimreflex` 命令;`pytest` 只扫 tests/;CI(未建)可直接复用此配置。

---

## 4.5 门面层小结

四个文件没有任何独立计算,但有三处"哲学落地"值得记录:①CLI 无 schema 时的诚实声明(不伪装量纲分析);②launcher 的 sys.path 兜底与双击等待(可用性工程);③pyproject 的 package-dir 扁平映射(用打包配置而不是目录搬家解决"根即包")。
