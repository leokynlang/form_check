"""cli.py 测试:命令行入口的契约。

用法一:dimreflex_app <csv> <目标列> <特征列...>  → 分析文件,退出码 0
用法二:dimreflex_app(无参数)                  → 内置双工况演示,退出码 0
用法错误                                       → 退出码 2,打印用法
"""

import json
from pathlib import Path

import pytest
from dimreflex import cli

PROJECT_ROOT = Path(__file__).resolve().parent.parent
HEAT_CSV = PROJECT_ROOT / "heat_exchanger_data.csv"


def _write_heat_csv_and_schema(tmp_path, csv_name="heat.csv"):
    csv = tmp_path / csv_name
    csv.write_bytes(HEAT_CSV.read_bytes())
    schema = {
        "system": "tube_side_heat_transfer",
        "target": "h",
        "variables": [
            {"name": "rho", "unit": "kg/m**3",
             "role": "material_property"},
            {"name": "mu", "unit": "Pa*s", "role": "material_property"},
            {"name": "cp", "unit": "J/(kg*K)",
             "role": "material_property"},
            {"name": "k", "unit": "W/(m*K)", "role": "material_property"},
            {"name": "D", "unit": "m", "role": "geometry"},
            {"name": "v", "unit": "m/s", "role": "state"},
            {"name": "h", "unit": "W/(m**2*K)", "role": "target"},
        ],
    }
    (tmp_path / (csv.stem + ".schema.json")).write_text(
        json.dumps(schema, ensure_ascii=False), encoding="utf-8")
    return csv


@pytest.fixture
def mini_csv(tmp_path):
    """单特征幂律小数据:Nu = 0.023·Re^0.8。"""
    csv = tmp_path / "mini.csv"
    csv.write_text(
        "Re,Nu\n"
        "5000,23.0\n8000,34.1\n12000,47.5\n20000,70.2\n"
        "35000,110.4\n60000,173.9\n90000,245.6\n",
        encoding="utf-8",
    )
    return csv


class TestCliMain:
    def test_csv_mode_reports_power_law(self, mini_csv, capsys):
        rc = cli.main([str(mini_csv), "Nu", "Re"])
        out = capsys.readouterr().out
        assert rc == 0
        assert "power_law" in out          # 幂律族应上榜
        assert "Nu" in out                 # 公式里应有目标列名
        assert "R²" in out                 # 报告应带精度

    def test_demo_mode_runs_all_families(self, capsys):
        rc = cli.main([])
        out = capsys.readouterr().out
        assert rc == 0
        assert "piecewise_power_law" in out  # 演示数据是双工况,分段族应出现
        assert "power_law" in out

    def test_single_arg_infers_roles(self, mini_csv, capsys):
        """只给 CSV 路径(=拖拽文件到 exe 上):最后一列作目标,其余作特征。"""
        rc = cli.main([str(mini_csv)])
        out = capsys.readouterr().out
        assert rc == 0
        assert "Nu" in out and "power_law" in out

    def test_single_arg_with_quoted_path(self, mini_csv, capsys):
        """拖拽进来的路径可能带引号,必须剥掉。"""
        rc = cli.main([f'"{mini_csv}"'])
        assert rc == 0

    def test_no_schema_mode_is_honest_generic(self, mini_csv, capsys):
        """没有 schema 时不得伪装成量纲分析;输出须显式声明裸拟合模式。"""
        rc = cli.main([str(mini_csv), "Nu", "Re"])
        out = capsys.readouterr().out
        assert rc == 0
        assert "非量纲感知模式" in out
        assert "power_law" in out

    def test_schema_sidecar_activates_dimensional_mode(self, tmp_path, capsys):
        csv = _write_heat_csv_and_schema(tmp_path)
        rc = cli.main([str(csv)])
        out = capsys.readouterr().out
        assert rc == 0
        assert "量纲感知模式" in out
        assert "power_law" in out

    def test_explicit_schema_flag_activates_dimensional_mode(
            self, tmp_path, capsys):
        csv = tmp_path / "data.csv"
        csv.write_bytes(HEAT_CSV.read_bytes())
        schema_path = tmp_path / "units.json"
        schema_path.write_text(
            json.dumps({
                "system": "heat", "target": "h",
                "variables": [
                    {"name": "rho", "unit": "kg/m**3",
                     "role": "material_property"},
                    {"name": "mu", "unit": "Pa*s",
                     "role": "material_property"},
                    {"name": "cp", "unit": "J/(kg*K)",
                     "role": "material_property"},
                    {"name": "k", "unit": "W/(m*K)",
                     "role": "material_property"},
                    {"name": "D", "unit": "m", "role": "geometry"},
                    {"name": "v", "unit": "m/s", "role": "state"},
                    {"name": "h", "unit": "W/(m**2*K)",
                     "role": "target"},
                ],
            }, ensure_ascii=False), encoding="utf-8")
        rc = cli.main(["--schema", str(schema_path), str(csv)])
        out = capsys.readouterr().out
        assert rc == 0
        assert "量纲感知模式" in out
        assert "power_law" in out

    def test_bad_schema_exits_1(self, tmp_path, capsys):
        csv = tmp_path / "heat.csv"
        csv.write_bytes(HEAT_CSV.read_bytes())
        schema_path = tmp_path / "heat.schema.json"
        schema_path.write_text("{", encoding="utf-8")
        rc = cli.main([str(csv)])
        assert rc == 1
        assert "schema" in capsys.readouterr().out

    def test_two_args_infers_features(self, mini_csv, capsys):
        """给路径+目标列:其余列自动作为特征。"""
        rc = cli.main([str(mini_csv), "Nu"])
        out = capsys.readouterr().out
        assert rc == 0
        assert "Re" in out  # 特征清单里应出现被自动纳入的 Re

    def test_single_column_csv_returns_usage_error(self, tmp_path, capsys):
        """只有一列的数据无从谈'特征→目标',属用法错误(rc=2)。"""
        csv = tmp_path / "one_col.csv"
        csv.write_text("x\n1\n2\n3\n", encoding="utf-8")
        rc = cli.main([str(csv)])
        assert rc == 2

    def test_interactive_enter_runs_demo(self, monkeypatch, capsys):
        """双击后直接回车 → 内置演示。"""
        import io
        monkeypatch.setattr(cli.sys, "stdin",
                            io.StringIO("\n"))
        monkeypatch.setattr(cli.sys.stdin, "isatty", lambda: True)
        rc = cli.main([])
        out = capsys.readouterr().out
        assert rc == 0
        assert "piecewise_power_law" in out

    def test_interactive_path_runs_analysis(self, monkeypatch, capsys,
                                            mini_csv):
        """双击后输入 CSV 路径 → 直接分析该文件。"""
        import io
        monkeypatch.setattr(cli.sys, "stdin",
                            io.StringIO(f"{mini_csv}\n"))
        monkeypatch.setattr(cli.sys.stdin, "isatty", lambda: True)
        rc = cli.main([])
        out = capsys.readouterr().out
        assert rc == 0
        assert "power_law" in out

    def test_piped_stdin_skips_prompt_runs_demo(self, capsys):
        """stdin 被重定向(自动化/管道调用)时不得停下来等输入。"""
        rc = cli.main([])
        out = capsys.readouterr().out
        assert rc == 0
        assert "piecewise_power_law" in out

    def test_interactive_two_bad_paths_exit_2(self, monkeypatch, capsys):
        """交互模式连续两次输入无效路径 → 打印用法,退出码 2。"""
        import io
        monkeypatch.setattr(cli.sys, "stdin",
                            io.StringIO("不存在的A.csv\n不存在的B.csv\n"))
        monkeypatch.setattr(cli.sys.stdin, "isatty", lambda: True)
        rc = cli.main([])
        assert rc == 2
        assert "用法" in capsys.readouterr().out

    def test_missing_file_reports_error(self, tmp_path, capsys):
        rc = cli.main([str(tmp_path / "不存在.csv"), "Nu", "Re"])
        out = capsys.readouterr().out
        assert rc == 1
        assert "不存在" in out or "找不到" in out or "Error" in out

    def test_gbk_console_does_not_crash(self, mini_csv, monkeypatch):
        """回归:中文 Windows 控制台默认 GBK,'✓' 等符号曾致崩溃。

        模拟 GBK stdout,报告必须仍能输出(换码/替换均可,绝不抛异常)。
        """
        import io

        gbk_out = io.TextIOWrapper(io.BytesIO(), encoding="gbk",
                                   line_buffering=True)
        monkeypatch.setattr(cli.sys, "stdout", gbk_out)
        rc = cli.main([str(mini_csv), "Nu", "Re"])
        assert rc == 0  # 能走到这里就说明没被 UnicodeEncodeError 炸掉
