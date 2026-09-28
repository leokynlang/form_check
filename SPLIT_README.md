# form_check — 拟合优度与函数正确性:必然因果不成立的证据仓

**整体目标(2026-09-20 修订,用户指令)**:不是为了证明发现了一个新的检测
架构,而是为了证明**函数结果与函数拟合匹配度本身不构成必然因果**——
拟合匹配度(R²、留出 CV、AICc、RMSE、参数标准误)高,不必然推出拟合
函数正确,也不构成正确性证书;匹配度低同样不必然推出会被现有检验查出。

- 主证据链(2026-09-20):**`no_causal_link.py` / `no_causal_link.md`(权威
  一键量化,12 项锚点核对,两侧断言操作化)**、`reset_fitting_quality_report.md`
  (合成 57,000 + SRSD 真实公式 240,000 数据集)与
  `fitting_quality_vs_function_error.md`(域内优度 vs 域外真值偏差直检);
- **form_check vs RESET 公平比较(本仓拆分时的原始用途)降位为测量仪器**:
  四门功效用于刻画「优度不含形式信息」的边界与机制(方向对齐度才是检出
  决定量),检测器不是贡献主体;
- 诚实口径:合成台成立且显著(CUBIC 域内 R²=0.9983、域外中位 44%;
  R²≥0.99 错形式 100% 域外偏差>10%),SRSD 真实公式部分成立(例外
  bonus.5;两端盲区 ii.11.28/i.41.16);偏差必须定义在域外/函数行为上,
  域内平均偏差与 R² 数学绑定。

以下为拆分与运行说明(原文保留):

2026-09-05 自 `D:\dimreflex`(head 554f861, 工作区 dirty)拆分。
本文件夹**自包含**:根目录即原 `experiments/form_check_vs_reset/`,
dimreflex 框架以 vendored 包形式位于 `./dimreflex/`(核心 `discover.py`
sha256 `bdf7b069…`,与拆分时主仓一致)。无任何跨文件夹引用。

## 独立运行

```bash
pip install -r requirements.txt        # numpy/pandas/scipy/pint
"$PYTHON" run_experiment.py --tag pilot --seeds 50     # 合成 19 场景
"$PYTHON" run_srsd.py --tag srsd_pilot                 # SRSD-Feynman 32 任务
"$PYTHON" analyze.py pilot
"$PYTHON" run_battery.py                               # 3-2-1 电池
"$PYTHON" validation/pilot_form_check_vs_reset.py      # 1h pilot(exit 0 预期)
cd tests && python -m pytest -q                        # 144 测试 + 电池测试
```

注意:
- 本机若 pip 安装过 dimreflex(editable),`_bootstrap.py` 的元路径查找器
  会强制使用本目录 vendored 版,不会被遮蔽。
- 裸导入 `scenarios`/`common` 命中**本目录**的实验版;vendored 核心包内
  也有同名 `scenarios.py`(多场景拟合,不同物),由 sys.path 顺序隔离。
- `capture_env()` 的 `git_head` 记空串(本文件夹非 git 仓库,git 命令非零退出未被检查),属预期。
- `SRSD` 上游浅克隆到 `%TEMP%/srsd-benchmark`(`srsd_data.py` 自动)。
- `paper/`、`artifacts/`、`report_parts/` 为论文与知识件;构建脚本
  (`build_hx_paper.py`、`artifacts/build_unified_paper.py`)已改为脚本相对路径。
- `runs/` 大文件为产物,可经上述 runner 重生成;结论件(report/summary)保留。
