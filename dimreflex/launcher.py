"""PyInstaller 打包入口:双击/命令行运行 dimreflex 的薄启动器。

扁平布局下,脚本位于包目录内部,直接运行时 Python 的模块搜索路径
不包含父目录,因此先把父目录插入 sys.path 再导入。
双击(无参数)运行时停留等待回车,避免结果窗口一闪而过;
终端里带参数运行则直接返回。
"""

import contextlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dimreflex.cli import main

if __name__ == "__main__":
    rc = main()
    if len(sys.argv) <= 1:
        with contextlib.suppress(EOFError):
            input("\n按回车键退出...")
    sys.exit(rc)
