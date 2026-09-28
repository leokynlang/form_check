# -*- coding: utf-8 -*-
"""render_v5_pdf.py — draft_v5_complete.md 的简版 PDF 渲染(纯阅读用)。

复用 render_draft_pdf.py 的管线,仅切换源文件与输出名。
"""
import os
import sys

PAPER_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PAPER_DIR)

import render_draft_pdf as base

base.MD_NAME = "draft_v5_complete.md"
base.OUT = os.path.join(PAPER_DIR, "draft_v5_complete.pdf")
base.main() if hasattr(base, "main") else None
