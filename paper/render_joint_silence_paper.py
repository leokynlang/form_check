# -*- coding: utf-8 -*-
"""render_draft_pdf.py — 论文 markdown 工作稿的简版 PDF 渲染(纯阅读用)。

与 build_pdf.py(投稿版:封面/TOC/图/代码清单)不同,本脚本把
draft_v3_complete.md 原样排成可读 PDF,含变更记录等过程性内容。
输出: draft_v3_complete.pdf
"""
import os
import re

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase.pdfmetrics import registerFontFamily
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                KeepTogether)

PAPER_DIR = os.path.dirname(os.path.abspath(__file__))
MD_NAME = os.path.join(os.path.expanduser("~"), "Desktop", "Silent_Acceptance_JointSilence_paper.md")
OUT = os.path.join(os.path.expanduser("~"), "Desktop", "Silent_Acceptance_JointSilence_paper.pdf")

FDIR = "C:/Windows/Fonts/"
pdfmetrics.registerFont(TTFont("Times New Roman", FDIR + "times.ttf"))
pdfmetrics.registerFont(TTFont("Times New Roman-Bold", FDIR + "timesbd.ttf"))
pdfmetrics.registerFont(TTFont("Times New Roman-Italic", FDIR + "timesi.ttf"))
pdfmetrics.registerFont(TTFont("Consolas", FDIR + "consola.ttf"))
registerFontFamily("Times New Roman", normal="Times New Roman",
                   bold="Times New Roman-Bold", italic="Times New Roman-Italic",
                   boldItalic="Times New Roman-Bold")

PAGE_W, PAGE_H = A4
MARGIN = 0.9 * inch
TEXT_W = PAGE_W - 2 * MARGIN

S = {
    "title": ParagraphStyle("title", fontName="Times New Roman-Bold",
                            fontSize=15, leading=19, spaceAfter=10),
    "h1": ParagraphStyle("h1", fontName="Times New Roman-Bold", fontSize=12.5,
                         leading=16, spaceBefore=14, spaceAfter=5),
    "h2": ParagraphStyle("h2", fontName="Times New Roman-Bold", fontSize=11,
                         leading=14, spaceBefore=10, spaceAfter=4),
    "body": ParagraphStyle("body", fontName="Times New Roman", fontSize=9.5,
                           leading=13, spaceAfter=5, alignment=4),
    "cell": ParagraphStyle("cell", fontName="Times New Roman", fontSize=8,
                           leading=10),
    "list": ParagraphStyle("list", fontName="Times New Roman", fontSize=9.5,
                           leading=13, spaceAfter=3, leftIndent=14,
                           alignment=4),
}


def esc(t):
    return t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def inline(t):
    t = esc(t)
    t = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", t)          # links -> text
    t = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", t)          # bold
    t = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<i>\1</i>", t)  # italic
    t = re.sub(r"`([^`]+)`", r"<font face='Consolas' size='8.5'>\1</font>", t)
    return t


def make_table(rows):
    ncol = max(len(r) for r in rows)
    rows = [r + [""] * (ncol - len(r)) for r in rows]
    # 列宽按各列最长单元格文本量加权(上限截断),让 Reading 类长文本列变宽
    weights = [min(max(len(r[c]) for r in rows), 60) + 4 for c in range(ncol)]
    tot = sum(weights)
    widths = [TEXT_W * w / tot for w in weights]
    data = [[Paragraph(inline(c), S["cell"]) for c in r] for r in rows]
    tbl = Table(data, colWidths=widths, repeatRows=1)
    tbl.setStyle([("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#999999")),
                  ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#eeeeee")),
                  ("VALIGN", (0, 0), (-1, -1), "TOP"),
                  ("LEFTPADDING", (0, 0), (-1, -1), 3),
                  ("RIGHTPADDING", (0, 0), (-1, -1), 3)])
    return tbl


def main():
    lines = open(os.path.join(PAPER_DIR, MD_NAME), encoding="utf-8").read()
    lines = lines.splitlines()
    story = []
    i = 0
    while i < len(lines):
        ln = lines[i]
        if not ln.strip() or ln.strip() == "---":
            i += 1
            continue
        if ln.startswith("# ") and not story:
            story.append(Paragraph(inline(ln[2:]), S["title"]))
            i += 1
            continue
        if ln.startswith("## "):
            story.append(Paragraph(inline(ln[3:]), S["h1"]))
            i += 1
            continue
        if ln.startswith("### "):
            story.append(Paragraph(inline(ln[4:]), S["h2"]))
            i += 1
            continue
        if ln.startswith("|"):
            rows = []
            while i < len(lines) and lines[i].startswith("|"):
                cells = [c.strip() for c in
                         lines[i].strip().strip("|").split("|")]
                if not all(re.fullmatch(r":?-{3,}:?", c) for c in cells):
                    rows.append(cells)
                i += 1
            story.append(Spacer(1, 4))
            story.append(make_table(rows))
            story.append(Spacer(1, 6))
            continue
        m = re.match(r"^(\d+)\.\s+(.*)$", ln)
        if m:
            story.append(Paragraph(f"{m.group(1)}.&nbsp;&nbsp;"
                                   + inline(m.group(2)), S["list"]))
            i += 1
            continue
        if ln.startswith("- "):
            story.append(Paragraph("•&nbsp;&nbsp;" + inline(ln[2:]),
                                   S["list"]))
            i += 1
            continue
        # 普通段落:连续非空行合并
        buf = [ln]
        i += 1
        while (i < len(lines) and lines[i].strip()
               and not re.match(r"^(#|\||\d+\.\s|- )", lines[i])):
            buf.append(lines[i])
            i += 1
        story.append(Paragraph(inline(" ".join(buf)), S["body"]))

    doc = SimpleDocTemplate(OUT, pagesize=A4,
                            leftMargin=MARGIN, rightMargin=MARGIN,
                            topMargin=0.85 * inch, bottomMargin=0.85 * inch,
                            title="Silent Acceptance — draft v5 (working)")
    doc.build(story)
    print("ok ->", OUT)


if __name__ == "__main__":
    main()
