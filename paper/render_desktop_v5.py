# -*- coding: utf-8 -*-
"""render_desktop_v5.py — 桌面 Silent_Acceptance_v5.md → v5.pdf 渲染。

字形修复:∅→Ø;∈/∉/≡ 经 SimHei;上/下标 unicode→标签(与 build_pdf 同表)。
独立命名,不与平行会话的 render_v5_pdf.py 冲突。
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
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table)

SRC = os.path.join(os.path.expanduser("~"), "Desktop",
                   "Silent_Acceptance_v5.md")
OUT = os.path.join(os.path.expanduser("~"), "Desktop",
                   "Silent_Acceptance_v5.pdf")

FDIR = "C:/Windows/Fonts/"
pdfmetrics.registerFont(TTFont("Times New Roman", FDIR + "times.ttf"))
pdfmetrics.registerFont(TTFont("Times New Roman-Bold", FDIR + "timesbd.ttf"))
pdfmetrics.registerFont(TTFont("Times New Roman-Italic", FDIR + "timesi.ttf"))
pdfmetrics.registerFont(TTFont("SimHei", FDIR + "simhei.ttf"))
pdfmetrics.registerFont(TTFont("Consolas", FDIR + "consola.ttf"))
registerFontFamily("Times New Roman", normal="Times New Roman",
                   bold="Times New Roman-Bold",
                   italic="Times New Roman-Italic",
                   boldItalic="Times New Roman-Bold")
registerFontFamily("SimHei", normal="SimHei", bold="SimHei")

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


def U(*codes):
    return "".join(chr(c) for c in codes)


def inline(t):
    t = esc(t)
    t = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", t)
    for pat, rep in [
        (U(0x207B, 0xB3, 0xB4), "<super>" + U(0x2212) + "34</super>"),
        (U(0x207B, 0xB2, 0xB3), "<super>" + U(0x2212) + "23</super>"),
        (U(0x207B, 0xB9, 0xB2), "<super>" + U(0x2212) + "12</super>"),
        (U(0x207B) + U(0xB9), "<super>" + U(0x2212) + "1</super>"),
        (U(0xB2), "<super>2</super>"),
        (U(0xB3), "<super>3</super>"),
        (U(0x1D62), "<sub>i</sub>"),
        (U(0x2C7C), "<sub>j</sub>"),
        (U(0x2081), "<sub>1</sub>"),
        (U(0x2082), "<sub>2</sub>"),
        (U(0x2096), "<sub>k</sub>"),
        (U(0x27F9), U(0x2192)),
    ]:
        t = t.replace(pat, rep)
    # 逐位上标兜底:多位指数(如 ⁻¹⁵)整体进 <super>,基线一致
    for d, name in [(0x2070, '0'), (0xB9, '1'), (0xB2, '2'), (0xB3, '3'),
                    (0x2074, '4'), (0x2075, '5'), (0x2076, '6'),
                    (0x2077, '7'), (0x2078, '8'), (0x2079, '9')]:
        t = t.replace(U(d), '<super>' + name + '</super>')
    t = t.replace(U(0x207B), '<super>' + U(0x2212) + '</super>')
    for ch in (U(0x2208), U(0x2209), U(0x2261)):
        t = t.replace(ch, "<font name='SimHei'>" + ch + "</font>")
    t = t.replace(U(0x2205), U(0x00D8))
    t = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", t)
    t = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<i>\1</i>", t)
    t = re.sub(r"`([^`]+)`", r"<font face='Consolas' size='8.5'>\1</font>", t)
    return t


def make_table(rows):
    ncol = max(len(r) for r in rows)
    rows = [r + [""] * (ncol - len(r)) for r in rows]
    maxlen = [max(len(r[c]) for r in rows) for c in range(ncol)]
    weights = [max(8.0, float(min(m, 55))) for m in maxlen]

    def longest(c):
        return max((len(w) for r in rows for w in r[c].split()), default=0)

    floors = [min(14, longest(c) + 2) for c in range(ncol)]
    weights = [max(f, w) for f, w in zip(floors, weights)]
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
    lines = open(SRC, encoding="utf-8").read().split("\n")
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
                            title="Silent Acceptance — v5",
                            author="Anonymous")
    doc.build(story)
    print("ok ->", OUT)


if __name__ == "__main__":
    main()
