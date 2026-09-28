# -*- coding: utf-8 -*-
"""make_mindmap.py — 论文思维导图(英文,v3 口径),嵌入 A4 PDF 用。
布局:中心节点 + 左右各 3 个分支;分支 = 标题块 + 叶子芯片列。
字号按嵌入比例换算:图宽 8.6in,PDF 内容宽约 6.3in → 缩放 0.73;
叶子 10.5pt → 有效 7.7pt,分支标题 13.5pt → 9.9pt。
每分支叶数 ≤ 4,叶文案经过宽度预算(芯片约 195pt)。
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

plt.rcParams["axes.unicode_minus"] = False

BRANCHES = [
    # (title, leaves, color, side, ytop)
    ("Research question", [
        "Good fit \u2260 correct form",
        "Metrics silently accept",
        "Breaks extrapolation",
    ], "#297189", "R", 0.015),
    ("Gates under test (A\u2013F)", [
        "A permutation curvature",
        "B RESET fitted-F",
        "C\u00b7D RESET-HC3 variants",
        "E\u00b7F CV-quad / AICc-quad",
    ], "#547ea9", "R", 0.30),
    ("Key findings", [
        "HC3 oversized 42\u201356%",
        "AICc-quad: best power",
        "Blind spots complementary",
        "Alarm enrichment 2\u20132.8\u00d7",
    ], "#ac5047", "R", 0.63),
    ("Protocol P1\u2013P5", [
        "Common \u03b1 = 0.05 basis",
        "Power at matched 5% FPR",
        "Preregistered labels",
        "Typed failure registration",
    ], "#489261", "L", 0.015),
    ("Six task suites", [
        "S1 synthetic \u00b7 19 scenarios",
        "S2 SRSD-Feynman \u00b7 32 tasks",
        "S3 family map \u00b7 S4 battery",
        "S5 gate cmp \u00b7 S6 real data",
    ], "#917b37", "L", 0.365),
    ("Battery & deployment", [
        "Viability set, not verdict",
        "Nulls coordinate-dependent",
        "Wild bootstrap repair",
        "OR union for recall",
    ], "#a58443", "L", 0.70),
]

FIG_W, FIG_H = 8.6, 6.0
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H), dpi=300)
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

CX0, CX1 = 0.370, 0.630
CY0, CY1 = 0.420, 0.602

L_X0, L_X1 = 0.015, 0.350
R_X0, R_X1 = 0.650, 0.985
LEAF_FS = 9.6
TH = 0.050
LH, GAP = 0.044, 0.015


def box(x0, y0, x1, y1, fc, ec, lw=1.2, r=0.016):
    p = FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                       boxstyle=f"round,pad=0,rounding_size={r}",
                       fc=fc, ec=ec, lw=lw, mutation_aspect=FIG_W / FIG_H)
    ax.add_patch(p)


def tint(hex_color, f=0.90):
    h = hex_color.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    mix = lambda c: int(c + (255 - c) * f)
    return f"#{mix(r):02x}{mix(g):02x}{mix(b):02x}"


# ── 中心节点 ──
box(CX0, CY0, CX1, CY1, tint("#297189", 0.86), "#297189", lw=2.0)
ax.text((CX0 + CX1) / 2, 0.556, "Form-Specification\nGates", ha="center",
        va="center", fontsize=11.0, fontweight="bold", color="#12262e")
ax.text((CX0 + CX1) / 2, 0.492, "Silent Acceptance", ha="center",
        va="center", fontsize=9.5, color="#297189", style="italic")
ax.text((CX0 + CX1) / 2, 0.460, "508,825 fitted datasets", ha="center",
        va="center", fontsize=8.8, color="#41616c")
ax.text((CX0 + CX1) / 2, 0.438, "6 gates \u00b7 6 suites", ha="center",
        va="center", fontsize=8.8, color="#41616c")

for title, leaves, color, side, ytop in BRANCHES:
    x0 = L_X0 if side == "L" else R_X0
    x1 = L_X1 if side == "L" else R_X1
    ty0 = 1 - ytop - TH
    box(x0, ty0, x1, ty0 + TH, tint(color, 0.80), color, lw=1.5)
    ax.text(x0 + (x1 - x0) / 2, ty0 + TH / 2, title, ha="center", va="center",
            fontsize=12.5, fontweight="bold", color="#222")
    hub_x = CX0 if side == "L" else CX1
    hub_y = (CY0 + CY1) / 2
    tip_x = x1 if side == "L" else x0
    ax.plot([hub_x, (hub_x + tip_x) / 2, tip_x],
            [hub_y, ty0 + TH / 2, ty0 + TH / 2],
            color=color, lw=1.6, solid_capstyle="round", zorder=0)
    # 叶子芯片
    fx0 = x0 + 0.022 if side == "L" else x0
    fx1 = x1 if side == "L" else x1 - 0.022
    spine_x = x0 + 0.010 if side == "L" else x1 - 0.010
    y = ty0 - GAP
    spine_top, spine_bot = y - 0.002, y - len(leaves) * (LH + GAP) + GAP + 0.002
    ax.plot([spine_x, spine_x], [spine_bot, spine_top],
            color=color, lw=1.0, zorder=0)
    for leaf in leaves:
        ly0 = y - LH
        box(fx0, ly0, fx1, ly0 + LH, tint(color, 0.94), color, lw=0.7, r=0.011)
        ax.text((fx0 + fx1) / 2, ly0 + LH / 2, leaf, ha="center", va="center",
                fontsize=LEAF_FS, color="#1a1a1a")
        lx = fx0 if side == "L" else fx1
        ax.plot([spine_x, lx], [ly0 + LH / 2, ly0 + LH / 2],
                color=color, lw=0.9, zorder=0)
        y = ly0 - GAP

fig.savefig("mindmap.png", dpi=300, bbox_inches="tight",
            facecolor="white", pad_inches=0.06)
print("saved mindmap.png")
