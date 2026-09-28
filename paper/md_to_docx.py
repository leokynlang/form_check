"""把论文 Markdown 转成格式化 Word 文档。

支持标题、粗体、斜体、行内代码、表格、注解块、分隔线和列表。
默认转换当前目录的 ``draft_v1_complete.md``；也可从命令行指定任意输入和输出。
"""

import argparse
import re
from pathlib import Path

from docx import Document
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt, RGBColor

HEADING_COLOR = RGBColor(0x1F, 0x3B, 0x5C)


def add_runs(paragraph, text: str) -> None:
    """把简单 Markdown 行内标记添加到 Word 段落。"""
    tokens = re.split(r"(\*\*.+?\*\*|\*[^*]+?\*|`[^`]+?`)", text)
    for token in tokens:
        if not token:
            continue
        if token.startswith("**") and token.endswith("**"):
            run = paragraph.add_run(token[2:-2])
            run.bold = True
        elif token.startswith("*") and token.endswith("*") and len(token) > 2:
            run = paragraph.add_run(token[1:-1])
            run.italic = True
        elif token.startswith("`") and token.endswith("`"):
            run = paragraph.add_run(token[1:-1])
            run.font.name = "Consolas"
            run.font.color.rgb = RGBColor(0x88, 0x22, 0x22)
        else:
            paragraph.add_run(token)


def shade(paragraph, color: str = "EFE7DD") -> None:
    """给段落添加底纹。"""
    properties = paragraph._p.get_or_add_pPr()
    shading = OxmlElement("w:shd")
    shading.set(qn("w:val"), "clear")
    shading.set(qn("w:fill"), color)
    properties.append(shading)


def add_heading(document, text: str, level: int) -> None:
    heading = document.add_heading(level=level)
    add_runs(heading, text)
    for run in heading.runs:
        run.font.color.rgb = HEADING_COLOR


def add_table(document, rows: list[list[str]]) -> None:
    if not rows:
        return

    column_count = max(len(row) for row in rows)
    table = document.add_table(rows=len(rows), cols=column_count)
    table.style = "Light Grid Accent 1"
    for row_index, row in enumerate(rows):
        for column_index in range(column_count):
            cell = table.cell(row_index, column_index)
            cell.text = ""
            text = row[column_index] if column_index < len(row) else ""
            add_runs(cell.paragraphs[0], text)
            for run in cell.paragraphs[0].runs:
                run.font.size = Pt(9)
                if row_index == 0:
                    run.bold = True
    document.add_paragraph()


def convert_markdown(source: Path, destination: Path) -> Path:
    """转换一个 Markdown 文件并返回生成的 Word 路径。"""
    document = Document()
    style = document.styles["Normal"]
    style.font.name = "Calibri"
    style.font.size = Pt(10.5)
    style.element.rPr.rFonts.set(qn("w:eastAsia"), "微软雅黑")

    lines = source.read_text(encoding="utf-8").splitlines()
    index = 0
    while index < len(lines):
        text = lines[index].strip()
        if not text:
            index += 1
            continue

        if (
            text.startswith("|")
            and index + 1 < len(lines)
            and lines[index + 1].strip().startswith("|")
        ):
            rows = []
            while index < len(lines) and lines[index].strip().startswith("|"):
                cells = [
                    cell.strip()
                    for cell in lines[index].strip().strip("|").split("|")
                ]
                if not all(re.fullmatch(r":?-{2,}:?", cell) for cell in cells):
                    rows.append(cells)
                index += 1
            add_table(document, rows)
            continue

        if text.startswith("# "):
            add_heading(document, text[2:], 0)
        elif text.startswith("## "):
            add_heading(document, text[3:], 1)
        elif text.startswith("### "):
            add_heading(document, text[4:], 2)
        elif text.startswith("---") and set(text) <= {"-"}:
            document.add_paragraph()
        elif text.startswith(("📝", "**Annotated draft")):
            paragraph = document.add_paragraph()
            shade(paragraph)
            add_runs(paragraph, text)
        elif text.startswith("- "):
            paragraph = document.add_paragraph(style="List Bullet")
            add_runs(paragraph, text[2:])
        elif re.match(r"^\d+\.\s", text):
            paragraph = document.add_paragraph(style="List Number")
            add_runs(paragraph, re.sub(r"^\d+\.\s", "", text))
        else:
            paragraph = document.add_paragraph()
            add_runs(paragraph, text)
        index += 1

    destination.parent.mkdir(parents=True, exist_ok=True)
    document.save(destination)
    return destination


def main(argv: list[str] | None = None) -> int:
    paper_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "source",
        nargs="?",
        type=Path,
        default=paper_dir / "draft_v1_complete.md",
        help="待转换的 Markdown 文件",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="输出 .docx 路径；默认与输入文件同名",
    )
    args = parser.parse_args(argv)
    output = args.output or args.source.with_suffix(".docx")
    generated = convert_markdown(args.source.resolve(), output.resolve())
    print(f"saved {generated}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
