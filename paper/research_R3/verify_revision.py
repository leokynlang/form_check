"""Verify R3's evidence links, preserved legacy tables and rendered PDF.

Run with the same PDF-capable Python runtime used by build_pdf.py.
Visual inspection is recorded separately, against the final PDF hash.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pdfplumber
from pypdf import PdfReader

PAPER = Path(__file__).resolve().parent
ROOT = PAPER.parent.parent
PDF = PAPER / "Silent_Acceptance_research_R3_2026-10-18.pdf"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def require(condition, message: str):
    if not condition:
        raise ValueError(message)


def main():
    manuscript = (PAPER / "manuscript.md").read_text(encoding="utf-8")
    legacy = (ROOT / "paper/research_R2/manuscript.md").read_text(encoding="utf-8")
    old_rows = [line for line in legacy.splitlines() if line.startswith("|")]
    new_rows = [line for line in manuscript.splitlines() if line.startswith("|")]
    # The historical row sequence must remain intact, including repeated headers.
    cursor = iter(new_rows)
    for row in old_rows:
        require(any(candidate == row for candidate in cursor), "R2 table row changed or reordered")
    facts = read_json(PAPER / "evidence/revision_facts.json")
    require("pending" not in facts["status"], "R3 evidence is still pending")
    for path, expected in facts["source_sha256"].items():
        require(digest(ROOT / path) == expected, f"Changed R3 evidence source: {path}")
    audits = ["runs/p02_fixed_domain/stage1/independent_audit.json",
              "runs/p02_fixed_domain/stage2/independent_audit.json",
              "runs/p02_fixed_domain/complete_grid/independent_audit.json",
              "runs/p03_physical_protocol/independent_selection_review.json"]
    for path in audits:
        require(read_json(ROOT / path)["passed"], f"Audit did not pass: {path}")
    snapshot = read_json(ROOT / "runs/p02_fixed_domain/stage2_preservation_before.json")
    for key, directory in [("stage1_files", ROOT / "runs/p02_fixed_domain/stage1"),
                           ("paper_R2_files", ROOT / "paper/research_R2")]:
        current = {p.relative_to(directory).as_posix(): digest(p)
                   for p in directory.rglob("*") if p.is_file()}
        require(current == snapshot[key], f"Preserved tree changed: {key}")

    reader = PdfReader(PDF)
    pdf_text = "\n".join(page.extract_text() or "" for page in reader.pages)
    for marker in ["P02_FINAL_AUDIT_STATUS", "TODO", "TBD"]:
        require(marker not in manuscript and marker not in pdf_text, f"Unresolved marker: {marker}")
    for text in ["246,240", "108,000", "Table 18", "Table 19", "Figure D.6", "Figure D.7"]:
        require(text in pdf_text, f"Missing R3 content in PDF: {text}")
    with pdfplumber.open(PDF) as document:
        for number, page in enumerate(document.pages, 1):
            require(abs(page.width - 595.28) < 2 and abs(page.height - 841.89) < 2,
                    f"Unexpected page size: {number}")
            for char in page.chars:
                require(char["x0"] >= -.5 and char["x1"] <= page.width + .5
                        and char["top"] >= -.5 and char["bottom"] <= page.height + .5,
                        f"Character outside PDF page: {number}")
                require(char["text"] not in "\ufffd\u25a0\u25a1\x00",
                        f"Suspect missing glyph on page: {number}")
    visual = read_json(PAPER / "evidence/visual_review.json")
    require(visual["pdf_sha256"] == digest(PDF) and visual["passed"],
            "Visual review is missing or belongs to a different PDF")
    result = {"passed": True, "pdf_sha256": digest(PDF),
              "manuscript_sha256": digest(PAPER / "manuscript.md"),
              "pages": len(reader.pages), "preserved_legacy_table_rows": len(old_rows),
              "evidence_source_hashes_verified": len(facts["source_sha256"]),
              "audits": {path: digest(ROOT / path) for path in audits},
              "preserved_trees": {key: len(snapshot[key])
                                  for key in ("stage1_files", "paper_R2_files")},
              "page_bounds_and_glyphs": "passed", "unresolved_markers": 0,
              "visual_review_sha256": digest(PAPER / "evidence/visual_review.json")}
    (PAPER / "evidence/revision_validation.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
