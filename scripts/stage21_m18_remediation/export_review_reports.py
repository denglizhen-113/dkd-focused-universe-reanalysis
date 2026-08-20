from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DOC_DIR = ROOT / "docs" / "stage21_m18_remediation"
EXPORTER_PATH = ROOT / "scripts" / "stage21_m18_remediation" / "export_m18_documents.py"


def load_exporter():
    spec = importlib.util.spec_from_file_location("m18_exporter_for_reports", EXPORTER_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(EXPORTER_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    exporter = load_exporter()
    reports = [
        DOC_DIR / "TWO_REVIEWER_RESPONSE_MATRIX.md",
        DOC_DIR / "双审稿报告严格整改总报告.md",
    ]
    for markdown in reports:
        docx = markdown.with_suffix(".docx")
        pdf = markdown.with_suffix(".pdf")
        exporter.to_docx(markdown, docx)
        exporter.to_pdf(docx, pdf)
        print(f"exported={pdf}")


if __name__ == "__main__":
    main()
