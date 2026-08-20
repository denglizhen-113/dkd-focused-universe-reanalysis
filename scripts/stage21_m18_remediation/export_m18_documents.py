from __future__ import annotations

import csv
import hashlib
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MANUSCRIPT_DIR = ROOT / "manuscript_ready" / "stage21_m18_remediation"
PACKAGE_DIR = ROOT / "submission_package" / "stage21_m18_remediation"
SUPPLEMENT_DIR = PACKAGE_DIR / "supplementary"
# Use the verified manuscript copy shipped in this release as the Pandoc
# reference document for Word styles, headers, page numbers, and line numbers.
REFERENCE_DOC = ROOT / "manuscript" / "scientific_reports_remediated_manuscript.docx"
PDF_EXPORTER = ROOT / "scripts" / "stage21_m17_assembly" / "export_m17_pdf.ps1"
PANDOC = Path(r"C:\ProgramData\anaconda3\Library\bin\pandoc.exe")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def to_docx(md: Path, docx: Path) -> None:
    subprocess.run(
        [
            str(PANDOC), str(md), "--from=gfm", "--to=docx", "--standalone",
            f"--reference-doc={REFERENCE_DOC}",
            f"--resource-path={SUPPLEMENT_DIR};{PACKAGE_DIR};{ROOT}",
            "--output", str(docx),
        ],
        cwd=ROOT,
        check=True,
    )


def to_pdf(docx: Path, pdf: Path) -> None:
    subprocess.run(
        [
            "powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File",
            str(PDF_EXPORTER), "-DocxPath", str(docx), "-PdfPath", str(pdf),
        ],
        cwd=ROOT,
        check=True,
    )


def build_supplement_with_figures() -> Path:
    base = (SUPPLEMENT_DIR / "supplementary_information.md").read_text(encoding="utf-8")
    inserts = """

## Supplementary Figure S1 image

![Complete three-study estimates](Supplementary_Figure_S1.png)

## Supplementary Figure S2 image

![Unsupervised PCA diagnostics](Supplementary_Figure_S2.png)

## Supplementary Figure S3 image

![Within-dataset sample correlations](Supplementary_Figure_S3.png)
"""
    output = SUPPLEMENT_DIR / "supplementary_information_with_figures.md"
    output.write_text(base + inserts, encoding="utf-8")
    return output


def main() -> None:
    manuscript_md = MANUSCRIPT_DIR / "scientific_reports_remediated_manuscript.md"
    manuscript_docx = MANUSCRIPT_DIR / "scientific_reports_remediated_manuscript.docx"
    manuscript_pdf = MANUSCRIPT_DIR / "scientific_reports_remediated_manuscript.pdf"
    to_docx(manuscript_md, manuscript_docx)
    to_pdf(manuscript_docx, manuscript_pdf)
    for path in (manuscript_docx, manuscript_pdf):
        target = PACKAGE_DIR / path.name
        target.write_bytes(path.read_bytes())

    supplement_md = build_supplement_with_figures()
    supplement_docx = SUPPLEMENT_DIR / "supplementary_information.docx"
    supplement_pdf = SUPPLEMENT_DIR / "supplementary_information.pdf"
    to_docx(supplement_md, supplement_docx)
    to_pdf(supplement_docx, supplement_pdf)

    cover_md = PACKAGE_DIR / "cover_letter_required_fields_template.md"
    cover_docx = PACKAGE_DIR / "cover_letter_required_fields_template.docx"
    cover_pdf = PACKAGE_DIR / "cover_letter_required_fields_template.pdf"
    to_docx(cover_md, cover_docx)
    to_pdf(cover_docx, cover_pdf)

    files = [path for path in PACKAGE_DIR.rglob("*") if path.is_file()]
    manifest_path = PACKAGE_DIR / "package_file_manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "size_bytes", "sha256"])
        writer.writeheader()
        for path in sorted(files):
            if path == manifest_path:
                continue
            writer.writerow(
                {
                    "path": path.relative_to(PACKAGE_DIR).as_posix(),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )
    print(f"manuscript_pdf={manuscript_pdf}")
    print(f"supplement_pdf={supplement_pdf}")


if __name__ == "__main__":
    main()
