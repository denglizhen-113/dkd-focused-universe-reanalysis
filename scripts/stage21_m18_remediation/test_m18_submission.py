from __future__ import annotations

import re
import unittest
from pathlib import Path

import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[2]
PACKAGE_DIR = ROOT / "submission_package" / "stage21_m18_remediation"
MANUSCRIPT = ROOT / "manuscript_ready" / "stage21_m18_remediation" / "scientific_reports_remediated_manuscript.md"


class TestM18Submission(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.text = MANUSCRIPT.read_text(encoding="utf-8")
        cls.manifest = pd.read_csv(PACKAGE_DIR / "figure_table_manifest.csv").fillna("")

    def test_primary_claims_updated(self) -> None:
        self.assertIn("none met meta-analysis FDR<0.05", self.text)
        self.assertIn("94/141", self.text)
        self.assertIn("39/145", self.text)
        self.assertIn("8/9", self.text)

    def test_prohibited_internal_labels_absent(self) -> None:
        prohibited = [
            r"\[ledger:", r"\[M\d+:", r"\[evidence:", r"Stage 21", r"Phase 3",
            r"M17", r"strict empirical", r"Multi-Cohort", r"independently eligible",
            r"multiple independent validation cohorts", r"independent validation studies",
            r"[A-Za-z]:\\",
        ]
        for pattern in prohibited:
            self.assertIsNone(re.search(pattern, self.text, flags=re.IGNORECASE), pattern)

    def test_known_mojibake_absent(self) -> None:
        for token in ("I虏", "Nordb酶", "L盲hnemann", "Cs谩rdi", "�"):
            self.assertNotIn(token, self.text)

    def test_single_nucleus_is_not_primary_claim(self) -> None:
        title_abstract = self.text.split("## Keywords", 1)[0]
        self.assertNotIn("single-nucleus", title_abstract.lower())
        self.assertIn("Single-nucleus results were removed", self.text)

    def test_manifest_items_exist(self) -> None:
        self.assertEqual(len(self.manifest), 25)
        for row in self.manifest.itertuples(index=False):
            for field in ("file_png", "file_pdf"):
                value = getattr(row, field)
                if value:
                    self.assertTrue((PACKAGE_DIR / value).exists(), value)
            self.assertTrue(str(row.source_data).strip())
            self.assertTrue(str(row.cited_in).strip())

    def test_main_figure_dimensions(self) -> None:
        for number in range(1, 5):
            path = PACKAGE_DIR / "main_figures" / f"Figure_{number}.png"
            with Image.open(path) as image:
                width, height = image.size
            self.assertGreaterEqual(width, 2400)
            self.assertGreaterEqual(height, 1400)

    def test_release_blockers_explicit(self) -> None:
        required = (PACKAGE_DIR / "AUTHOR_AND_RELEASE_FIELDS_REQUIRED.md").read_text(encoding="utf-8")
        self.assertIn("Full author names", required)
        self.assertIn("Public code repository URL", required)
        self.assertIn("archival DOI", required)
        self.assertIn("Final human visual approval", required)


if __name__ == "__main__":
    unittest.main()
