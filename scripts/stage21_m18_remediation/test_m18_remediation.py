from __future__ import annotations

import json
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = ROOT / "tables" / "stage21_m18_remediation"
DOC_DIR = ROOT / "docs" / "stage21_m18_remediation"


class TestM18Remediation(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.focused = pd.read_csv(TABLE_DIR / "focused_universe_145.csv")
        cls.per_cohort = pd.read_csv(TABLE_DIR / "per_cohort_effect_sizes_145.csv")
        cls.primary = pd.read_csv(TABLE_DIR / "primary_reml_hk_145.csv")
        cls.external = pd.read_csv(TABLE_DIR / "external_gse30528_gene_validation_145.csv")
        cls.pathways = pd.read_csv(TABLE_DIR / "external_gse30528_pathway_permutation.csv")
        cls.tub = pd.read_csv(TABLE_DIR / "gse30529_donor_disjoint_gene_validation_145.csv")
        cls.tub_pathways = pd.read_csv(TABLE_DIR / "gse30529_donor_disjoint_pathway_permutation.csv")

    def test_complete_145_gene_family(self) -> None:
        self.assertEqual(len(self.focused), 145)
        self.assertEqual(self.focused["gene_symbol"].nunique(), 145)
        self.assertEqual(len(self.primary), 145)
        self.assertEqual(len(self.external), 145)

    def test_four_cohort_effect_table(self) -> None:
        self.assertEqual(len(self.per_cohort), 145 * 4)
        self.assertEqual(
            set(self.per_cohort["cohort"]),
            {"GSE142025_C1", "GSE96804", "GSE30528", "GSE30529"},
        )

    def test_primary_model_and_family(self) -> None:
        self.assertTrue(
            self.primary["model"].eq(
                "REML random effects with modified Hartung-Knapp primary t inference"
            ).all()
        )
        self.assertTrue(self.primary["fdr_bh_145_family"].between(0, 1).all())
        self.assertTrue(self.primary["p_value_for_145_family"].between(0, 1).all())
        incomplete = ~self.primary["complete_three_study_mapping"]
        self.assertTrue(self.primary.loc[incomplete, "p_value_for_145_family"].eq(1.0).all())

    def test_hk_and_prediction_intervals(self) -> None:
        complete = self.primary["complete_three_study_mapping"]
        frame = self.primary.loc[complete]
        self.assertTrue(np.isfinite(frame["ci_95_low_hk"]).all())
        self.assertTrue(np.isfinite(frame["ci_95_high_hk"]).all())
        self.assertTrue(np.isfinite(frame["prediction_interval_95_low"]).all())
        self.assertTrue(np.isfinite(frame["prediction_interval_95_high"]).all())
        self.assertTrue(np.isfinite(frame["prediction_interval_95_low_modified_hk"]).all())
        self.assertTrue(np.isfinite(frame["prediction_interval_95_high_modified_hk"]).all())
        self.assertTrue((frame["ci_95_low_hk"] <= frame["pooled_effect"]).all())
        self.assertTrue((frame["ci_95_high_hk"] >= frame["pooled_effect"]).all())

    def test_external_scope(self) -> None:
        self.assertTrue(
            self.external["validation_scope"].eq(
                "one external source study, prespecified glomerular compartment"
            ).all()
        )
        self.assertTrue(self.external["external_fdr_bh_145_family"].between(0, 1).all())

    def test_pathway_permutation(self) -> None:
        self.assertEqual(len(self.pathways), 9)
        self.assertTrue(self.pathways["permutations"].eq(10_000).all())
        self.assertTrue(self.pathways["one_sided_permutation_p"].between(0, 1).all())
        self.assertTrue(self.pathways["fdr_bh_9_pathways"].between(0, 1).all())
        self.assertTrue(
            self.pathways["correlation_handling"].str.contains("recomputed jointly").all()
        )

    def test_inputs_unchanged(self) -> None:
        hashes = pd.read_csv(TABLE_DIR / "input_hash_verification.csv")
        self.assertTrue(hashes["unchanged"].all())

    def test_donor_disjoint_secondary_scope(self) -> None:
        self.assertEqual(len(self.tub), 145)
        self.assertTrue(self.tub["n_case"].eq(5).all())
        self.assertTrue(self.tub["n_control"].eq(8).all())
        self.assertTrue(self.tub["analysis_scope"].str.contains("not an additional independent").all())
        self.assertEqual(len(self.tub_pathways), 9)
        self.assertTrue(self.tub_pathways["permutations"].eq(10_000).all())

    def test_diagnostics_and_outlier_sensitivity(self) -> None:
        sample_qc = pd.read_csv(TABLE_DIR / "sample_correlation_outlier_diagnostics.csv")
        sensitivity = pd.read_csv(TABLE_DIR / "primary_outlier_exclusion_sensitivity_145.csv")
        self.assertEqual(len(sensitivity), 145)
        self.assertTrue(sample_qc["predefined_outlier_flag"].isin([True, False]).all())
        self.assertEqual(int(sensitivity["fdr_lt_0_05"].sum()), 0)

    def test_decision_language(self) -> None:
        decision = (DOC_DIR / "M18_STATISTICAL_DECISION.md").read_text(encoding="utf-8")
        self.assertIn("analysis-defined focused universe", decision)
        self.assertIn("not a preregistered universe", decision)
        self.assertIn("Single-nucleus findings are not part of the primary", decision)
        self.assertNotIn("strict empirical", decision.lower())

    def test_settings(self) -> None:
        settings = json.loads((TABLE_DIR / "analysis_settings.json").read_text(encoding="utf-8"))
        self.assertEqual(settings["analysis_family"], 145)
        self.assertEqual(settings["primary_gse30122_estimand"], "GSE30528 glomerular compartment")
        self.assertEqual(settings["single_nucleus_role"], "supplementary provisional context only")


if __name__ == "__main__":
    unittest.main()
