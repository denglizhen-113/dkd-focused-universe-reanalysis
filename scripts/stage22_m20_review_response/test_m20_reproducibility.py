#!/usr/bin/env python3
"""Numerical regression checks for the M20 primary robustness analysis."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "tables" / "stage22_m20_review_response"
INPUT = ROOT / "data_processed" / "m20_primary_reproduction"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    manifest = pd.read_csv(INPUT / "input_manifest.csv")
    for row in manifest.itertuples(index=False):
        path = INPUT / row.file
        require(path.exists(), f"missing frozen input: {path}")
        require(sha256(path) == row.sha256, f"hash mismatch: {path.name}")

    coverage = pd.read_csv(OUT / "primary_pathway_measurement_coverage.csv")
    common = dict(zip(coverage["reactome_name"], coverage["common_three_cohort_members"]))
    expected = {
        "Complement cascade": 51, "Coagulation pathway": 52,
        "Cell surface interactions at the vascular wall": 125,
        "Chemokine receptors bind chemokines": 46,
        "Extracellular matrix organization": 255,
        "Cellular response to hypoxia": 47,
        "Signaling by TGF-beta Receptor Complex": 85,
    }
    require(common == expected, f"common-gene coverage changed: {common}")

    results = pd.read_csv(OUT / "primary_pathway_studentized_maxT_results.csv")
    require(len(results) == 42, "expected 7 pathways x 3 cohorts x 2 estimands")
    primary = results.loc[results["estimand"].eq("common_three_cohort_intersection")]
    require(len(primary) == 21, "expected 21 primary pathway rows")
    require((primary["studentized_maxT_fwer_p"] > 0).all(), "exact/Monte Carlo P values must not be zero")
    allocations = dict(zip(primary["cohort"], primary["permutation_allocations"]))
    require(allocations["GSE96804"] == 100000, "GSE96804 Monte Carlo allocation count changed")
    require(allocations["GSE30528"] == 497420, "GSE30528 exact allocation count changed")
    require(allocations["GSE104948_H7"] == 480700, "H7 exact allocation count changed")

    rep = pd.read_csv(OUT / "primary_pathway_operational_replication.csv")
    calls = set(rep.loc[rep["operational_two_of_three_call"], "reactome_name"])
    expected_calls = {
        "Complement cascade", "Cell surface interactions at the vascular wall",
        "Chemokine receptors bind chemokines", "Extracellular matrix organization",
    }
    require(calls == expected_calls, f"operational call set changed: {calls}")
    coag = rep.loc[rep["reactome_name"].eq("Coagulation pathway")].iloc[0]
    require(int(coag["sources_studentized_maxT_lt_0_05"]) == 0, "coagulation result changed")

    loo = pd.read_csv(OUT / "gene_meta_leave_one_source_out.csv")
    require(set(loo["omitted_source"]) == {"NONE", "GSE96804", "GSE30528", "GSE104948_H7"}, "leave-one-source-out strata incomplete")
    require(loo["gene_symbol"].nunique() == 783, "leave-one-source-out family must retain 783 genes")

    summary = json.loads((OUT / "m20_robustness_summary.json").read_text(encoding="utf-8"))
    require(set(summary["replicated_pathways"]) == expected_calls, "summary and table disagree")
    print("M20_NUMERICAL_REGRESSION=PASS")


if __name__ == "__main__":
    main()
