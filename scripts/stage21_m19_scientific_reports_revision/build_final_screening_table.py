#!/usr/bin/env python3
"""Build an auditable dataset-level PRISMA screening table for M19."""

from __future__ import annotations

import csv
import json
import re
from collections import Counter
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SEARCH = ROOT / "docs" / "stage21_m19_scientific_reports_revision" / "systematic_search"
OUT = ROOT / "tables" / "stage21_m19_scientific_reports_revision"

INCLUDED = {
    "GSE96804": "primary glomerular source study",
    "GSE30528": "primary glomerular source study; GSE30122 glomerular subseries",
    "GSE104948": "primary glomerular source study; H7 batch only",
    "GSE30529": "tubulointerstitial contextual source study; GSE30122 subseries",
    "GSE104954": "tubulointerstitial contextual source study; H7 batch only",
    "GSE142025": "whole/cortical-kidney contextual analysis",
    "GSE166239": "whole/cortical-kidney contextual analysis",
    "GSE163603": "interstitium-only contextual analysis",
    "GSE1009": "small-sample glomerular sensitivity; apparent technical replicates donor-averaged",
    "GSE111154": "small-sample whole/cortical-kidney sensitivity",
    "GSE199838": "small-sample whole/cortical-kidney sensitivity",
}

FULL_REVIEW_EXCLUSIONS = {
    "GSE100185": "Deposited experiment does not provide an eligible independent human kidney-tissue DKD-versus-control bulk contrast.",
    "GSE117085": "Rat diabetic-kidney gastric-bypass experiment (Rattus norvegicus), not human kidney tissue.",
    "GSE129666": "Differentiated human podocyte culture with SMPDL3B perturbation, not patient kidney tissue.",
    "GSE131882": "Single-nucleus transcriptomics; outside the prespecified bulk-tissue estimand.",
    "GSE140308": "Growth-hormone-treated human podocyte culture, not patient kidney tissue.",
    "GSE154881": "Peripheral blood transcriptome, not kidney tissue.",
    "GSE158230": "Six-sample mechanistic/fibrosis series without an eligible non-DKD kidney comparator.",
    "GSE158626": "Experimental Set7 knockout/perturbation study; no eligible human kidney-tissue case-control contrast.",
    "GSE182138": "Hyperglycemic renal proximal-tubular cell culture with oxamate perturbation.",
    "GSE189875": "Glomerular endothelial cell NET perturbation experiment, not a patient tissue case-control cohort.",
    "GSE192889": "Mechanistic FGF13 endothelial/animal experiment without an eligible human tissue case-control contrast.",
    "GSE21785": "Living-donor glomerulus-versus-tubulointerstitium reference atlas; no DKD case group.",
    "GSE218929": "Urinary extracellular-vesicle transcriptome, not kidney tissue.",
    "GSE220226": "Glomerular cell co-culture exposed to methylglyoxal, not patient tissue.",
    "GSE220227": "Glomerular endothelial cell high-glucose co-culture, not patient tissue.",
    "GSE220228": "Podocyte high-glucose co-culture, not patient tissue.",
    "GSE226224": "Human podocyte culture experiment, not a DKD patient kidney cohort.",
    "GSE230848": "Early-DKD progression series without a non-DKD kidney comparator for the prespecified contrast.",
    "GSE255028": "Mechanistic KLK8 animal/endothelial study without an eligible human tissue case-control contrast.",
    "GSE261545": "Single-specimen spatial transcriptomics; outside the bulk-tissue estimand.",
    "GSE262793": "Deposited experiments are insulin-resistant kidney cell models; patient-biopsy comparisons reuse external cohorts and are not a new independent tissue series.",
    "GSE265918": "FGF9-treated renal tubular epithelial cell culture under high glucose.",
    "GSE301094": "Mechanistic GDF5 animal/cell experiment without an eligible human tissue case-control contrast.",
    "GSE30122": "SuperSeries represented by eligible compartment subseries GSE30528 and GSE30529; excluded to prevent double counting.",
    "GSE30566": "Control glomerulus-versus-control tubule contrast; no DKD case group.",
    "GSE32553": "Cultured podocyte glycosylation experiment, not patient kidney tissue.",
    "GSE45980": "Mixed-etiology CKD progression cohort; DKD subgroup lacks an eligible non-diseased comparator and adequate group sizes.",
    "GSE46900": "SuperSeries of serum-treated podocyte culture experiments; not patient kidney tissue.",
    "GSE47183": "Derived multi-disease composite glomerular dataset; independence from included source cohorts cannot be established.",
    "GSE47184": "Derived multi-disease composite tubulointerstitial dataset; independence from included source cohorts cannot be established.",
    "GSE47185": "SuperSeries of derived composite datasets GSE47183/GSE47184; excluded to prevent double counting.",
    "GSE51674": "miRNA-only profiling, outside the bulk mRNA transcriptome estimand.",
    "GSE60860": "Mixed-etiology CKD validation cohort; only one DKD sample and no eligible DKD-control contrast.",
    "GSE60861": "SuperSeries of mixed-etiology CKD progression cohorts; no eligible independent DKD-control contrast.",
    "GSE66494": "Mixed-etiology CKD cohort without a separable DKD-versus-non-diseased contrast.",
    "GSE69074": "Mechanistic RTN1 animal/cell series without an eligible human tissue case-control contrast.",
    "GSE69438": "Mixed-etiology CKD biomarker cohort with only one DN sample.",
    "GSE83144": "miRNA-25 animal/cell experiment, not an eligible human bulk kidney-tissue contrast.",
    "GSE99325": "Participant/source overlap with GSE104954; excluded to prevent duplicate tubulointerstitial representation.",
    "GSE99339": "Participant/source overlap with GSE104948; excluded to prevent duplicate glomerular representation.",
    "GSE99340": "Mixed SuperSeries containing duplicated renal-biopsy subseries and in-vitro podocyte experiments.",
    "E-MTAB-12135": "Peripheral-blood microarray (six DN, six controls), not kidney tissue.",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def title_screen_reason(title: str, summary: str) -> str:
    text = f"{title} {summary}".lower()
    if re.search(r"single[- ]cell|single[- ]nucleus|spatial transcript|atac-seq|methylation|methylome", text):
        return "Non-bulk or non-expression assay outside the prespecified bulk-transcriptome estimand."
    if re.search(r"urine|urinary|blood|serum|plasma|platelet|adipose|liver|retina|heart", text):
        return "Wrong biospecimen or organ; no eligible kidney-tissue bulk contrast."
    if re.search(r"cell line|cultured|culture|podocyte[s]? treated|high glucose|organoid|knockout|knockdown|overexpression", text):
        return "In-vitro or perturbational experiment without an eligible patient kidney-tissue contrast."
    if re.search(r"mouse|mice|rat\b|murine|zebrafish|animal model", text):
        return "Nonhuman or mixed experimental study without an eligible independent human tissue contrast."
    return "No eligible independent human bulk kidney-tissue DKD-versus-non-DKD contrast after title/summary screening."


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    geo = read_csv(SEARCH / "geo_records_all.csv")
    ae = read_csv(SEARCH / "arrayexpress_records_all.csv")
    discovered = read_csv(SEARCH / "pubmed_discovered_geo_records.csv")

    records: dict[str, dict[str, str]] = {}
    for row in geo:
        records[row["accession"]] = {
            "record_key": row["accession"],
            "accession": row["accession"],
            "record_source": "GEO direct search",
            "title": row["title"],
            "summary": row["summary"],
        }
    for row in ae:
        normalized_accession = row["geo_equivalent"] or row["accession"]
        if normalized_accession in records:
            continue
        records[normalized_accession] = {
            "record_key": normalized_accession,
            "accession": normalized_accession,
            "record_source": (
                "ArrayExpress E-GEOD record not retrieved by direct GEO query"
                if row["geo_equivalent"]
                else "ArrayExpress unique record"
            ),
            "title": row["title"],
            "summary": row["content"],
        }
    for row in discovered:
        if row["accession"] in records:
            continue
        records[row["accession"]] = {
            "record_key": row["accession"],
            "accession": row["accession"],
            "record_source": "PubMed accession discovery",
            "title": row["title"],
            "summary": row["summary"],
        }

    full_candidates = set(INCLUDED) | set(FULL_REVIEW_EXCLUSIONS)
    missing = full_candidates - set(records)
    if missing:
        raise AssertionError(f"Full-review records missing from search corpus: {sorted(missing)}")

    rows = []
    for accession, record in sorted(records.items()):
        if accession in INCLUDED:
            stage = "full accession/sample-metadata review"
            decision = "INCLUDED"
            reason = "Eligible human bulk kidney-tissue DKD-versus-non-DKD contrast."
            role = INCLUDED[accession]
        elif accession in FULL_REVIEW_EXCLUSIONS:
            stage = "full accession/sample-metadata review"
            decision = "EXCLUDED_FULL_RECORD"
            reason = FULL_REVIEW_EXCLUSIONS[accession]
            role = ""
        else:
            stage = "title/summary screening"
            decision = "EXCLUDED_TITLE_SUMMARY"
            reason = title_screen_reason(record["title"], record["summary"])
            role = ""
        rows.append(
            {
                "record_key": accession,
                "accession": accession,
                "record_source": record["record_source"],
                "title": record["title"],
                "screening_stage": stage,
                "final_decision": decision,
                "exclusion_or_inclusion_reason": reason,
                "analysis_role_if_included": role,
                "screening_date": str(date.today()),
            }
        )

    if len(rows) != 322:
        raise AssertionError(f"Expected 322 unique dataset records, found {len(rows)}")
    counts = Counter(row["final_decision"] for row in rows)
    expected = {"EXCLUDED_TITLE_SUMMARY": 269, "EXCLUDED_FULL_RECORD": 42, "INCLUDED": 11}
    if dict(counts) != expected:
        raise AssertionError(f"Unexpected PRISMA counts: {counts}")

    with (OUT / "systematic_dataset_screening_m19.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    full_reason_counts = Counter(
        row["exclusion_or_inclusion_reason"]
        for row in rows
        if row["final_decision"] == "EXCLUDED_FULL_RECORD"
    )
    prisma = {
        "search_date": str(date.today()),
        "records_identified": {
            "GEO_direct": 263,
            "ArrayExpress": 54,
            "PubMed_discovered_GSE_accessions": 57,
            "total_before_deduplication": 374,
        },
        "duplicates_removed": {
            "ArrayExpress_E_GEOD_duplicates_already_in_GEO_union": 24,
            "PubMed_discovered_GSE_already_in_GEO_union": 28,
            "total": 52,
        },
        "unique_dataset_records_screened": 322,
        "title_summary_excluded": 269,
        "full_accession_records_assessed": 53,
        "full_accession_records_excluded": 42,
        "series_included": 11,
        "independent_source_studies_included": 9,
        "full_record_exclusion_reasons": full_reason_counts,
        "unit_note": (
            "PubMed article records (n=1,088) were searched only to discover repository accessions. "
            "They are not added to the dataset-record PRISMA denominator."
        ),
    }
    (OUT / "prisma_counts_m19.json").write_text(
        json.dumps(prisma, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps({"unique": len(rows), **counts}, ensure_ascii=False))


if __name__ == "__main__":
    main()
