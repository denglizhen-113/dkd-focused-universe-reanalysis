#!/usr/bin/env python3
"""Reproducible ArrayExpress/BioStudies search for the M19 DKD reanalysis."""

from __future__ import annotations

import csv
import json
import math
import re
import time
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "stage21_m19_scientific_reports_revision" / "systematic_search"
BASE = "https://www.ebi.ac.uk/biostudies/api/v1/search"
QUERIES = {
    "AE_Q1": '"diabetic nephropathy"',
    "AE_Q2": '"diabetic kidney disease"',
}


def request_json(params: dict[str, str | int]) -> dict:
    url = BASE + "?" + urllib.parse.urlencode(params)
    last_error: Exception | None = None
    for attempt in range(5):
        try:
            request = urllib.request.Request(
                url, headers={"User-Agent": "M19-DKD-systematic-search/1.0"}
            )
            with urllib.request.urlopen(request, timeout=90) as response:
                return json.load(response)
        except Exception as error:  # transient TLS resets occur at this endpoint
            last_error = error
            time.sleep(2**attempt)
    raise RuntimeError(f"BioStudies request failed after retries: {url}") from last_error


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    all_records: dict[str, dict] = {}
    log_rows = []
    for query_id, query in QUERIES.items():
        first = request_json(
            {"query": query, "collection": "ArrayExpress", "pageSize": 100, "page": 1}
        )
        pages = max(1, math.ceil(int(first["totalHits"]) / 100))
        hits = list(first.get("hits", []))
        for page in range(2, pages + 1):
            result = request_json(
                {"query": query, "collection": "ArrayExpress", "pageSize": 100, "page": page}
            )
            hits.extend(result.get("hits", []))
        log_rows.append(
            {
                "query_id": query_id,
                "query": query,
                "collection": "ArrayExpress",
                "search_date": str(date.today()),
                "total_hits": int(first["totalHits"]),
                "retrieved_hits": len(hits),
                "expanded_efo_terms": first.get("expandedEfoTerms", []),
                "expanded_synonyms": first.get("expandedSynonyms", []),
                "endpoint": BASE,
            }
        )
        for hit in hits:
            accession = hit.get("accession", "")
            if not accession:
                continue
            row = all_records.setdefault(
                accession,
                {
                    "accession": accession,
                    "title": hit.get("title", ""),
                    "author": hit.get("author", ""),
                    "release_date": hit.get("release_date", ""),
                    "content": hit.get("content", ""),
                    "query_ids": [],
                },
            )
            row["query_ids"].append(query_id)

    output_rows = []
    for accession, row in sorted(all_records.items()):
        text = f"{row['title']} {row['content']}".lower()
        human = "homo sapiens" in text or "human" in text
        renal = bool(re.search(r"kidney|renal|glomerul|tubul|interstiti|nephrectom|biopsy", text))
        disease = bool(re.search(r"diabetic nephropathy|diabetic kidney disease|diabetic renal|\bdkd\b", text))
        tissue = bool(re.search(r"kidney tissue|renal tissue|glomeruli|glomerular|tubulointerstiti|renal biopsy|kidney biopsy", text))
        experimental_only = bool(re.search(r"cell line|podocyte[s]? treated|cultured|mouse|mice|rat\b", text)) and not tissue
        geo_match = re.fullmatch(r"E-GEOD-(\d+)", accession)
        output_rows.append(
            {
                "accession": accession,
                "title": row["title"],
                "release_date": row["release_date"],
                "query_ids": ";".join(sorted(set(row["query_ids"]))),
                "human_signal": human,
                "renal_signal": renal,
                "disease_signal": disease,
                "kidney_tissue_signal": tissue,
                "experimental_only_signal": experimental_only,
                "priority_screen": human and renal and disease and not experimental_only,
                "geo_equivalent": f"GSE{geo_match.group(1)}" if geo_match else "",
                "database_duplicate_of_GEO": bool(geo_match),
                "content": row["content"],
            }
        )

    fieldnames = list(output_rows[0])
    with (OUT / "arrayexpress_records_all.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)
    priority = [row for row in output_rows if row["priority_screen"]]
    with (OUT / "arrayexpress_records_priority_screen.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(priority)
    (OUT / "arrayexpress_search_log.json").write_text(
        json.dumps(
            {
                "search_date": str(date.today()),
                "database": "BioStudies ArrayExpress collection",
                "queries": log_rows,
                "unique_records": len(output_rows),
                "priority_records": len(priority),
                "method_note": (
                    "All pages were retrieved. E-GEOD accessions were recorded as database duplicates "
                    "of their GSE equivalents before eligibility screening."
                ),
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    print(json.dumps({"unique_records": len(output_rows), "priority_records": len(priority)}))


if __name__ == "__main__":
    main()
