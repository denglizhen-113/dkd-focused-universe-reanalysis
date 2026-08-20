#!/usr/bin/env python3
"""Reproducible NCBI GEO search for the M19 DKD transcriptomic review.

The script intentionally separates retrieval from final human eligibility review.
It queries the official NCBI E-utilities service, records the Entrez query
translation, retrieves all GSE document summaries, and assigns only a broad
priority-screen flag. Final inclusion and exclusion decisions belong in the
manually audited screening table.
"""

from __future__ import annotations

import csv
import json
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "stage21_m19_scientific_reports_revision" / "systematic_search"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "dkd_m19_prisma"
EMAIL = "3070116993@qq.com"

QUERIES = {
    "Q1_disease_phrases": (
        'Homo sapiens[Organism] AND gse[Entry Type] AND '
        '("diabetic nephropathy"[All Fields] OR '
        '"diabetic kidney disease"[All Fields] OR '
        '"diabetic glomerulosclerosis"[All Fields])'
    ),
    "Q2_broad_diabetes_kidney": (
        'Homo sapiens[Organism] AND gse[Entry Type] AND '
        'diabet*[All Fields] AND '
        '(kidney[All Fields] OR renal[All Fields] OR '
        'glomerul*[All Fields] OR tubulointersti*[All Fields])'
    ),
}

EXPRESSION_TYPES = {
    "Expression profiling by array",
    "Expression profiling by high throughput sequencing",
}


def request_json(endpoint: str, params: dict[str, str], retries: int = 4) -> dict:
    url = f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(params)}"
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": f"{TOOL}/1.0"})
            with urllib.request.urlopen(req, timeout=60) as response:
                return json.load(response)
        except Exception as exc:  # network/API retry boundary
            last_error = exc
            time.sleep(1.0 + attempt)
    raise RuntimeError(f"NCBI request failed after {retries} attempts: {url}") from last_error


def esearch(term: str) -> dict:
    payload = request_json(
        "esearch.fcgi",
        {
            "db": "gds",
            "term": term,
            "retmode": "json",
            "retmax": "10000",
            "tool": TOOL,
            "email": EMAIL,
        },
    )["esearchresult"]
    if int(payload["count"]) != len(payload["idlist"]):
        raise RuntimeError(
            f"ESearch returned {len(payload['idlist'])} IDs for count={payload['count']}"
        )
    return payload


def esummary(ids: list[str], chunk_size: int = 100) -> list[dict]:
    records: list[dict] = []
    for start in range(0, len(ids), chunk_size):
        chunk = ids[start : start + chunk_size]
        payload = request_json(
            "esummary.fcgi",
            {
                "db": "gds",
                "id": ",".join(chunk),
                "retmode": "json",
                "tool": TOOL,
                "email": EMAIL,
            },
        )["result"]
        for uid in payload.get("uids", []):
            records.append(payload[str(uid)])
        time.sleep(0.4)
    return records


def priority_screen(record: dict) -> tuple[bool, str]:
    text = f"{record.get('title', '')} {record.get('summary', '')}".lower()
    disease = any(
        term in text
        for term in (
            "diabetic nephropathy",
            "diabetic kidney disease",
            "diabetic glomerulosclerosis",
            "diabetic renal",
            "diabetic kidney",
        )
    )
    renal_tissue = any(
        term in text
        for term in (
            "kidney biopsy",
            "kidney tissue",
            "renal biopsy",
            "renal tissue",
            "glomeruli",
            "glomerular",
            "tubulointerstit",
        )
    )
    expression = record.get("gdstype", "") in EXPRESSION_TYPES
    flags = []
    if not disease:
        flags.append("disease wording not explicit")
    if not renal_tissue:
        flags.append("renal tissue wording not explicit")
    if not expression:
        flags.append("not a bulk-expression GEO type")
    return disease and renal_tissue and expression, "; ".join(flags) or "priority human review"


def flatten(record: dict, query_ids: list[str]) -> dict[str, object]:
    priority, reason = priority_screen(record)
    samples = record.get("samples", []) or []
    pubmed = record.get("pubmedids", []) or []
    return {
        "uid": record.get("uid", ""),
        "accession": record.get("accession", ""),
        "title": record.get("title", ""),
        "summary": record.get("summary", ""),
        "organism": record.get("taxon", ""),
        "entry_type": record.get("entrytype", ""),
        "gds_type": record.get("gdstype", ""),
        "platform_ids": record.get("gpl", ""),
        "publication_date": record.get("pdat", ""),
        "n_samples": record.get("n_samples", len(samples)),
        "pubmed_ids": ";".join(map(str, pubmed)),
        "bioproject": record.get("bioproject", ""),
        "ftp_link": record.get("ftplink", ""),
        "matched_queries": ";".join(sorted(query_ids)),
        "priority_screen": priority,
        "priority_reason": reason,
        "final_decision": "PENDING_HUMAN_REVIEW",
        "final_exclusion_reason": "",
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write an empty screening table: {path}")
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    retrieved_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()

    searches: dict[str, dict] = {}
    uid_queries: dict[str, list[str]] = {}
    for query_id, term in QUERIES.items():
        result = esearch(term)
        searches[query_id] = {
            "submitted_query": term,
            "count": int(result["count"]),
            "query_translation": result.get("querytranslation", ""),
            "uids": result["idlist"],
        }
        for uid in result["idlist"]:
            uid_queries.setdefault(uid, []).append(query_id)
        time.sleep(0.4)

    records = esummary(sorted(uid_queries, key=int))
    rows = [flatten(record, uid_queries[str(record["uid"])]) for record in records]
    rows.sort(key=lambda row: str(row["accession"]))

    search_log = {
        "retrieved_at_utc": retrieved_at,
        "database": "NCBI GEO DataSets (gds), restricted to GSE entry type",
        "service": "NCBI E-utilities ESearch and ESummary",
        "tool": TOOL,
        "queries": searches,
        "union_count": len(rows),
        "screening_note": (
            "Keyword priority_screen is not an eligibility decision. Every union record "
            "requires human title/summary review, and potentially full GEO/paper review."
        ),
    }
    (OUT / "geo_search_log.json").write_text(
        json.dumps(search_log, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    write_csv(OUT / "geo_records_all.csv", rows)
    write_csv(OUT / "geo_records_priority_screen.csv", [r for r in rows if r["priority_screen"]])

    print(f"Retrieved {len(rows)} unique GEO Series records")
    print(f"Priority-screen records: {sum(bool(r['priority_screen']) for r in rows)}")
    print(f"Outputs: {OUT}")


if __name__ == "__main__":
    main()
