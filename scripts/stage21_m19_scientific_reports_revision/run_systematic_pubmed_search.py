#!/usr/bin/env python3
"""Reproducible PubMed accession-discovery search for the M19 DKD reanalysis."""

from __future__ import annotations

import csv
import json
import re
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "stage21_m19_scientific_reports_revision" / "systematic_search"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
COMMON = (
    '(transcriptom*[Title/Abstract] OR "gene expression"[Title/Abstract] OR '
    'microarray[Title/Abstract] OR RNA-seq[Title/Abstract] OR sequencing[Title/Abstract]) AND '
    '(kidney[Title/Abstract] OR renal[Title/Abstract] OR glomerul*[Title/Abstract] OR '
    'tubulointerstiti*[Title/Abstract]) AND (human*[Title/Abstract] OR patient*[Title/Abstract])'
)
QUERIES = {
    "PM_Q1": f'("diabetic nephropathy"[Title/Abstract]) AND {COMMON}',
    "PM_Q2": f'("diabetic kidney disease"[Title/Abstract]) AND {COMMON}',
}


def get(endpoint: str, params: dict[str, str | int]) -> bytes:
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(params)
    request = urllib.request.Request(url, headers={"User-Agent": "M19-DKD-systematic-search/1.0"})
    with urllib.request.urlopen(request, timeout=90) as response:
        return response.read()


def text_content(node: ET.Element | None) -> str:
    return "" if node is None else "".join(node.itertext()).strip()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pmids_by_query: dict[str, list[str]] = {}
    logs = []
    for query_id, query in QUERIES.items():
        root = ET.fromstring(
            get("esearch.fcgi", {"db": "pubmed", "term": query, "retmax": 100000, "retmode": "xml"})
        )
        pmids = [node.text or "" for node in root.findall(".//IdList/Id")]
        pmids_by_query[query_id] = pmids
        logs.append(
            {
                "query_id": query_id,
                "query": query,
                "search_date": str(date.today()),
                "count": int(root.findtext("Count", "0")),
                "retrieved": len(pmids),
                "query_translation": root.findtext("QueryTranslation", ""),
            }
        )

    all_pmids = sorted(set().union(*map(set, pmids_by_query.values())), key=int)
    articles: dict[str, dict] = {}
    for start in range(0, len(all_pmids), 200):
        batch = all_pmids[start : start + 200]
        root = ET.fromstring(
            get("efetch.fcgi", {"db": "pubmed", "id": ",".join(batch), "retmode": "xml"})
        )
        for article in root.findall(".//PubmedArticle"):
            pmid = article.findtext(".//MedlineCitation/PMID", "")
            title = text_content(article.find(".//ArticleTitle"))
            abstract = " ".join(text_content(node) for node in article.findall(".//Abstract/AbstractText"))
            journal = article.findtext(".//Journal/Title", "")
            year = article.findtext(".//PubDate/Year", "") or article.findtext(".//PubDate/MedlineDate", "")
            doi = ""
            for article_id in article.findall(".//ArticleId"):
                if article_id.attrib.get("IdType") == "doi":
                    doi = article_id.text or ""
            combined = f"{title} {abstract}"
            accessions = sorted(set(re.findall(r"\b(?:GSE|E-MTAB-|E-GEOD-)\d+\b", combined, flags=re.I)))
            articles[pmid] = {
                "pmid": pmid,
                "year": year,
                "journal": journal,
                "title": title,
                "doi": doi,
                "query_ids": ";".join(
                    query_id for query_id, ids in pmids_by_query.items() if pmid in set(ids)
                ),
                "repository_accessions_in_title_or_abstract": ";".join(accessions),
                "abstract": abstract,
            }

    rows = [articles[pmid] for pmid in all_pmids if pmid in articles]
    fieldnames = list(rows[0]) if rows else [
        "pmid", "year", "journal", "title", "doi", "query_ids",
        "repository_accessions_in_title_or_abstract", "abstract",
    ]
    with (OUT / "pubmed_records_all.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    # Resolve every GSE accession mentioned in titles/abstracts back to an
    # official GEO DataSets record. This makes the supplementary discovery arm
    # screenable even when the broad repository query did not retrieve it.
    discovered_accessions = sorted(
        {
            accession.upper()
            for row in rows
            for accession in row["repository_accessions_in_title_or_abstract"].split(";")
            if re.fullmatch(r"GSE\d+", accession, flags=re.I)
        }
    )
    geo_union_path = OUT / "geo_records_all.csv"
    geo_union = set()
    if geo_union_path.exists():
        with geo_union_path.open("r", encoding="utf-8-sig", newline="") as handle:
            geo_union = {row["accession"] for row in csv.DictReader(handle)}
    accession_query = " OR ".join(f"{accession}[ACCN]" for accession in discovered_accessions)
    search_root = ET.fromstring(
        get(
            "esearch.fcgi",
            {"db": "gds", "term": f"({accession_query}) AND gse[Entry Type]", "retmax": 10000, "retmode": "xml"},
        )
    )
    gds_uids = [node.text or "" for node in search_root.findall(".//IdList/Id")]
    discovered_geo = []
    for start in range(0, len(gds_uids), 100):
        payload = json.loads(
            get(
                "esummary.fcgi",
                {"db": "gds", "id": ",".join(gds_uids[start : start + 100]), "retmode": "json"},
            )
        )["result"]
        for uid in payload.get("uids", []):
            record = payload[uid]
            accession = record.get("accession", "")
            text = f"{record.get('title', '')} {record.get('summary', '')}".lower()
            discovered_geo.append(
                {
                    "accession": accession,
                    "title": record.get("title", ""),
                    "summary": record.get("summary", ""),
                    "organism": ";".join(record.get("taxon", [])) if isinstance(record.get("taxon"), list) else record.get("taxon", ""),
                    "gds_type": record.get("gdstype", ""),
                    "n_samples": record.get("n_samples", ""),
                    "already_in_primary_geo_union": accession in geo_union,
                    "priority_screen": bool(
                        re.search(r"human|homo sapiens", text)
                        and re.search(r"kidney|renal|glomerul|tubulointerstiti|biopsy", text)
                        and re.search(r"diabetic nephropathy|diabetic kidney disease|\bdkd\b", text)
                    ),
                }
            )
    discovered_geo.sort(key=lambda row: int(re.sub(r"\D", "", row["accession"])))
    discovered_fields = list(discovered_geo[0]) if discovered_geo else [
        "accession", "title", "summary", "organism", "gds_type", "n_samples",
        "already_in_primary_geo_union", "priority_screen",
    ]
    with (OUT / "pubmed_discovered_geo_records.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=discovered_fields)
        writer.writeheader()
        writer.writerows(discovered_geo)
    (OUT / "pubmed_search_log.json").write_text(
        json.dumps(
            {
                "database": "PubMed",
                "search_date": str(date.today()),
                "queries": logs,
                "unique_records": len(rows),
                "records_with_repository_accession_in_title_or_abstract": sum(
                    bool(row["repository_accessions_in_title_or_abstract"]) for row in rows
                ),
                "unique_gse_accessions_discovered": len(discovered_accessions),
                "gse_accessions_not_in_primary_geo_union": sum(
                    not row["already_in_primary_geo_union"] for row in discovered_geo
                ),
                "purpose": "Supplementary accession discovery; repository records remain the unit of eligibility.",
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    print(json.dumps({"unique_records": len(rows)}))


if __name__ == "__main__":
    main()
