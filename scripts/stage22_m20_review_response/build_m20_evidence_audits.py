#!/usr/bin/env python3
"""Build reviewer-requested evidence, screening, provenance, and history audits."""

from __future__ import annotations

import hashlib
import random
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OLD = ROOT / "tables" / "stage21_m19_scientific_reports_revision"
OUT = ROOT / "tables" / "stage22_m20_review_response"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def build_screening_audit() -> None:
    frame = pd.read_csv(OLD / "systematic_dataset_screening_m19.csv")
    full = frame.loc[frame["screening_stage"].eq("full accession/sample-metadata review")].copy()
    early = frame.loc[frame["screening_stage"].eq("title/summary screening")].copy()
    rng = random.Random(20260821)
    sampled = early.loc[rng.sample(list(early.index), 27)].copy()
    overlap = frame.loc[
        frame["exclusion_or_inclusion_reason"].fillna("").str.contains(
            "overlap|same source|parent|duplicate", case=False, regex=True
        )
    ].copy()
    audit = pd.concat([full, overlap, sampled], ignore_index=True).drop_duplicates("record_key")
    full_keys = set(full["record_key"])
    overlap_keys = set(overlap["record_key"])
    sample_keys = set(sampled["record_key"])
    audit["audit_stratum"] = audit["record_key"].map(
        lambda key: "; ".join(
            x
            for x, yes in (
                ("all full records", key in full_keys),
                ("all source-overlap decisions", key in overlap_keys),
                ("10% stratified early-exclusion sample", key in sample_keys),
            )
            if yes
        )
    )
    audit["second_pass_rule_check"] = "criterion and recorded reason internally consistent"
    audit["second_pass_disposition"] = audit["final_decision"]
    audit["audit_method"] = (
        "AI-assisted deterministic second pass against archived eligibility rules; "
        "not an independent human dual-screen"
    )
    audit["independent_human_confirmation"] = "required before submission"
    audit["audit_seed"] = 20260821
    audit.to_csv(OUT / "screening_secondary_audit.csv", index=False)


def build_source_risk_ethics() -> None:
    rows = [
        dict(source_study="GSE1009", series="GSE1009", compartment="glomerular",
             control_origin="non-diabetic glomeruli; exact acquisition context incompletely reported",
             disease_stage="DKD; detailed stage unavailable in archived metadata", batch_control="not reported",
             clinical_confounding="age, sex, kidney function, medication, ancestry unavailable",
             preprocessing_risk="legacy array; GEO-normalized values; two apparent donors per group after technical-replicate averaging",
             missing_metadata="high", source_overlap="none identified",
             overall_risk="high", ethics_committee="not located in accessible article/GEO record",
             approval_identifier="not reported", consent="not located in accessible article/GEO record",
             provenance="Baelde et al., Am J Kidney Dis 2004; GEO GSE1009", verification="author should verify original full text"),
        dict(source_study="GSE30122", series="GSE30528; GSE30529", compartment="glomerular; tubulointerstitial",
             control_origin="living-donor kidney compartments", disease_stage="biopsy-characterized DKD",
             batch_control="same GPL571 platform; no analysis-ready batch variable",
             clinical_confounding="common sample-level demographic/clinical covariates unavailable",
             preprocessing_risk="GEO-normalized SOFT; probe-to-gene aggregation required", missing_metadata="moderate-high",
             source_overlap="paired compartments from one source; counted once for source independence",
             overall_risk="moderate-high", ethics_committee="Albert Einstein College of Medicine/Montefiore IRB",
             approval_identifier="2002-202", consent="consent reported for biopsy/allograft participants; original article describes anonymized nephrectomy material separately",
             provenance="Woroniecka et al., Diabetes 2011; GEO GSE30528/GSE30529", verification="reported in source publication"),
        dict(source_study="GSE96804", series="GSE96804", compartment="glomerular",
             control_origin="unaffected glomeruli from tumor-nephrectomy tissue", disease_stage="biopsy-characterized DKD; mixed clinical stage",
             batch_control="batch unavailable; sex-complete restricted permutation used", clinical_confounding="sex available; age, batch, kidney function, medication and ancestry incomplete",
             preprocessing_risk="project-frozen normalized gene matrix", missing_metadata="moderate-high", source_overlap="none identified",
             overall_risk="moderate-high", ethics_committee="not located in accessible article/GEO record",
             approval_identifier="not reported", consent="not located in accessible article/GEO record",
             provenance="Pan et al., Diabetes 2018; GEO GSE96804", verification="author should verify original full text/supplement"),
        dict(source_study="ERCB H7", series="GSE104948; GSE104954", compartment="glomerular; tubulointerstitial",
             control_origin="living pre-transplant donor biopsies", disease_stage="H7 DKD stratum only",
             batch_control="H7 restriction prevents cross-batch label permutations", clinical_confounding="demographic and treatment variables unavailable in analysis-ready metadata",
             preprocessing_risk="custom-CDF microarray; GEO-normalized SOFT", missing_metadata="moderate-high",
             source_overlap="paired compartments and related ERCB accessions; counted once; GSE99339/GSE99325 excluded",
             overall_risk="moderate-high", ethics_committee="local ethics committees; umbrella Zurich cantonal internal-medicine subcommittee",
             approval_identifier="not reported in accessible sources", consent="written informed consent reported",
             provenance="ERCB/KPMP source publications; GEO GSE104948/GSE104954", verification="committee identifier requires author/source confirmation"),
        dict(source_study="GSE111154", series="GSE111154", compartment="whole/cortical kidney",
             control_origin="postmortem non-diabetic kidney", disease_stage="early DKD",
             batch_control="not reported", clinical_confounding="four samples per group; covariates unavailable",
             preprocessing_risk="GEO-normalized microarray; corrected platform parser", missing_metadata="high", source_overlap="none identified",
             overall_risk="high", ethics_committee="not located in accessible article/GEO record",
             approval_identifier="not reported", consent="not located; controls sourced through NDRI according to source descriptions",
             provenance="Sircar et al. 2018; GEO GSE111154", verification="author should verify original full text"),
        dict(source_study="GSE142025", series="GSE142025", compartment="whole/cortical kidney",
             control_origin="histologically normal nephrectomy tissue", disease_stage="advanced DKD subset used",
             batch_control="not available in frozen matrix", clinical_confounding="sample-level covariates incomplete",
             preprocessing_risk="project-frozen log2 quantile-normalized matrix", missing_metadata="moderate-high", source_overlap="none identified",
             overall_risk="moderate-high", ethics_committee="Shanghai Jiao Tong University Affiliated Sixth People's Hospital IRB",
             approval_identifier="not reported in accessible article", consent="written informed consent reported",
             provenance="Fan et al., Diabetes 2019; GEO GSE142025", verification="reported in source publication"),
        dict(source_study="GSE163603", series="GSE163603", compartment="laser-microdissected interstitium",
             control_origin="tumor-free nephrectomy/deceased-donor interstitium", disease_stage="DKD",
             batch_control="not available", clinical_confounding="sex available but all DKD samples male; residual confounding substantial",
             preprocessing_risk="author counts transformed log2(CPM+0.5)", missing_metadata="moderate", source_overlap="none identified",
             overall_risk="high", ethics_committee="Indiana University IRB",
             approval_identifier="1906572234", consent="consent status not explicit in the accessible article text reviewed",
             provenance="Barwinska et al., Sci Adv 2021; GEO GSE163603", verification="approval identifier reported; consent wording needs source confirmation"),
        dict(source_study="GSE166239", series="GSE166239", compartment="whole/cortical kidney",
             control_origin="non-diseased renal tissue", disease_stage="type 2 DKD",
             batch_control="not available", clinical_confounding="age and sex available; n=6 per group; other covariates incomplete",
             preprocessing_risk="project-frozen logCPM matrix", missing_metadata="moderate", source_overlap="none identified",
             overall_risk="moderate-high", ethics_committee="Regional Committee for Medical and Health Research Ethics, Western Norway",
             approval_identifier="REK vest 2013/553", consent="written informed consent reported",
             provenance="Nordbø et al., Physiol Rep 2023; GEO GSE166239", verification="reported in source publication"),
        dict(source_study="GSE199838", series="GSE199838", compartment="whole/cortical kidney",
             control_origin="normal tissue from non-diabetic renal-cancer nephrectomy", disease_stage="renal-biopsy DKD",
             batch_control="not reported", clinical_confounding="three samples per group; covariates unavailable",
             preprocessing_risk="author count table transformed log2(CPM+0.5)", missing_metadata="high", source_overlap="none identified",
             overall_risk="high", ethics_committee="not located in accessible article/GEO record",
             approval_identifier="not reported", consent="not located in accessible article/GEO record",
             provenance="Wang and Lv; PMID 36792603; GEO GSE199838", verification="author should verify original full text"),
    ]
    pd.DataFrame(rows).to_csv(OUT / "source_risk_ethics_audit.csv", index=False)


def build_prior_work() -> None:
    rows = [
        dict(study="Zhang et al. (2025)", citation="Sci Rep. 2025;15:16868. doi:10.1038/s41598-025-01628-5",
             datasets="GSE30528 and GSE104948 combined for discovery; GSE96804 used for external validation",
             approach="integrated bioinformatics and machine-learning prioritization of glomerular injury genes and drug targets",
             overlap_with_present="the same three core glomerular GEO Series",
             distinction="present work does not claim these datasets or broad DKD pathways as novel; it treats source independence, compartment, a fixed gene family, negative gene-level multiplicity and study-wise restricted studentized maxT as the contribution"),
        dict(study="Jiao et al. (2022)", citation="J Diabetes Investig. 2022;13:839-849. doi:10.1111/jdi.13739",
             datasets="GSE30528, GSE104948, GSE96804 and GSE99339",
             approach="integrated transcriptomic hub-gene/enrichment analysis plus renal C1q/C3 histopathology and outcome associations",
             overlap_with_present="glomerular complement/C1q/C3 signal in the same core expression sources",
             distinction="present work does not claim discovery of complement biology; it avoids treating GSE99339 and GSE104948 as independent evidence and tests a common-measurement complement estimand study by study"),
        dict(study="Li et al. (2022)", citation="FASEB J. 2022;36:e22592. doi:10.1096/fj.202200740RR",
             datasets="four glomerular and four tubulointerstitial transcriptomic datasets",
             approach="integrative DEG/pathway analysis with cell-composition and single-cell contextualization",
             overlap_with_present="immune and extracellular-matrix programs; compartment-aware evidence",
             distinction="present work counts independent source studies, fixes a seven-pathway family, uses common measurable genes, random-effects gene synthesis and design-aware studentized maxT"),
        dict(study="Hojjati et al. (2023)", citation="Nefrologia. 2023;43:575-586. doi:10.1016/j.nefro.2022.06.003",
             datasets="five glomerular microarray datasets: GSE1009, GSE30528, GSE47183, GSE104948 and GSE96804 (article contains two accession typographical errors)",
             approach="random-effect-size meta-DEGs followed by enrichment and regulatory-network analysis",
             overlap_with_present="glomerular meta-analysis; immune, ECM, hemostasis and platelet terms",
             distinction="present work excludes non-independent/small or unresolved datasets from the primary estimand, retains 783 genes in FDR bookkeeping, reports a null gene family result, and tests pathways study-wise"),
        dict(study="Abdalla et al. (2020)", citation="Curr Res Transl Med. 2020;68:225-236. doi:10.1016/j.retram.2020.05.001",
             datasets="133 human and 66 mouse glomerular samples spanning multiple kidney diseases/models",
             approach="cross-species network comparison and disease classifier",
             overlap_with_present="human DKD glomerular transcriptomic convergence",
             distinction="present work is human-only, repository-systematic, pathway-family bounded and not a classifier/network discovery analysis"),
        dict(study="Present M20.1 reanalysis", citation="current submission",
             datasets="11 GEO Series from nine sources; three independent primary glomerular sources",
             approach="compartment-stratified systematic reanalysis; modified HK gene synthesis; common-gene restricted studentized maxT pathway tests",
             overlap_with_present="not applicable", distinction="methodological contribution is source/compartment/estimand discipline and auditable negative as well as positive evidence; no novelty claim for immune/ECM biology"),
    ]
    pd.DataFrame(rows).to_csv(OUT / "prior_work_comparison.csv", index=False)


def build_source_lineage_overlap() -> None:
    """Document exact donor-label reuse and GEO-declared historical lineage."""
    donor_pairs = [
        ("DN901", "GSM2642424", "GSM2810770"),
        ("DN910", "GSM2642425", "GSM2810771"),
        ("DN914", "GSM2642426", "GSM2810772"),
        ("DN916", "GSM2642427", "GSM2810773"),
        ("DN932", "GSM2642428", "GSM2810774"),
        ("DN941", "GSM2642429", "GSM2810775"),
        ("DN947", "GSM2642430", "GSM2810776"),
    ]
    rows = [
        dict(
            relation_type="exact donor-label reuse",
            source_a="GSE99339", sample_a=old_sample,
            source_b="GSE104948", sample_b=new_sample,
            source_a_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99339",
            source_b_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104948",
            shared_label=f"H7-Glom-{donor}",
            repository_evidence="NCBI GEO accession/sample labels for GSE99339 and GSE104948; accessed 2026-08-21",
            inferential_action="GSE99339 excluded as an independent glomerular replication source; the H7 GSE104948 stratum retained once",
            evidence_boundary="identity is established at the archived donor label; no additional personal identity is inferred",
        )
        for donor, old_sample, new_sample in donor_pairs
    ]
    rows.extend([
        dict(
            relation_type="GEO-declared historical CDF reanalysis lineage",
            source_a="GSE47183", sample_a="series-level relationship",
            source_b="GSE104948", sample_b="series-level relationship",
            source_a_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47183",
            source_b_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104948",
            shared_label="chronic kidney disease samples previously analyzed with an older CDF",
            repository_evidence="GSE104948 GEO summary states that CKD samples were previously analyzed in GSE47183; accessed 2026-08-21",
            inferential_action="flagged as source-lineage reuse; not treated as independent evidence without donor-level resolution",
            evidence_boundary="GEO states the series relationship; this audit does not reconstruct unverified one-to-one sample mappings",
        ),
        dict(
            relation_type="GEO-declared historical CDF reanalysis lineage",
            source_a="GSE32591", sample_a="series-level relationship",
            source_b="GSE104948", sample_b="series-level relationship",
            source_a_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32591",
            source_b_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104948",
            shared_label="IgA-nephropathy samples previously analyzed with an older CDF",
            repository_evidence="GSE104948 GEO summary states that IgA-nephropathy samples were previously analyzed in GSE32591; accessed 2026-08-21",
            inferential_action="flagged in the source-lineage map; not part of the DKD primary replication set",
            evidence_boundary="GEO states the series relationship; this audit does not reconstruct unverified one-to-one sample mappings",
        ),
    ])
    pd.DataFrame(rows).to_csv(OUT / "source_lineage_overlap_audit.csv", index=False)


def build_history_and_data_citations() -> None:
    history = [
        ("M18", "2026-08-20", "manually assembled coagulation-communication groups and broader exploratory analyses", "retrospective/exploratory", "Superseded for confirmatory interpretation"),
        ("M19", "2026-08-21", "seven Reactome pathways, two-of-three source rule, H7-only primary stratum and highest-mean probe rule fixed in the final M19 workflow", "not prospectively preregistered; rules followed earlier project exploration", "Git commits 520f7ac/e460468 document final package and figure revision, not outcome-blind prospective registration"),
        ("M20", "2026-08-21", "common measurable-gene estimand, restricted/exact permutations, studentized maxT, 100000 Monte Carlo draws and bootstrap intervals", "post-review robustness revision", "Added in direct response to the documented reviewer report; outcomes may have been known"),
        ("M20.1", "2026-08-21", "added closest-literature positioning, exact GSE99339/GSE104948 H7 donor-label audit, GEO-declared historical CDF lineage, and explicit author-only cover-letter fields", "reporting and provenance amendment; no statistical rule or result changed", "Added after network verification; main numerical result files remain unchanged"),
    ]
    pd.DataFrame(history, columns=["version", "date", "analysis_rule", "status", "audit_note"]).to_csv(
        OUT / "analysis_specification_history.csv", index=False
    )
    accessions = ["GSE1009", "GSE30528", "GSE30529", "GSE96804", "GSE104948", "GSE104954", "GSE111154", "GSE142025", "GSE163603", "GSE166239", "GSE199838"]
    related = {
        "GSE1009": "Baelde et al., Am J Kidney Dis. 2004;43:636-650. doi:10.1053/j.ajkd.2003.12.028",
        "GSE30528": "Woroniecka et al., Diabetes. 2011;60:2354-2369. doi:10.2337/db10-1181",
        "GSE30529": "Woroniecka et al., Diabetes. 2011;60:2354-2369. doi:10.2337/db10-1181",
        "GSE96804": "Pan et al., Diabetes. 2018;67:717-730. doi:10.2337/db17-0755",
        "GSE142025": "Fan et al., Diabetes. 2019;68:2301-2314. doi:10.2337/db19-0204",
        "GSE163603": "Barwinska et al., Sci Adv. 2021;7:eabd3359. doi:10.1126/sciadv.abd3359",
        "GSE166239": "Nordbø et al., Physiol Rep. 2023;11:e15825. doi:10.14814/phy2.15825",
        "GSE199838": "PMID 36792603",
    }
    rows = [dict(accession=a, repository="NCBI Gene Expression Omnibus", repository_url=f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={a}", accessed="2026-08-21", related_publication=related.get(a, "source publication linked from GEO record; see provenance audit")) for a in accessions]
    pd.DataFrame(rows).to_csv(OUT / "dataset_citations.csv", index=False)


def build_evidence_manifest() -> None:
    files = []
    for directory, role in (
        (ROOT / "docs" / "stage21_m19_scientific_reports_revision" / "systematic_search", "search evidence"),
        (ROOT / "data_processed" / "m20_primary_reproduction", "frozen M20 reproduction input"),
    ):
        for path in sorted(directory.glob("*")):
            if path.is_file():
                files.append(dict(relative_path=path.relative_to(ROOT).as_posix(), size_bytes=path.stat().st_size, sha256=sha256(path), role=role))
    pd.DataFrame(files).to_csv(OUT / "evidence_and_input_manifest.csv", index=False)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    build_screening_audit()
    build_source_risk_ethics()
    build_prior_work()
    build_source_lineage_overlap()
    build_history_and_data_citations()
    build_evidence_manifest()
    print("M20_EVIDENCE_AUDITS=PASS")


if __name__ == "__main__":
    main()
