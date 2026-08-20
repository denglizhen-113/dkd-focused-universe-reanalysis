# Compartment-Stratified DKD Transcriptomic Reanalysis (M20)

This repository is the reproducibility package for the manuscript *Compartment-Stratified Systematic Reanalysis of Complement, Coagulation, and Matrix Transcriptional Programs in Diabetic Kidney Disease*.

## Release v1.2.1 contents

The submission-cited network-verified snapshot is the [`v1.2.1` tag](https://github.com/denglizhen-113/dkd-focused-universe-reanalysis/tree/v1.2.1). The original M20 reviewer-response release remains under `v1.2.0`; historical M19 materials remain under `v1.1.0`–`v1.1.2`; M18 remains under `v1.0.0`–`v1.0.2`.

- Reproducible GEO, ArrayExpress/BioStudies, and PubMed accession-discovery searches.
- Dataset-level screening decisions for 322 unique records.
- A compartment-stratified primary analysis across three independent glomerular sources.
- Seven fixed Reactome pathways and a 783-gene canonical family.
- REML plus modified Hartung-Knapp gene synthesis with all 783 genes retained in FDR bookkeeping.
- A common-measurable-gene primary pathway estimand, sex-restricted or exact label allocations, studentized maxT family-wise inference, 100,000 Monte Carlo allocations where needed, and 5,000-resample bootstrap intervals.
- Leave-one-source-out gene synthesis, Monte Carlo precision reporting, source-level confounding/ethics provenance, a complete PRISMA 2020 checklist, and a transparent screening second pass.
- Explicit comparison with the closest 2025 *Scientific Reports* reuse and the 2022 complement/C1q/C3 study.
- Exact mapping of seven repeated H7 donor labels between GSE99339 and GSE104948, plus GEO-declared GSE47183/GSE32591 historical-CDF lineage boundaries.
- Publication-standard Figures 1–4, 28 supplementary tables, and a self-contained source-code ZIP validated from extraction under the locked scientific-computing environment.
- A deliberately marked cover-letter draft requiring author-supplied reviewer, exclusion and Editorial Board Member declarations before upload.

## Data provenance

The 11 included public Series are available from the NCBI Gene Expression Omnibus:

- GSE1009, GSE30528, GSE30529, GSE96804, GSE104948, GSE104954, GSE111154, GSE142025, GSE163603, GSE166239, and GSE199838.

The public repository records remain authoritative for primary data. Canonical pathway inputs are identified by official Reactome/GOA download URLs and SHA-256 hashes.

## Reproduction

The reviewer-requested primary robustness analysis is self-contained. Run from the repository root:

```powershell
python scripts/stage22_m20_review_response/run_m20_robustness.py
python scripts/stage22_m20_review_response/test_m20_reproducibility.py
```

The M20 tables are under `tables/stage22_m20_review_response/`; the manuscript, package, audit documents and upload-ready files are in the corresponding M20 directories. `submission_ready_scientific_reports_m20/Source_Code_M20.zip` contains the same analysis, frozen inputs and tests for extraction-only replay.

## Important scope statement

No kidney compartments are pooled. The primary synthesis is restricted to GSE96804, GSE30528 and GSE104948_H7. No individual gene met BH FDR<0.05 across 783 genes. Under the M20 common-measurement studentized maxT analysis, complement, vascular-wall interaction, chemokine-receptor binding and extracellular-matrix organization met an operational two-of-three-source criterion; coagulation did not. These are observational transcript associations conditional on the stated exchangeability restrictions, not causal or biochemical validation.

## License

Code in this repository is released under the MIT License. Public GEO-derived data remain subject to their source repositories' terms and attribution requirements.
