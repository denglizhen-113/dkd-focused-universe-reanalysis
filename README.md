# Compartment-Stratified DKD Transcriptomic Reanalysis (M19)

This repository is the reproducibility package for the manuscript *Compartment-Stratified Systematic Reanalysis of Complement, Coagulation, and Matrix Transcriptional Programs in Diabetic Kidney Disease*.

## Release v1.1.1 contents

The submission-cited M19 snapshot is the [`v1.1.1` tag](https://github.com/denglizhen-113/dkd-focused-universe-reanalysis/tree/v1.1.1). Tag `v1.1.0` is the initial M19 publication; `v1.1.1` updates the manuscript and cover letter to cite the public M19 snapshot. Historical M18 materials remain available under `v1.0.2`.

- Reproducible GEO, ArrayExpress/BioStudies, and PubMed accession-discovery searches.
- Dataset-level screening decisions for 322 unique records.
- A compartment-stratified primary analysis across three independent glomerular sources.
- Seven fixed Reactome pathways and a 783-gene canonical family.
- REML plus modified Hartung-Knapp gene synthesis and joint-label maxT pathway testing.
- Measured-covariate, leave-one-gene, leave-one-sample, correlation-QC, and probe-aggregation sensitivities.
- Complete M19 manuscript, figures, source tables, checksums, and upload-ready package.

## Data provenance

The 11 included public Series are available from the NCBI Gene Expression Omnibus:

- GSE1009, GSE30528, GSE30529, GSE96804, GSE104948, GSE104954, GSE111154, GSE142025, GSE163603, GSE166239, and GSE199838.

The public repository records remain authoritative for primary data. Canonical pathway inputs are identified by official Reactome/GOA download URLs and SHA-256 hashes.

## Reproduction

Run the M19 scripts from the repository root in the order below:

```powershell
python scripts/stage21_m19_scientific_reports_revision/run_systematic_geo_search.py
python scripts/stage21_m19_scientific_reports_revision/run_systematic_arrayexpress_search.py
python scripts/stage21_m19_scientific_reports_revision/run_systematic_pubmed_search.py
python scripts/stage21_m19_scientific_reports_revision/build_final_screening_table.py
python scripts/stage21_m19_scientific_reports_revision/run_m19_compartment_analysis.py
python scripts/stage21_m19_scientific_reports_revision/build_m19_submission_package.py
python scripts/stage21_m19_scientific_reports_revision/test_m19_submission_package.py
```

The primary tables are under `tables/stage21_m19_scientific_reports_revision/`; the manuscript, package, and upload-ready files are in the corresponding M19 directories.

## Important scope statement

No kidney compartments are pooled. The primary synthesis is restricted to GSE96804, GSE30528, and GSE104948_H7. No individual gene met BH FDR<0.05 across 783 genes. Complement, chemokine-receptor binding, and extracellular-matrix organization met the fixed pathway replication rule; coagulation did not.

## License

Code in this repository is released under the MIT License. Public GEO-derived data remain subject to their source repositories' terms and attribution requirements.
