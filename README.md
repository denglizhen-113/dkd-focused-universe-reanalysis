# Focused-Universe DKD Transcriptomic Reanalysis (M18)

This repository is the reproducibility package for the manuscript *Focused-Universe Reanalysis of Complement, Coagulation, and Extracellular-Matrix Transcriptional Programs in Diabetic Kidney Disease*.

## Release v1.0.0 contents

- A complete 145-gene REML random-effects analysis with modified Hartung-Knapp inference.
- Sample-label permutation testing for nine predefined pathways.
- Donor-disjoint GSE30529 sensitivity analyses.
- Main figure-generation and document-export scripts.
- Complete M18 result tables and supplementary data tables.
- The processed inputs required by the M18 workflow.

## Data provenance

The public source datasets are available from the NCBI Gene Expression Omnibus:

- GSE142025: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142025
- GSE96804: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96804
- GSE30528: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30528
- GSE30529: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30529

The repository includes processed matrices and cached GEO SOFT files used by the M18 reanalysis. The public GEO records remain the authoritative source for the underlying study data.

## Reproduction

Create the recorded Conda environment, then run the remediation script from the repository root:

```powershell
conda env create --file environment.yml
conda activate dkd-focused-universe-m18
python scripts/stage21_m18_remediation/run_m18_remediation.py
python scripts/stage21_m18_remediation/run_m18_qc.py
python scripts/stage21_m18_remediation/build_m18_submission.py
```

The primary tables are written under `tables/stage21_m18_remediation/`. The manuscript source and rendered deliverables are under `manuscript/`.

## Important scope statement

GSE30528 and GSE30529 are two compartments from source study GSE30122 and are not treated as separate independent source studies. The primary external gene-level analysis is restricted to GSE30528.

## License

Code in this repository is released under the MIT License. Public GEO-derived data remain subject to their source repositories' terms and attribution requirements.
