# M20 reviewer-response matrix

Date: 21 August 2026  
Scope: internal response to the supplied Scientific Reports-style review report. This is not represented as correspondence from or to the journal.

| Review issue | M20 action | Evidence | Status |
| --- | --- | --- | --- |
| Observational-label exchangeability | GSE96804 labels are restricted within sex; GSE30528 is exactly enumerated at fixed case count; ERCB is restricted to H7 and exactly enumerated. Residual unavailable covariates are stated in every output and the manuscript. | Supplementary Table S6; Methods | Addressed within available metadata; residual confounding remains |
| Platform-dependent pathway membership | The primary estimand is the fixed intersection measurable across all three sources. Cohort-specific mapped sets are sensitivity estimands only. | Supplementary Tables S6 and S8 | Complete |
| Unstudentized raw maxT | Primary inference now centers and scales each pathway by its permutation-null mean and SD before two-sided maxT. Raw-mean maxT is retained only as sensitivity. | Supplementary Table S6 | Complete |
| Monte Carlo precision | GSE96804 uses 100,000 restricted allocations with plus-one P values, MC SE and Wilson intervals. Both smaller cohorts use exhaustive enumeration. | Supplementary Tables S6 and S9 | Complete |
| Pathway uncertainty | Figure 4 and S6 report 5,000-resample stratified-bootstrap 95% intervals. The legend states that these intervals are descriptive and not multiplicity-adjusted. | Figure 4; Supplementary Table S6 | Complete |
| Leave-one-source-out | All 783 family members are recomputed with each primary source omitted. Two-source fits are labelled diagnostic. | Supplementary Table S10 | Complete |
| Complete PRISMA reporting | The short checklist is replaced with every PRISMA 2020 item and sub-item, including 10b, 13b–f, 16b, 20a–d, 23a–d and 24a–c. Figure 1 is explicitly adapted to dataset-record reporting. | PRISMA_2020_checklist.pdf; Figure 1 | Complete |
| Independent screening audit | All 53 full records, all source-overlap decisions and a seeded 10% sample of early exclusions underwent a second-pass rule audit. It is explicitly labelled AI-assisted and not an independent human dual screen. | Supplementary Table S22 | Partially addressed; human confirmation remains required |
| Formal risk of bias/confounding | Nine source studies are assessed for control origin, disease stage, batch, clinical confounding, preprocessing, missing metadata, overlap and ethics provenance. No opaque numeric score is used. | Supplementary Table S21 | Complete |
| Self-contained code | The code ZIP now includes analysis, regression tests, frozen expression/design/pathway inputs, pinned requirements and a README. | Source_Code_M20.zip | Complete |
| Clean reproduction | The ZIP was extracted into a new directory and run without project-external helpers or downloads. Seven primary output files were byte-identical to the workspace results. | Clean-run log summarized in M20_FINAL_SELF_AUDIT.md | Complete |
| Ethics and consent provenance | Source-specific committee, identifier and consent statements are reported when located; unresolved items are marked not located/not reported. No non-compliance is inferred. | Supplementary Table S21; Methods | Partially addressed; author must verify unresolved sources and local institutional determination |
| Closely related literature | Direct comparison now includes Zhang et al. 2025 in *Scientific Reports*, Jiao et al. 2022 on glomerular C1q/C3, Li et al. 2022, Hojjati et al. 2023 and Abdalla et al. 2020. Novelty is restricted to the inferential and reproducibility architecture. | Introduction, Discussion; Supplementary Table S23 | Complete |
| Exact source/donor overlap | The seven repeated H7 labels DN901, DN910, DN914, DN916, DN932, DN941 and DN947 are mapped from GSE99339 GSM records to their GSE104948 counterparts. GSE99339 is therefore not a fourth independent replication source. | Results; Supplementary Table S28 | Complete |
| Historical source lineage | The GSE104948 GEO-declared older-CDF relationships with GSE47183 and GSE32591 are recorded without inventing unverified one-to-one mappings. | Supplementary Table S28 | Complete within GEO evidence boundary |
| Formal GEO data citations | Eleven GEO records are individually cited in references 20–30 and linked to related source articles where available. | References; Supplementary Table S25 | Complete |
| Reference sequence | References are numbered 1–30 in first-use sequence; unused legacy references were removed. | Manuscript | Complete |
| Confirmatory language | The manuscript says M19-defined but not prospectively preregistered; M20 is post-review. “Replicated” is replaced by “operational criterion met” where it could overstate evidence. | Abstract, Methods, Discussion; Supplementary Table S24 | Complete |
| Control heterogeneity in main table | Primary control origin, allocation restriction and unavailable covariates are shown in Table 1. | Main Table 1 | Complete |
| H7 decision history | The H7 rule and lack of provable outcome-blind prospective fixation are disclosed. | Results; Supplementary Table S24 | Complete |
| Figure 1 | Adapted to the PRISMA 2020 flow fields and explicitly identifies dataset records as the reporting unit. | Figure 1 | Complete |
| Figure 3 | Title now states that the 12 genes are selected by lowest unadjusted P and that none survives multiplicity correction. | Figure 3 | Complete |
| Figure 4 | Replaced the heatmap with a study-wise forest display showing effect intervals, source identity, significance by fill and the operational count. | Figure 4 | Complete |
| Supplementary Figure S2 accessibility | Direction change is encoded by both shape and color. | Supplementary Figure S2 | Complete |
| Plotting software | Python, NumPy, pandas, SciPy and Matplotlib versions are stated. | Methods: Software and reproducibility | Complete |
| Supplement first page | Full study title and author identity are present. | supplementary_information.pdf | Complete |
| AI disclosure | Use is described by task and bounded by human accountability; the text states that AI did not autonomously determine retained inferential conclusions. | Methods: Software and reproducibility | Complete |
| Search-cutoff date | GEO was retrieved at 2026-08-20 17:41 UTC, equal to 2026-08-21 01:41 in Asia/Shanghai; ArrayExpress and PubMed logs also state 2026-08-21. Manuscript, Fig. 1 and cover letter consistently use 21 August 2026. | Archived search logs; Figure 1; manuscript; cover letter | Complete; no future-dated cutoff in the current 21 August 2026 environment |
| Cover letter and portal fields | Cover letter contains explicit placeholders for suggested reviewers with institutional contacts, exclusions and any prior Editorial Board Member discussion. It is deliberately named `cover_letter_AUTHOR_COMPLETION_REQUIRED.pdf` because those facts cannot be inferred safely. | cover letter; submission checklist | Blocking author action required |
| DOI-bearing archive | No DOI is claimed. A versioned GitHub tag is supplied; DOI deposition remains recommended. | Code availability | Optional author action |

## Result changes that must not be obscured

- Complement: common-gene mean in GSE96804 is 0.104 (bootstrap 95% CI −0.049 to 0.266; studentized maxT P=0.699). Complement meets the operational rule only through GSE30528 and ERCB H7.
- Vascular-wall interaction: now meets the operational two-source rule under the common-gene/studentized analysis.
- Chemokine-receptor binding: meets the rule in two sources.
- Extracellular-matrix organization: meets the rule in all three sources.
- Coagulation: meets the rule in zero sources and remains directionally mixed.
- Gene family: 0/783 genes meet BH FDR<0.05.

## Residual submission blockers

1. A human must independently confirm the screening/source-overlap sample if the editor requires true dual screening.
2. The author/institution must decide and document whether this public de-identified secondary analysis requires an institutional determination or exemption letter.
3. Unresolved original-source ethics/consent fields should be checked against full articles/supplements or clarified with data generators if the editor requests documentation.
4. The three bracketed declarations in `cover_letter_AUTHOR_COMPLETION_REQUIRED.pdf` must be replaced with author-confirmed information; the filename should then be changed to `cover_letter.pdf`.
5. A Zenodo or other DOI-bearing archive is optional but advisable.
