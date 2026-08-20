# Compartment-Stratified Systematic Reanalysis of Complement, Coagulation, and Matrix Transcriptional Programs in Diabetic Kidney Disease

Lizhen Deng¹*  
¹College of Life Science and Technology, Huazhong University of Science and Technology, Wuhan, Hubei, China  
*Correspondence: Lizhen Deng; 3070116993@qq.com; ORCID 0009-0003-2428-8176

## Abstract

Public diabetic kidney disease (DKD) transcriptomes differ in tissue compartment, platform, control source, and clinical annotation, complicating claims of molecular replication. We systematically searched GEO, ArrayExpress, and PubMed through 21 August 2026 and screened 322 unique dataset records. Eleven GEO Series representing nine source studies met eligibility criteria. The primary analysis was restricted to three independent glomerular sources and a fixed seven-pathway Reactome family comprising 783 human genes; other compartments were analyzed separately. Random-effects synthesis used restricted maximum likelihood and modified Hartung–Knapp inference, with Benjamini–Hochberg correction across all 783 genes. No gene met FDR<0.05. Joint two-sided sample-label permutations with maxT family-wise correction identified reproducible positive net effects for complement cascade, chemokine-receptor binding, and extracellular-matrix organization in at least two of three glomerular sources. Coagulation, hypoxia, vascular-wall interaction, and TGF-β signaling did not meet the replication rule. Sex/age adjustment, leave-one-sample/gene analyses, and probe-aggregation sensitivities qualified but did not overturn the principal conclusions. These results support compartment-bounded transcriptional convergence of selected inflammatory and matrix programs, not causal pathway activation or a replicated coagulation signature.

## Keywords

diabetic kidney disease; glomerulus; complement; extracellular matrix; transcriptomics; meta-analysis

## Introduction

Diabetic kidney disease (DKD) combines glomerular, tubular, vascular, inflammatory, and fibrotic injury [1]. Human kidney transcriptomic studies have reported recurrent immune and matrix-associated expression changes, but the available cohorts differ in renal compartment, disease stage, control tissue, assay, and preprocessing [2–5]. Pooling these settings without a defined biological estimand can produce a precise number whose biological meaning is unclear.

Earlier versions of this project used manually assembled gene groups and reused datasets across discovery and synthesis. The present analysis replaces those groups with a fixed, externally versioned Reactome family [14], conducts a reproducible repository and literature search following PRISMA reporting principles [15], counts source studies rather than accessions when donor or parent-study overlap is present, and prohibits cross-compartment pooling. We asked two bounded questions: whether individual genes in the canonical family show multiplicity-controlled evidence across three independent glomerular sources, and whether canonical pathway mean effects recur under within-study joint-label tests with family-wise error control.

## Results

### Systematic search and eligible evidence

The searches identified 263 GEO records, 54 ArrayExpress records, and 57 GSE accessions through 1,088 PubMed records used solely for accession discovery. After 52 database duplicates were removed, 322 unique dataset records were screened; 269 were excluded at title/summary screening and 53 underwent full accession and sample-metadata review. Forty-two were excluded, leaving 11 Series representing nine source studies (Fig. 1; Supplementary Table S1). GSE30528/GSE30529 are compartments of GSE30122; GSE104948/GSE104954 are compartments of ERCB H7. GSE99339/GSE99325 were excluded because their participant/source labels overlap GSE104948/GSE104954.

### Compartment-specific design

The primary glomerular synthesis contained GSE96804 (41 DKD/20 controls), GSE30528 (9/13), and the H7 stratum of GSE104948 (7/18). The H7 restriction avoided mixing batches. Tubulointerstitial, whole/cortical-kidney, interstitium-only, and small-sample cohorts were evaluated separately (Fig. 2; Table 1). No effect was pooled across compartments.

### Primary gene-level synthesis

The seven canonical pathways contained 783 unique GOA-filtered human symbols. Complete estimates from all three primary glomerular sources were available for 582 genes. REML random-effects models with modified Hartung–Knapp inference yielded no gene meeting BH FDR<0.05 across the 783-member family (0/783; Fig. 3; Supplementary Table S4). This null multiplicity-controlled result limits individual-gene claims; it does not establish absence of all expression differences.

### Canonical pathway replication

Each cohort used the mean unaligned Hedges' g over mapped genes. Case/control labels were reassigned jointly across samples, preserving the observed gene-correlation structure. Exact enumeration was used when there were no more than 10,000 allocations; otherwise 10,000 Monte Carlo permutations and plus-one P values were used. maxT controlled family-wise error across the fixed seven pathways. Complement cascade, chemokine receptors bind chemokines, and extracellular matrix organization were positive in all three primary sources and had maxT P<0.05 in two sources, satisfying the prespecified replication rule (Fig. 4; Table 2). Coagulation was mixed in direction and significant in no primary source after maxT correction. Vascular-wall interactions were positive in all three but significant in only one; hypoxia and TGF-β were directionally heterogeneous.

### Sensitivity analyses

Sex adjustment in GSE96804 did not change any pathway direction. Age-plus-sex adjustment in GSE166239 changed one of seven directions; sex adjustment in GSE163603 changed none, although all DKD samples there were male and residual confounding remains. Leave-one-sample-out analysis produced no direction changes in the three primary glomerular cohorts; eight changes occurred only in contextual cohorts. Highest-mean versus median-probe aggregation correlated strongly across canonical-union gene effects (Pearson 0.913 in GSE30528 and 0.901 in GSE30529); only the nonreplicated coagulation pathway changed direction in GSE30528. Correlation flags were retained as screening flags, not automatic sample exclusions, because all could reflect disease-associated global structure rather than technical failure.

### Table 1. Eligible Series and analysis roles

| Series/analysis | Source study | Compartment | DKD/control | Role |
| --- | --- | --- | --- | --- |
| GSE96804 | GSE96804 | glomerular | 41 / 20 | primary glomerular meta-analysis |
| GSE30528 | GSE30122 | glomerular | 9 / 13 | primary glomerular meta-analysis |
| GSE30529 | GSE30122 | tubulointerstitial | 10 / 12 | compartment-specific contextual analysis |
| GSE104948_H7 | ERCB_H7 | glomerular | 7 / 18 | primary glomerular meta-analysis |
| GSE104954_H7 | ERCB_H7 | tubulointerstitial | 7 / 18 | compartment-specific contextual analysis |
| GSE142025_advanced | GSE142025 | whole/cortical kidney | 21 / 9 | compartment-specific contextual analysis |
| GSE166239 | GSE166239 | whole/cortical kidney | 6 / 6 | compartment-specific contextual analysis |
| GSE163603_interstitium | GSE163603 | interstitium only | 6 / 9 | compartment-specific contextual analysis |
| GSE1009_donor_averaged | GSE1009 | glomerular | 2 / 2 | small-sample glomerular sensitivity |
| GSE111154 | GSE111154 | whole/cortical kidney | 4 / 4 | compartment-specific contextual analysis |
| GSE199838 | GSE199838 | whole/cortical kidney | 3 / 3 | compartment-specific contextual analysis |

### Table 2. Primary glomerular pathway replication

| Reactome pathway | Positive sources | Sources with maxT P<0.05 | Primary call |
| --- | --- | --- | --- |
| Complement cascade | 3 | 2 | Replicated |
| Coagulation pathway | 2 | 0 | Not replicated |
| Cell surface interactions at the vascular wall | 3 | 1 | Not replicated |
| Chemokine receptors bind chemokines | 3 | 2 | Replicated |
| Extracellular matrix organization | 3 | 2 | Replicated |
| Cellular response to hypoxia | 1 | 1 | Not replicated |
| Signaling by TGF-beta Receptor Complex | 1 | 0 | Not replicated |

## Discussion

The principal result is a deliberate separation of gene-level and pathway-level evidence. No member of the 783-gene canonical union survived the small-k, multiplicity-controlled glomerular meta-analysis. Nevertheless, three broad programs showed study-wise, family-wise-controlled positive mean effects in at least two independent glomerular sources. Aggregation can stabilize weak, distributed effects, but it does not identify which genes are causal, whether proteins or biochemical pathway activity change, or which cell type produces the signal.

Complement and matrix remodeling are biologically compatible with DKD inflammation and scarring [1–5]. The chemokine-receptor result similarly supports a distributed inflammatory transcript pattern. By contrast, the coagulation pathway did not replicate: its net direction differed across primary studies and none passed maxT correction. The title retains coagulation because it was a fixed target of the reanalysis, but the conclusion is explicitly negative for a replicated coagulation signature.

Compartment discipline materially changed the interpretation. Glomeruli, tubulointerstitium, whole biopsy/cortex, and isolated interstitium have different cellular compositions and control sources. GSE30528 and GSE30529 cannot be counted as independent evidence, and GSE104948/GSE104954 likewise represent paired compartments from one source collection. Separate contextual results are therefore descriptive rather than components of a universal DKD effect.

The study remains limited by archived observational data. Sample-level age, sex, renal function, medication, batch, and ancestry were incompletely recorded. Label permutation is exact or Monte Carlo-valid only under exchangeability of the archived labels and cannot remove confounding. Control tissues include living-donor and nephrectomy sources. Platforms and upstream preprocessing differ, and the analysis begins from processed matrices or normalized GEO values rather than uniform raw-file reprocessing. Only three independent sources inform the primary gene synthesis, making heterogeneity and interval estimates imprecise [9,10]. The pathway family was fixed for this revision but not prospectively preregistered. Systematic-search screening was performed by one reviewer. Finally, overlapping pathway membership means that a positive mean effect is not evidence of independent mechanisms.

## Methods

### Search strategy and eligibility

GEO DataSets was searched through NCBI E-utilities with two field-qualified human GSE queries covering diabetic nephropathy/DKD and broader diabetes-plus-kidney terminology. The ArrayExpress collection in BioStudies was searched for “diabetic nephropathy” and “diabetic kidney disease”; its synonym expansion returned the same 54 records. PubMed was searched with two title/abstract expressions combining DKD terms, transcriptomic terms, kidney terms, and human/patient terms. All query strings, translations, timestamps, returned identifiers, summaries, and screening decisions are archived. PubMed articles were used only to discover repository accessions, so article counts were not mixed with dataset-record counts. Eligible studies contained human bulk kidney-tissue expression data with separable DKD and non-DKD groups. Cell culture, animal-only, blood/urine, miRNA-only, single-cell/nucleus, spatial, noncase-control, derived duplicate, and unresolvable-overlap records were excluded. Small studies were retained as sensitivities rather than excluded by an arbitrary size threshold [15].

### Canonical pathway family

Seven Reactome pathways were fixed before the M19 calculations: Complement cascade (R-HSA-166658), Coagulation pathway (R-HSA-9769740), Cell surface interactions at the vascular wall (R-HSA-202733), Chemokine receptors bind chemokines (R-HSA-380108), Extracellular matrix organization (R-HSA-1474244), Cellular response to hypoxia (R-HSA-1234174), and Signaling by TGF-beta Receptor Complex (R-HSA-170834). The 21 June 2026 Reactome GMT was filtered to symbols in the 28 May 2026 GOA human GAF, yielding 783 unique genes. Source URLs, retrieval dates, SHA-256 hashes, and memberships are archived [14]. Historical manually curated sets have no confirmatory role.

### Expression processing and study effects

Expression matrices were analyzed within cohorts. Microarray probes were mapped to uppercase gene symbols; the primary label-independent aggregation selected the probe with highest all-sample mean, with median-probe aggregation as a sensitivity. Author count tables were transformed to log2(CPM+0.5). GSE1009 apparent a/b technical replicates were averaged at the donor-label level. Hedges' g represented DKD minus control standardized mean difference [8], and two-sided Welch P values were retained for within-study description.

### Gene meta-analysis

The primary estimand was the mean glomerular DKD-control standardized difference across GSE96804, GSE30528, and GSE104948 H7. Each gene used REML between-study variance and modified Hartung–Knapp t inference with k−1 degrees of freedom [9,10]. Nonestimable genes remained in the 783-member family with P=1 for BH correction [11]. No cross-compartment estimate was calculated.

### Pathway permutation and replication rule

For each pathway and cohort, the statistic was the arithmetic mean of observed gene-level Hedges' g without aligning signs to an outcome-derived direction. All gene effects were recomputed jointly after reallocating the fixed number of case labels. When the number of allocations was ≤10,000 every allocation was enumerated; otherwise 10,000 allocations were sampled with seed values derived from 20260821. Two-sided plus-one P values and maxT FWER were reported over seven pathways. A pathway was called replicated only if at least two of the three independent primary glomerular sources had maxT FWER P<0.05 with the same net direction.

### Robustness and quality control

Measured-covariate OLS sensitivities used disease plus sex for GSE96804 and GSE163603, and disease plus standardized age plus sex for GSE166239. These models were not generalized to datasets lacking comparable metadata. Leave-one-gene and leave-one-sample analyses assessed pathway sign stability. Sample QC used median Pearson correlation over the 1,000 most variable genes and flagged robust z<−3. Probe aggregation was compared in GSE30528/GSE30529. Flags were not exclusion rules without independent technical evidence.

### Reproducibility and AI assistance

Scripts, search logs, fixed pathway files, source tables, checksums, and an environment specification accompany the submission. OpenAI Codex was used under author supervision for code generation, statistical comparison, consistency checking, figures, tables, and language drafting. No generative image model was used; figures were generated from archived numerical tables. The author is responsible for verifying all code, results, references, and prose. AI was not used to generate or alter primary data and is not an author.

## Data availability

All input expression data are publicly available in GEO under GSE1009, GSE30528, GSE30529, GSE96804, GSE104948, GSE104954, GSE111154, GSE142025, GSE163603, GSE166239, and GSE199838. Accession-level eligibility, sample definitions, processed effects, and hashes are supplied in Supplementary Tables S1–S19. No new primary data were generated.

## Code availability

The exact M19 code, search logs, tables, manuscript, and submission assets are publicly archived at https://github.com/denglizhen-113/dkd-focused-universe-reanalysis/tree/v1.1.1 and are also supplied in Source_Code_M19.zip and Source_Data_M19.zip. The Git tag provides a versioned public snapshot; no DOI has been assigned.

## References

1. Alicic, R. Z., Rooney, M. T. & Tuttle, K. R. Diabetic kidney disease: challenges, progress, and possibilities. *Clin. J. Am. Soc. Nephrol.* **12**, 2032–2045 (2017). https://doi.org/10.2215/CJN.11491116
2. Woroniecka, K. I. et al. Transcriptome analysis of human diabetic kidney disease. *Diabetes* **60**, 2354–2369 (2011). https://doi.org/10.2337/db10-1181
3. Pan, Y. et al. Dissection of glomerular transcriptional profile in patients with diabetic nephropathy: SRGAP2a protects podocyte structure and function. *Diabetes* **67**, 717–730 (2018). https://doi.org/10.2337/db17-0755
4. Fan, Y. et al. Comparison of kidney transcriptomic profiles of early and advanced diabetic nephropathy reveals potential new mechanisms for disease progression. *Diabetes* **68**, 2301–2314 (2019). https://doi.org/10.2337/db19-0204
5. Nordbø, O. P. et al. Transcriptomic analysis reveals partial epithelial-mesenchymal transition and inflammation as common pathogenic mechanisms in hypertensive nephrosclerosis and type 2 diabetic nephropathy. *Physiol. Rep.* **11**, e15825 (2023). https://doi.org/10.14814/phy2.15825
6. Squair, J. W. et al. Confronting false discoveries in single-cell differential expression. *Nat. Commun.* **12**, 5692 (2021). https://doi.org/10.1038/s41467-021-25960-2
7. Lähnemann, D. et al. Eleven grand challenges in single-cell data science. *Genome Biol.* **21**, 31 (2020). https://doi.org/10.1186/s13059-020-1926-6
8. Hedges, L. V. Distribution theory for Glass's estimator of effect size and related estimators. *J. Educ. Stat.* **6**, 107–128 (1981). https://doi.org/10.2307/1164588
9. IntHout, J., Ioannidis, J. P. A. & Borm, G. F. The Hartung-Knapp-Sidik-Jonkman method for random effects meta-analysis. *BMC Med. Res. Methodol.* **14**, 25 (2014). https://doi.org/10.1186/1471-2288-14-25
10. Röver, C., Knapp, G. & Friede, T. Hartung-Knapp-Sidik-Jonkman approach and its modification for random-effects meta-analysis with few studies. *BMC Med. Res. Methodol.* **15**, 99 (2015). https://doi.org/10.1186/s12874-015-0091-1
11. Benjamini, Y. & Hochberg, Y. Controlling the false discovery rate. *J. R. Stat. Soc. B* **57**, 289–300 (1995). https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
12. Langfelder, P. & Horvath, S. WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics* **9**, 559 (2008). https://doi.org/10.1186/1471-2105-9-559
13. Ritchie, M. E. et al. limma powers differential expression analyses. *Nucleic Acids Res.* **43**, e47 (2015). https://doi.org/10.1093/nar/gkv007
14. Gillespie, M. et al. The Reactome pathway knowledgebase 2022. *Nucleic Acids Res.* **50**, D687–D692 (2022). https://doi.org/10.1093/nar/gkab1028
15. Page, M. J. et al. The PRISMA 2020 statement. *BMJ* **372**, n71 (2021). https://doi.org/10.1136/bmj.n71

## Acknowledgements

The author acknowledges the investigators who generated and publicly shared the reanalyzed datasets.

## Funding

The author received no specific funding for this work.

## Ethics statement

This secondary analysis used only publicly available, de-identified transcriptomic data and involved no new recruitment, intervention, specimen collection, participant contact, or access to identifiable private information. No new participant consent was sought. Ethics and consent procedures for original collection remain the responsibility of the source studies and are described in their repository records and publications. No animal data were analyzed.

## Author contributions

Lizhen Deng: conceptualization, methodology, software, formal analysis, investigation, data curation, visualization, writing—original draft, writing—review and editing, and project administration.

## Competing interests

The author declares no competing interests.

## Figure legends

**Figure 1. Systematic dataset identification and eligibility.** Dataset records, rather than PubMed articles, are the screening unit. Full decisions are in Supplementary Table S1.

**Figure 2. Compartment-specific study roles.** Bars show DKD and control samples. Compartments and accessions from the same source were not treated as independent effects.

**Figure 3. Primary glomerular gene synthesis.** The 15 lowest unadjusted modified Hartung–Knapp P values are displayed for visualization; none of 783 genes met BH FDR<0.05.

**Figure 4. Canonical pathway evidence in three primary glomerular sources.** Values are unaligned mean Hedges' g and maxT-FWER P. Asterisks indicate P<0.05. Replication required concordant significance in at least two sources.
