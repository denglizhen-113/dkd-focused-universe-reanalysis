# Compartment-Stratified Systematic Reanalysis of Complement, Coagulation, and Matrix Transcriptional Programs in Diabetic Kidney Disease

Lizhen Deng¹*  
¹College of Life Science and Technology, Huazhong University of Science and Technology, Wuhan, Hubei, China  
*Correspondence: Lizhen Deng; 3070116993@qq.com; ORCID 0009-0003-2428-8176

## Abstract

Public diabetic kidney disease (DKD) transcriptomes differ in tissue compartment, control source, platform, and clinical annotation. We systematically screened 322 unique dataset records and included 11 GEO Series from nine source studies. The primary analysis comprised three independent glomerular sources and seven fixed Reactome pathways (783 unique genes); other compartments were analyzed separately. No gene met Benjamini–Hochberg FDR<0.05 under restricted-maximum-likelihood synthesis with modified Hartung–Knapp inference. Pathway tests used genes measurable in all three cohorts, cohort-specific restricted or exact label allocations, and studentized maxT family-wise correction. Complement cascade, vascular-wall interaction, chemokine-receptor binding, and extracellular-matrix organization met an operational two-of-three-source criterion; coagulation met it in zero sources. Complement was weak in GSE96804 after the common-gene restriction. Bootstrap intervals and leave-one-source-out analyses exposed substantial uncertainty. These observational transcript associations remain conditional on archived-label exchangeability and cannot establish biochemical activation, causality, or a universal DKD signature. The contribution is a source-aware, compartment-bounded and cleanly reproducible assessment that reports negative multiplicity-controlled gene evidence alongside qualified pathway-level associations.

## Keywords

diabetic kidney disease; glomerulus; complement; extracellular matrix; transcriptomics; meta-analysis

## Introduction

Diabetic kidney disease (DKD) comprises glomerular, tubular, vascular, inflammatory and fibrotic injury, while public human-kidney transcriptomes vary markedly in disease stage, tissue acquisition and control origin [1]. These differences make a single cross-study molecular effect difficult to define.

Closely related analyses already exist. Li et al. integrated four glomerular and four tubular datasets with single-cell contextualization and reported considerable cross-dataset inconsistency before prioritizing TEKT2 and PIAS2 [2]. Hojjati et al. combined five glomerular microarray datasets by random-effects meta-analysis and reported 1,364 meta-DEGs followed by immune, extracellular-matrix and hemostasis enrichment [3]. Abdalla et al. compared human and mouse glomerular networks across several kidney diseases [4]. Zhang et al. used GSE30528 and GSE104948 for discovery and GSE96804 for external validation in a 2025 *Scientific Reports* machine-learning analysis [5]. Jiao et al. integrated GSE30528, GSE104948, GSE96804 and GSE99339 and specifically reported glomerular C1q/C3 findings [6]. Thus, neither reuse of the three core datasets nor complement, immune or matrix biology is claimed here as novel. The intended increment is methodological: dataset-level systematic identification, explicit source de-duplication, compartment separation, a fixed Reactome family [7], complete multiplicity bookkeeping, a common-measurement pathway estimand and design-aware study-wise inference.

Following PRISMA 2020 reporting principles [8], we asked two bounded questions. First, do any genes from the seven-pathway family show multiplicity-controlled evidence across three independent glomerular sources? Second, do pathway-average transcriptional associations meet an operational two-of-three-source rule when every source uses the same measurable members and cohort-appropriate label allocations? The analysis was not prospectively preregistered; version history and post-review amendments are reported explicitly.

## Results

### Dataset identification and evidence architecture

The searches identified 263 GEO records, 54 ArrayExpress records and 57 GSE accessions through 1,088 PubMed records used only for accession discovery. After removing 52 duplicates, 322 unique dataset records were screened; 269 were excluded at title/summary stage and 53 were assessed using accession and sample metadata. Forty-two were excluded, leaving 11 GEO Series from nine source studies (Fig. 1; Supplementary Table S1). The reporting unit in Fig. 1 is a *dataset record*, an adaptation of the PRISMA flow template. GSE30528/GSE30529 are paired compartments of GSE30122; GSE104948/GSE104954 are paired compartments of ERCB H7. GSE99339 contains the same archived H7 glomerular donor labels DN901, DN910, DN914, DN916, DN932, DN941 and DN947 that reappear in GSE104948, so it was not counted as a fourth independent replication source. The GSE104948 record also declares historical older-CDF analyses in GSE47183 and GSE32591; these series-level relationships were treated as source-lineage warnings rather than independent evidence (Supplementary Table S28). Related or overlapping accessions were not counted as independent sources.

The primary synthesis comprised GSE96804 (41 DKD/20 controls), GSE30528 (9/13) and the H7 stratum of GSE104948 (7/18). H7 was retained as a within-batch stratum, but this choice was not prospectively blinded to all earlier project results (Supplementary Table S24). Tubulointerstitial, whole/cortical-kidney, interstitium-only and very small glomerular cohorts remained contextual or sensitivity analyses (Fig. 2). Control origins included living-donor and tumor-nephrectomy tissue, and clinical covariates were incomplete (Table 1; Supplementary Tables S2 and S21).

### Gene-level synthesis and source sensitivity

The fixed seven-pathway family contained 783 unique GOA-filtered human symbols; 582 had estimates in all three primary sources. Restricted-maximum-likelihood random-effects models with modified Hartung–Knapp inference yielded 0/783 genes at Benjamini–Hochberg FDR<0.05 (Fig. 3; Supplementary Table S4). The 12 visually displayed lowest unadjusted P values are descriptive selections, not discoveries. Leave-one-source-out results are reported for all 783 family members (Supplementary Table S10); with only two retained sources per omission, these fits diagnose source dependence rather than provide stable alternative inference.

### Common-measurement, design-aware pathway tests

The three-source common measurable intersections contained 51 complement, 52 coagulation, 125 vascular-wall, 46 chemokine-receptor, 255 extracellular-matrix, 47 hypoxia and 85 TGF-β-pathway genes (Supplementary Table S8). These fixed intersections defined the primary pathway estimand; cohort-specific mapped sets were retained only as sensitivity analyses.

GSE96804 used 100,000 Monte Carlo allocations restricted within archived sex strata. GSE30528 used all 497,420 case-count-preserving allocations. GSE104948 H7 used all 480,700 allocations within the H7 stratum. The seven pathway statistics were standardized by their permutation-null mean and standard deviation before two-sided maxT correction. Complement, vascular-wall interaction, chemokine-receptor binding and extracellular-matrix organization met the operational criterion of studentized maxT P<0.05 with concordant positive direction in at least two sources (Fig. 4; Table 2; Supplementary Tables S6–S9). Extracellular matrix met the threshold in three sources; each of the other three met it in two. Coagulation met it in zero sources and was directionally mixed. Hypoxia met it in one source and TGF-β in none.

The common-gene restriction materially qualified the former complement presentation: GSE96804 complement mean Hedges' g was 0.104 (bootstrap 95% CI −0.049 to 0.266; studentized maxT P=0.699), whereas GSE30528 and H7 remained positive (0.482 and 0.483; maxT P=0.000623 and 0.000907). Thus, the operational complement call rests on two sources, not uniform strength across all three. Bootstrap intervals are descriptive within-cohort uncertainty intervals and are not multiplicity-adjusted.

### Robustness, bias and certainty boundaries

Probe aggregation, measured-covariate models, leave-one-gene/sample analyses and sample-correlation flags remain reported as sensitivities (Supplementary Tables S11–S17). They cannot repair unrecorded age, kidney function, medication, ancestry or batch variables. The formal source audit judged the three primary sources at moderate-high risk of confounding because control acquisition and clinical metadata differ; GSE1009, GSE111154, GSE163603 and GSE199838 were high risk in at least one domain (Supplementary Table S21). Reporting bias and a formal certainty grade were not estimated because only three heterogeneous primary sources were available and the dataset-search setting does not support conventional funnel-plot diagnostics.

### Table 1. Primary glomerular sources, controls and design restrictions

{{TABLE1}}

### Table 2. Operational primary pathway criterion

{{TABLE2}}

## Discussion

The most secure result is negative and family-bounded: no individual gene met the stated multiplicity threshold. Four broad pathway averages nevertheless met a study-wise operational criterion. This pattern is compatible with distributed shifts that are individually weak, but it does not identify causal genes, cell types, protein abundance or biochemical pathway activation.

The common-measurement analysis clarifies why pathway claims require more restraint than a heatmap of cohort-specific mapped sets suggests. Complement coverage varied from 59 to 105 measurable members, whereas only 51 were shared. Restricting the estimand weakened complement in the largest cohort, although two smaller sources retained positive, family-wise-significant associations. Vascular-wall interaction newly met the two-source rule under studentization and common-gene restriction. These changes are reported rather than concealed because they demonstrate sensitivity to the estimand and null scaling.

The results overlap prior DKD integrations but answer a different question. Li et al. emphasized cross-dataset inconsistency and cell-type context [2]; Hojjati et al. sought pooled meta-DEGs and subsequent enrichment [3]; Abdalla et al. sought cross-species network commonality [4]. Zhang et al. used the same three core GEO Series for feature and drug-target prioritization [5], while Jiao et al. reported glomerular complement/C1q/C3 findings after integrating those Series with the donor-overlapping GSE99339 record [6]. Here, source studies rather than accessions are the independence unit, exact repeated donor labels are exposed, compartments are never pooled, nonestimable genes stay in the FDR family, and pathway evidence is required to recur study by study. The novelty is therefore an auditable inferential architecture—not dataset novelty or discovery of previously unknown complement, inflammatory or matrix biology (Supplementary Tables S23 and S28).

Several limitations remain irreducible. Archived labels are exchangeable only conditionally on recorded design restrictions. Sex-stratified permutation addresses one variable in GSE96804, and the H7 restriction addresses a known batch boundary, but unmeasured clinical and acquisition confounding remains. Living-donor and tumor-nephrectomy controls are not biologically interchangeable. Three primary sources give imprecise heterogeneity estimates. Starting from processed matrices prevents uniform raw-data normalization. Pathway overlap creates dependent hypotheses even though maxT controls their joint family. The pathway family and two-of-three rule were fixed in M19 but were neither prospectively registered nor independent of earlier exploratory work; the common-gene/studentized analysis is explicitly post-review. Screening was performed by one author. An AI-assisted deterministic second pass covered all 53 full records, all recorded overlap decisions and a seeded 10% sample of early exclusions, but it is not falsely presented as an independent human dual screen; human confirmation remains an author task (Supplementary Table S22).

In conclusion, complement, vascular-wall interaction, chemokine-receptor binding and extracellular-matrix organization met the stated study-wise criterion across three independent glomerular sources, while coagulation did not. These are compartment-specific transcriptional associations conditional on common measurement and archived-label exchangeability—not causal mechanisms, protein-level activation or a universal DKD signature. The principal value of the analysis is the transparent boundary it places around both positive pathway evidence and negative gene-level evidence.

## Methods

### Search, eligibility and selection

GEO DataSets was searched through NCBI E-utilities with two field-qualified human GSE queries covering diabetic nephropathy/DKD and broader diabetes-plus-kidney terminology. The ArrayExpress collection in BioStudies was searched for “diabetic nephropathy” and “diabetic kidney disease”. PubMed used two title/abstract expressions combining DKD, transcriptomic, kidney and human/patient terms. Searches ended 21 August 2026. Exact queries, translations, timestamps, identifiers, summaries and hashes are archived. PubMed records were used only to discover accessions, so article and dataset-record counts were not mixed.

Eligible records contained human bulk kidney-tissue expression data with separable DKD and non-DKD groups. Cell culture, animal-only, blood/urine, miRNA-only, single-cell/nucleus, spatial, non-case-control, derived duplicate and unresolvable-overlap records were excluded. Small studies were assigned to sensitivity analyses rather than removed by an arbitrary size rule. One author performed selection. The structured second pass described above checks consistency with the archived rules but does not satisfy independent human dual screening. The complete PRISMA 2020 checklist and every full-record exclusion are supplied.

### Risk of bias, controls and ethics provenance

Each source was assessed for control selection, disease stage, batch handling, clinical confounding, preprocessing, missing metadata and source overlap. Judgements were qualitative because no validated risk-of-bias instrument is designed for retrospective public renal-transcriptome reanalysis. The domain-level rationale is reported rather than collapsed into a numeric score (Supplementary Table S21).

This work used only public, de-identified expression data and involved no recruitment, intervention, participant contact or access to identifiable private information. It does not claim a new ethics approval or waiver on behalf of the author's institution. Original-study committee, approval and consent statements were traced where accessible; unresolved fields are marked “not located” or “not reported”, not inferred. Ethics provenance was reviewed against the source reports for GSE1009, GSE30122, GSE96804, GSE142025, GSE163603 and GSE166239 [14–19]. Reported examples include Albert Einstein College of Medicine/Montefiore IRB 2002-202 for GSE30122, Shanghai Jiao Tong University Affiliated Sixth People's Hospital review and written consent for GSE142025, Indiana University IRB 1906572234 for GSE163603, and REK vest 2013/553 with written consent for GSE166239. The author must confirm whether a local institutional determination is required before submission. Absence of a located statement is a provenance gap, not an allegation of non-compliance.

### Pathway family and expression processing

Seven Reactome pathways were fixed in M19: Complement cascade (R-HSA-166658), Coagulation pathway (R-HSA-9769740), Cell surface interactions at the vascular wall (R-HSA-202733), Chemokine receptors bind chemokines (R-HSA-380108), Extracellular matrix organization (R-HSA-1474244), Cellular response to hypoxia (R-HSA-1234174), and Signaling by TGF-beta Receptor Complex (R-HSA-170834). The 21 June 2026 Reactome GMT was filtered to symbols in the 28 May 2026 GOA human GAF, yielding 783 unique genes [7]. Source hashes, retrieval dates and membership are archived. These rules were not prospectively preregistered (Supplementary Table S24).

Expression matrices were analyzed within cohorts. Microarray probes were mapped to uppercase gene symbols; the primary label-independent aggregation retained the probe with the highest all-sample mean and median-probe aggregation was a sensitivity. Author count tables were transformed as log2(CPM+0.5). Apparent GSE1009 a/b technical replicates were averaged by donor label. Hedges' g represented DKD minus control standardized mean difference [9]. Two-sided Welch P values were descriptive.

### Gene synthesis

For each of the 783 family members, cohort effects were combined by restricted-maximum-likelihood random-effects meta-analysis. Modified Hartung–Knapp intervals did not become narrower than conventional normal intervals [10,11]. Genes unavailable in one or more primary sources remained in the planned family with P=1 for Benjamini–Hochberg correction across all 783 members [12]. Leave-one-source-out estimates used the same model; two-source results were interpreted only as sensitivity diagnostics.

### Common-measurement pathway inference

For each pathway, the primary estimand was the arithmetic mean of gene-level Hedges' g over symbols measurable in all three primary cohorts. This prevents platform-specific membership from changing the quantity compared across sources. Cohort-specific mapped-set means were sensitivity estimands.

Disease labels were permuted jointly across samples so the observed cross-gene correlation was preserved [13]. GSE96804 allocations preserved case counts separately within archived sex strata. GSE30528 used exact case-count-preserving enumeration. GSE104948 was restricted to H7 before exact enumeration. For every pathway, the observed mean was centered and scaled by the permutation-null mean and standard deviation. The reference statistic was the maximum absolute studentized statistic across all seven pathways; its tail probability provided two-sided maxT family-wise control. Exact probabilities included the observed allocation and therefore could not be zero. Monte Carlo probabilities used the plus-one convention. GSE96804 used 100,000 allocations; its Monte Carlo standard errors and exact binomial intervals are in Supplementary Table S9. Raw-mean maxT is retained only as a sensitivity.

Within each cohort, 5,000 case/control-stratified bootstrap resamples estimated percentile 95% intervals for pathway means. These intervals describe sampling uncertainty conditional on the archived cohort and are not multiplicity-adjusted. A pathway met the operational rule only when at least two of three independent sources had studentized maxT P<0.05 in the same direction. “Criterion met” is not used as a synonym for causal validation.

### Software and reproducibility

The M20 primary robustness analysis used Python 3.12.12, NumPy 2.4.6, pandas 3.0.3 and SciPy 1.18.0; figures used Matplotlib 3.10.9. A frozen minimal input bundle contains only the three primary matrices, sample designs and authoritative canonical-pathway GMT with SHA-256 checksums. The code archive contains the analysis, regression tests, inputs and a clean-run README. It therefore does not depend on unbundled legacy helpers or future downloads of mutable reference files.

OpenAI Codex was used under author supervision for code generation, comparison of prespecified statistical options, consistency checking, figures, tables and language drafting. It also executed the explicitly labelled second-pass screening consistency audit. No generative image model was used. AI did not autonomously determine which inferential conclusions to retain or replace author scholarly judgement. The author retains responsibility for final verification of code, results, references, screening, interpretation and prose; AI is not an author and did not generate or alter primary data.

## Data availability

All expression data are publicly available in NCBI GEO: GSE1009, GSE30528, GSE30529, GSE96804, GSE104948, GSE104954, GSE111154, GSE142025, GSE163603, GSE166239 and GSE199838 [20–30]. Each repository record and, where available, its source publication are formally cited below and tabulated in Supplementary Table S25. No new primary data were generated.

## Code availability

The exact M20.1 code, frozen primary inputs, source tables, search logs, manuscript and checksums are supplied with the submission and archived in the versioned GitHub release https://github.com/denglizhen-113/dkd-focused-universe-reanalysis/tree/v1.2.1. No DOI has been assigned; a DOI-bearing archive is recommended but is not claimed.

## References

1. Alicic, R. Z., Rooney, M. T. & Tuttle, K. R. Diabetic kidney disease: challenges, progress, and possibilities. *Clin. J. Am. Soc. Nephrol.* **12**, 2032–2045 (2017). https://doi.org/10.2215/CJN.11491116
2. Li, Y. et al. Integrative transcriptome analysis reveals TEKT2 and PIAS2 involvement in diabetic nephropathy. *FASEB J.* **36**, e22592 (2022). https://doi.org/10.1096/fj.202200740RR
3. Hojjati, F., Roointan, A., Gholaminejad, A., Eshraghi, Y. & Gheisari, Y. Identification of key genes and biological regulatory mechanisms in diabetic nephropathy: meta-analysis of gene expression datasets. *Nefrología* **43**, 575–586 (2023). https://doi.org/10.1016/j.nefro.2022.06.003
4. Abdalla, M. et al. A common glomerular transcriptomic signature distinguishes diabetic kidney disease from other kidney diseases in humans and mice. *Curr. Res. Transl. Med.* **68**, 225–236 (2020). https://doi.org/10.1016/j.retram.2020.05.001
5. Zhang, L., Sun, Z., Yuan, Y. & Sheng, J. Integrating bioinformatics and machine learning to identify glomerular injury genes and predict drug targets in diabetic nephropathy. *Sci. Rep.* **15**, 16868 (2025). https://doi.org/10.1038/s41598-025-01628-5
6. Jiao, Y. et al. Activation of complement C1q and C3 in glomeruli might accelerate the progression of diabetic nephropathy: evidence from transcriptomic data and renal histopathology. *J. Diabetes Investig.* **13**, 839–849 (2022). https://doi.org/10.1111/jdi.13739
7. Gillespie, M. et al. The Reactome pathway knowledgebase 2022. *Nucleic Acids Res.* **50**, D687–D692 (2022). https://doi.org/10.1093/nar/gkab1028
8. Page, M. J. et al. The PRISMA 2020 statement: an updated guideline for reporting systematic reviews. *BMJ* **372**, n71 (2021). https://doi.org/10.1136/bmj.n71
9. Hedges, L. V. Distribution theory for Glass's estimator of effect size and related estimators. *J. Educ. Stat.* **6**, 107–128 (1981). https://doi.org/10.2307/1164588
10. IntHout, J., Ioannidis, J. P. A. & Borm, G. F. The Hartung–Knapp–Sidik–Jonkman method for random effects meta-analysis. *BMC Med. Res. Methodol.* **14**, 25 (2014). https://doi.org/10.1186/1471-2288-14-25
11. Röver, C., Knapp, G. & Friede, T. Hartung–Knapp–Sidik–Jonkman approach and its modification for random-effects meta-analysis with few studies. *BMC Med. Res. Methodol.* **15**, 99 (2015). https://doi.org/10.1186/s12874-015-0091-1
12. Benjamini, Y. & Hochberg, Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J. R. Stat. Soc. B* **57**, 289–300 (1995). https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
13. Winkler, A. M., Ridgway, G. R., Webster, M. A., Smith, S. M. & Nichols, T. E. Permutation inference for the general linear model. *NeuroImage* **92**, 381–397 (2014). https://doi.org/10.1016/j.neuroimage.2014.01.060
14. Baelde, H. J. et al. Gene expression profiling in glomeruli from human kidneys with diabetic nephropathy. *Am. J. Kidney Dis.* **43**, 636–650 (2004). https://doi.org/10.1053/j.ajkd.2003.12.028
15. Woroniecka, K. I. et al. Transcriptome analysis of human diabetic kidney disease. *Diabetes* **60**, 2354–2369 (2011). https://doi.org/10.2337/db10-1181
16. Pan, Y. et al. Dissection of glomerular transcriptional profile in patients with diabetic nephropathy: SRGAP2a protects podocyte structure and function. *Diabetes* **67**, 717–730 (2018). https://doi.org/10.2337/db17-0755
17. Fan, Y. et al. Comparison of kidney transcriptomic profiles of early and advanced diabetic nephropathy reveals potential new mechanisms for disease progression. *Diabetes* **68**, 2301–2314 (2019). https://doi.org/10.2337/db19-0204
18. Barwinska, D. et al. Molecular characterization of the human kidney interstitium in health and disease. *Sci. Adv.* **7**, eabd3359 (2021). https://doi.org/10.1126/sciadv.abd3359
19. Nordbø, O. P. et al. Transcriptomic analysis reveals partial epithelial–mesenchymal transition and inflammation as common pathogenic mechanisms in hypertensive nephrosclerosis and type 2 diabetic nephropathy. *Physiol. Rep.* **11**, e15825 (2023). https://doi.org/10.14814/phy2.15825
20. National Center for Biotechnology Information. Gene Expression Omnibus, GSE1009: Gene expression profiling in glomeruli from human kidneys with diabetic nephropathy (2004). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1009
21. National Center for Biotechnology Information. Gene Expression Omnibus, GSE30528: Diabetic nephropathy glomeruli (2011). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30528
22. National Center for Biotechnology Information. Gene Expression Omnibus, GSE30529: Diabetic nephropathy tubuli (2011). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30529
23. National Center for Biotechnology Information. Gene Expression Omnibus, GSE96804: Glomerular transcriptional profiles in diabetic nephropathy (2018). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96804
24. National Center for Biotechnology Information. Gene Expression Omnibus, GSE104948: Human glomerular transcriptomes from kidney disease cohorts (2017). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104948
25. National Center for Biotechnology Information. Gene Expression Omnibus, GSE104954: Human tubulointerstitial transcriptomes from kidney disease cohorts (2017). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104954
26. National Center for Biotechnology Information. Gene Expression Omnibus, GSE111154: RNA expression data for early diabetic nephropathy (2018). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111154
27. National Center for Biotechnology Information. Gene Expression Omnibus, GSE142025: Kidney transcriptomic profiles of early and advanced diabetic nephropathy (2019). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142025
28. National Center for Biotechnology Information. Gene Expression Omnibus, GSE163603: Molecular characterization of kidney interstitium in health and disease (2020). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163603
29. National Center for Biotechnology Information. Gene Expression Omnibus, GSE166239: Renal transcriptomes in nephrosclerosis and diabetic nephropathy (2021). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166239
30. National Center for Biotechnology Information. Gene Expression Omnibus, GSE199838: RNA sequencing of kidney tissue in diabetic kidney disease and controls (2023). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199838

## Acknowledgements

The author acknowledges the investigators and participants whose data were made publicly available.

## Funding

The author received no specific funding for this work.

## Author contributions

Lizhen Deng: conceptualization, methodology, software, formal analysis, investigation, data curation, visualization, writing—original draft, writing—review and editing, and project administration.

## Competing interests

The author declares no competing interests.

## Figure legends

**Figure 1. PRISMA 2020 flow adapted to dataset records.** The reporting unit is a repository dataset record, not a publication. PubMed records were used only to discover accessions. Counts show identification, removal before screening, screening, full accession/sample-metadata assessment and inclusion. All record-level decisions are in Supplementary Table S1.

**Figure 2. Compartment-specific study architecture.** **a,** Mirrored bars show control and DKD group sizes; color denotes renal compartment. **b,** Each row is an independent source study and each symbol an analyzed GEO Series. Shape denotes primary, contextual or sensitivity role; area scales with sample size. Paired compartments from GSE30122 and ERCB H7 were not counted as independent source effects.

**Figure 3. Descriptive gene-level display; none survives multiplicity correction.** **a,** Pooled Hedges' g and 95% modified Hartung–Knapp intervals for the 12 complete genes with the lowest unadjusted P values; selection is descriptive. **b,** All 582 complete three-source estimates. Benjamini–Hochberg correction retained the full 783-gene family, assigning P=1 to nonestimable members; 0/783 met FDR<0.05.

**Figure 4. Common-measurement pathway associations and operational study-wise criterion.** **a,** Mean Hedges' g over genes measurable in all three primary cohorts with 5,000-resample bootstrap 95% intervals. Color and marker identify source; filled markers denote two-sided studentized maxT family-wise P<0.05 and open markers P≥0.05. Intervals are descriptive and not multiplicity-adjusted. **b,** Number of sources meeting maxT P<0.05. The dashed line marks the operational two-source threshold. Meeting the rule denotes a recurring transcript association under stated assumptions, not causal or biochemical validation.
