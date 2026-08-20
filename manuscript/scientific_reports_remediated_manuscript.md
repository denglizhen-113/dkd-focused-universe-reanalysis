# Focused-Universe Reanalysis of Complement, Coagulation, and Extracellular-Matrix Transcriptional Programs in Diabetic Kidney Disease

Lizhen Deng

College of Life Science and Technology, Huazhong University of Science and Technology, Wuhan, Hubei, China

Corresponding author: Lizhen Deng (3070116993@qq.com)

ORCID: 0009-0003-2428-8176

## Abstract

Outcome-dependent gene filtering can invalidate conventional confirmation claims when discovery datasets are reused in meta-analysis. We reanalyzed nine biologically defined gene sets comprising 145 unique genes in human diabetic kidney disease transcriptomes. Discovery direction was estimated from GSE142025 and GSE96804. A selection-independent three-study synthesis added the prespecified glomerular GSE30528 compartment of source study GSE30122 and used restricted maximum likelihood with modified Hartung–Knapp inference, prediction intervals, and Benjamini–Hochberg correction across all 145 genes. Complete three-study estimates were available for 125/145 genes; none met meta-analysis FDR<0.05 (0/145). In the single external source study, 141/145 genes mapped, 94/141 followed the discovery direction, and 39/145 met external Welch-test FDR<0.05. Sample-label permutation, which preserved gene correlation, supported 8/9 predefined pathways after pathway-level FDR correction; leukocyte adhesion did not meet the threshold. These results support pathway-level transcriptional concordance within one external glomerular study but do not establish independent multi-study confirmation, causality, pathway activity, or cell-type localization. The analysis illustrates how broader, selection-independent families and small-sample meta-analytic inference materially alter conclusions drawn from a filtered gene subset.

## Keywords

Diabetic kidney disease; complement; coagulation; extracellular matrix; transcriptomics; meta-analysis; permutation testing

## Introduction

Diabetic kidney disease (DKD) develops through interacting glomerular, tubular, vascular, inflammatory, and fibrotic processes [1]. Human kidney transcriptomic studies have reported complement-associated, immune, vascular, and extracellular-matrix expression changes across DKD stages and tissue compartments [2–5]. Cross-study synthesis can help distinguish recurrent signals from dataset-specific findings, but its inferential validity depends on keeping discovery, filtering, and validation roles explicit.

A central difficulty arises when genes are selected using results from particular datasets and those same datasets are then reused to calculate ordinary meta-analytic P values for the selected subset. The resulting P values are conditional on an outcome-dependent filter and cannot be interpreted as if the gene family had been fixed independently of the observed results. A second difficulty occurs when two tissue compartments from the same parent study, including measurements from shared donors, are treated as independent cohorts.

We addressed both issues by returning to the complete 145-gene union encoded by nine biologically defined focused gene sets, before the outcome-dependent filter that retained 61 genes. Because the project lacks an external registration timestamp, this union is described as analysis-defined rather than preregistered. We prespecified one GSE30122 compartment—GSE30528 glomeruli—as the external estimand, used small-sample random-effects inference, and evaluated pathway-level concordance with sample-label permutation. Single-nucleus analyses were excluded from the primary claim because the available workflow lacked reference mapping, doublet detection, expanded quality control, and a complete donor-level pseudobulk analysis.

## Results

### Analysis family and study roles

The focused universe contained 145 unique genes distributed across nine overlapping biological pathways. The earlier outcome-dependent workflow retained 61 genes using GSE142025 significance and GSE96804 directional consistency. In the present analysis, all 145 genes remained in the inferential family. GSE142025 and GSE96804 defined discovery direction. GSE30528, the glomerular compartment of source study GSE30122, was fixed as the external compartment before the new result tables were inspected. GSE30529 was treated as a secondary tubulointerstitial dataset from the same source study and was not entered as an additional independent study (Fig. 1).

### Selection-independent three-study synthesis

Complete Hedges' g estimates from GSE142025, GSE96804, and GSE30528 were available for 125/145 genes. Random-effects models used restricted maximum likelihood and modified Hartung–Knapp t inference with two degrees of freedom. Non-estimable members were retained as P=1 for correction over the complete 145-gene family. No gene met meta-analysis FDR<0.05 (0/145; Fig. 2 and Supplementary Table S1). The unmodified HKSJ intervals and both unmodified and modified prediction intervals are supplied as sensitivities. Prediction intervals were wide because only three source studies contributed, underscoring that a stable mean effect across future settings could not be established for individual genes.

### Gene-level evaluation in one external source study

In GSE30528 glomeruli, 141/145 focused genes mapped to the platform and had estimable Welch tests. Of these, 94/141 followed the discovery direction, and 39/145 met FDR<0.05 after correction over the complete focused family (Fig. 3 and Supplementary Table S2). These are external results relative to the GSE142025/GSE96804 discovery direction, but they derive from one source study and therefore do not establish replication across independent source studies.

### Donor-disjoint tubulointerstitial sensitivity

Nine cross-accession donor pairs were inferred from sample titles and the shared GSE30122 parent study. Removing the corresponding GSE30529 samples left 5 DKD and 8 control tubulointerstitial samples. In this donor-disjoint sensitivity, 141/145 genes mapped, 86/141 followed the discovery direction, and 20/145 met within-compartment FDR<0.05. Sample-label permutation supported 8/9 pathways after pathway-level FDR correction (Supplementary Tables S9–S11). These results demonstrate non-uniform gene-level behavior across compartments: 7 genes met FDR<0.05 in both compartments, 32 only in GSE30528, and 13 only in donor-disjoint GSE30529. They do not convert the two compartments into independent studies.

### Pathway-level sample-label permutation

To retain gene-gene correlation, disease labels were permuted among the 22 GSE30528 samples while preserving the 9-DKD versus 13-control allocation, and all genes were recomputed jointly in each of 10,000 permutations. Eight of nine pathways met pathway-level FDR<0.05 (Fig. 4; Table 1). Complement activation, endothelial activation, macrophage recruitment, fibrosis, extracellular-matrix remodeling, coagulation, hypoxia, and transforming-growth-factor-beta signaling were supported. Leukocyte adhesion was directionally positive but did not meet the adjusted threshold (FDR=0.0998).

### Unsupervised structure, outlier, and WGCNA audits

PCA and sample-correlation diagnostics were run separately within each primary dataset. The predefined median-correlation rule flagged 8 samples; because all flagged samples were DKD samples and their biological versus technical origin could not be distinguished, they were retained in the primary analysis. Excluding them in a sensitivity analysis did not change the individual-gene conclusion (0/145 meta-analysis FDR<0.05; Supplementary Figs. S2–S3 and Supplementary Tables S14–S16).

WGCNA was retained only as an internal GSE96804 description. Of the 145 focused genes, 58 mapped to the DKD-associated black module. Nineteen of these also met GSE30528 FDR<0.05, whereas 39 black-module genes did not; conversely, 20 GSE30528 FDR-positive genes were outside the black module or outside the WGCNA top-variance universe (Supplementary Table S12). This partial overlap is compatible with compartment, cohort, variance-filter, and network-construction differences and is not independent validation.

### Table 1. External pathway-level permutation results

| Pathway | Mapped genes | Aligned Hedges' g | Permutation P | BH FDR |
| --- | --- | --- | --- | --- |
| COMPLEMENT_ACTIVATION | 21 | 1.02 | 9.999e-05 | 0.00018 |
| ENDOTHELIAL_ACTIVATION | 16 | 0.83 | 9.999e-05 | 0.00018 |
| MACROPHAGE_RECRUITMENT | 16 | 0.8147 | 9.999e-05 | 0.00018 |
| FIBROSIS | 12 | 1.073 | 9.999e-05 | 0.00018 |
| ECM_REMODELING | 15 | 0.6398 | 9.999e-05 | 0.00018 |
| COAGULATION_CASCADE | 21 | 0.4985 | 0.0003 | 0.00045 |
| HYPOXIA | 10 | 0.4689 | 0.0022 | 0.002828 |
| TGF_BETA_SIGNALING | 14 | 0.3593 | 0.0314 | 0.03532 |
| LEUKOCYTE_ADHESION | 11 | 0.2161 | 0.09979 | 0.09979 |

## Discussion

The principal finding is narrower than the earlier filtered-subset interpretation. When the full 145-gene family was analyzed with REML and modified Hartung–Knapp inference, no individual gene met meta-analysis FDR<0.05. The result does not imply that all underlying expression differences are absent. It shows that the available three-study evidence is insufficient for multiplicity-adjusted, small-sample random-effects claims at individual-gene resolution.

The independent evidence is stronger at pathway level within GSE30528. Eight predefined pathways were concordant under a sample-label permutation that retained the observed cross-gene structure. This supports an exploratory interpretation in which complement, coagulation, vascular, inflammatory, hypoxia, transforming-growth-factor-beta, and extracellular-matrix programs recur in the glomerular compartment of one external DKD study. Prior human kidney studies provide biological context for these themes [2–5], but the present reanalysis does not establish pathway activation, protein abundance, upstream regulation, or causality.

The contrast between 0/145 meta-analytic gene discoveries and 39/145 external Welch-test discoveries is not contradictory. The former uses random-effects inference over only three heterogeneous source studies with a t reference distribution and full-family multiplicity control. The latter tests one external compartment and does not estimate between-study heterogeneity. The pathway permutation then aggregates direction-aligned effects and gains stability while preserving gene correlation. These estimands answer different questions and should not be combined into a single “confirmed gene” count.

Model multiplicity was addressed by fixing one primary model—REML with modified Hartung–Knapp inference—and one 145-gene correction family. Unmodified HKSJ intervals, prediction intervals, sex adjustment, and correlation-outlier exclusion are labeled sensitivities and were not searched to choose a favorable result. A Bonferroni division by the six historical meta-analysis variants is therefore not applied to the new primary model; the historical variants no longer constitute six co-primary opportunities for discovery.

GSE30528 and GSE30529 should not be presented as separate sources of validation. Both are compartments of GSE30122, and nine cross-accession donor pairs can be inferred from sample titles and the shared parent study. Restricting the primary external estimand to GSE30528 avoids introducing an unestimated cross-compartment covariance. The earlier filtered 52-gene analysis reported 39/52 direction concordance in donor-disjoint GSE30529 but did not meet its matched-background criterion. That statistic is retained only as historical context because the candidate set and null construction were outcome-dependent; it is not combined with the new 145-gene tests.

The new full-family sensitivity clarifies, but does not erase, the compartment difference. Donor removal reduced GSE30529 from 10 DKD and 12 control samples to 5 and 8, respectively, making single-gene estimates less precise and more sensitive to individual samples. The two microdissected compartments also contain different cell mixtures and may carry genuinely different DKD responses: glomerular endothelial, mesangial, and podocyte programs need not track tubulointerstitial epithelial, stromal, vascular, and immune signals gene by gene. Manual microdissection, residual vascular content, probe behavior, and incomplete clinical metadata provide additional technical or compositional explanations. Consistent with this heterogeneity, only 7 genes met full-family FDR<0.05 in both compartments, whereas 32 were GSE30528-only and 13 were donor-disjoint-GSE30529-only. At pathway level, however, 8/9 sets were supported in each compartment, although the unsupported set differed. This pattern suggests that broad biological programs may recur while their contributing genes and effect magnitudes vary by compartment. It does not prove a compartment-specific mechanism, because donor pairing was inferred from titles, covariates were unavailable, and the retained GSE30529 sample was small. The correct conclusion is therefore same-study, compartment-sensitive concordance—not universal DKD replication and not evidence from two independent validation cohorts.

Across the nine paired donors, the median gene-wise cross-compartment expression correlation was only 0.219. This descriptive quantity is not the correlation between cohort-level disease-effect estimates and therefore cannot identify the rho parameter used in the retired correlated-compartment meta-analysis. Rather than choosing rho from these nine inferred pairs, the primary analysis removed rho dependence by representing GSE30122 with one prespecified compartment. The former rho-grid models are retained only as historical sensitivity analyses and are not used to select genes or conclusions.

Single-nucleus results were removed from the title, abstract, primary results, and principal conclusion. Existing marker-score labels may be shown separately as provisional methodological context only. Cell-type localization should be reconsidered after external reference mapping, doublet detection, expanded nucleus-level quality control, donor-aware integration, and donor-level pseudobulk analysis [12,13].

## Methods

### Focused universe and discovery boundary

The analysis-defined universe was the union of nine pathway gene sets covering complement activation, coagulation, endothelial activation, macrophage recruitment, leukocyte adhesion, fibrosis, extracellular-matrix remodeling, hypoxia, and transforming-growth-factor-beta signaling. The union contained 145 unique genes. Repository lineage placed this union before the result-dependent 145-to-61 filter, but no external preregistration or trusted third-party timestamp was available. All 145 genes were retained for multiplicity adjustment.

GSE142025 and GSE96804 were the discovery datasets. For each gene, discovery direction was the sign of the inverse-variance fixed-effect mean of the two Hedges' g estimates. GSE30528 glomeruli were the prespecified external compartment from parent study GSE30122. GSE30529 tubulointerstitium was excluded from the primary estimator because it belongs to the same source study and includes repeated donors.

### Expression contrasts and effect sizes

Human kidney datasets were analyzed separately rather than merged across platforms. GSE142025 compared 21 advanced-DKD samples with 9 controls. GSE96804 compared 41 DKD glomerular samples with 20 controls. GSE30528 compared 9 DKD with 13 control glomerular samples. For GEO microarrays, platform gene symbols were uppercased and, when multiple probes mapped to one gene, the probe with the highest mean expression across all samples was selected without using disease labels. Each study-level effect was Hedges' g for the case-minus-control standardized mean difference [8]. External gene-level P values used two-sided Welch tests.

### Random-effects synthesis and multiplicity

For each gene with all three effects estimable, between-study variance was estimated by restricted maximum likelihood. The pooled effect used inverse total-variance weights. Because only three studies contributed and their precisions differed, primary uncertainty used modified Hartung–Knapp with scale max(1, q), a t reference distribution, and k−1 degrees of freedom; unmodified HKSJ results were retained as a sensitivity [9,10]. Prediction intervals used k−2 degrees of freedom. Genes without complete three-study mapping were assigned P=1 for family-wise bookkeeping. Benjamini–Hochberg correction was applied to exactly 145 P values [7].

### External gene and pathway evaluation

External gene-level Welch P values in GSE30528 were corrected across all 145 genes, with unmapped genes retained as P=1. Direction concordance compared the GSE30528 Hedges' g sign with the discovery direction.

For pathway testing, the statistic was the mean discovery-direction-aligned Hedges' g among mapped genes in each pathway. Disease labels were permuted 10,000 times among GSE30528 samples while retaining nine case labels. All genes were recomputed jointly for every permutation, retaining gene correlation induced by the observed expression matrix. One-sided P values used the plus-one correction, and the nine pathway P values were BH-adjusted.

### Donor-disjoint and diagnostic sensitivities

Donor numbers were parsed from GSE30528 and GSE30529 sample titles. The nine numbers appearing in both accessions were treated as inferred repeated donors. Those nine GSE30529 samples were removed before the donor-disjoint tubulointerstitial sensitivity, leaving 5 DKD and 8 control samples. Welch tests and BH correction were repeated across all 145 genes, and the nine pathway statistics were recomputed in 10,000 joint sample-label permutations with a fixed 5-versus-8 allocation. Gene-wise Pearson correlations across the nine paired compartment measurements were calculated descriptively; they were not treated as estimates of cross-study effect covariance.

Within each primary dataset, PCA and Pearson sample-correlation matrices used the 1,000 genes with highest within-dataset variance. Potential correlation outliers were defined before inspecting downstream sensitivity results as samples with a robust z score below −3 for median correlation to all other samples. Flagged samples were not automatically removed because disease biology and technical artifacts could not be distinguished. The primary meta-analysis was repeated after their exclusion as a sensitivity. GSE96804 additionally permitted an intercept-plus-disease-plus-sex linear model; no common adjusted model was estimable across datasets.

WGCNA module membership was imported from the existing GSE96804 network analysis [11]. Overlap with the 145-gene universe and GSE30528 gene-level results was summarized descriptively. Because the network and one discovery direction component both came from GSE96804, no WGCNA P value was interpreted as external confirmation.

### Reproducibility

The analysis archived the focused universe, per-study effects, discovery directions, complete random-effects output, external gene tests, pathway permutation results, donor-disjoint sensitivities, sample-level metadata, diagnostics, analysis settings, input hashes, and an output manifest. The execution environment used Python 3.13.9, NumPy 2.3.5, pandas 2.3.3, SciPy 1.16.3, statsmodels 0.14.5, and Matplotlib 3.10.6. M18 recomputed statistical results from frozen processed matrices and cached normalized GEO SOFT values; it did not independently repeat every upstream raw-read or CEL-file preprocessing step. The public versioned code-and-results repository is available at https://github.com/denglizhen-113/dkd-focused-universe-reanalysis (version tag v1.0.2).

### Software and archival limitations

Historical discovery-stage limma calls did not explicitly archive the `trend` and `robust` arguments and therefore depended on the recorded package version's defaults. Some older Python stages recorded major/minor but not patch versions. The prior single-nucleus phase did not archive a final h5ad object and was reconstructed from intermediate files; those results are excluded from the present claims. The present reanalysis records its full Python patch versions and does not use the historical limma P values or reconstructed single-nucleus labels for primary inference. These changes improve current reproducibility but do not retrospectively make every upstream historical step raw-data reproducible.

## Limitations

The 145-gene universe has internal pipeline lineage but no external registration timestamp; it should not be described as preregistered. Only three independent source studies contributed to the primary synthesis, making heterogeneity estimates and prediction intervals imprecise even with small-sample methods. GSE30528 provides only one external source study, so its two renal compartments cannot establish multiple-study replication. The paired-donor expression correlations cannot identify correlation between compartment-level disease-effect estimates; the primary model therefore avoids a rho parameter rather than claiming to estimate it.

Individual-level batch, age, treatment, renal-function, and other covariates were incomplete or absent in the public annotations. GSE96804 recorded sex, but the other primary datasets did not supply a common covariate set. Residual confounding and technical structure therefore cannot be excluded. Correlation diagnostics flagged eight DKD samples, but these were not assumed to be technical failures. Probe selection by highest all-sample mean was label-independent but may favor highly expressed isoforms and does not resolve transcript-level effects. The primary computations start from processed expression matrices or normalized GEO values, not a uniform raw-data reprocessing pipeline.

The nine pathways overlap in gene membership. The sample-label permutation preserved within-dataset gene correlation, but pathway-level BH correction did not provide a max-statistic family-wise test across overlapping sets. The results are transcript-level associations and do not establish protein abundance, pathway activity, molecular regulation, or causality.

The single-nucleus workflow available for this project lacked reference mapping, doublet detection, expanded UMI and mitochondrial quality thresholds, donor-aware integration, and a complete donor-level pseudobulk analysis. It therefore cannot support cell-type localization claims in the present manuscript.

## Data availability

All analyzed datasets are publicly available from the Gene Expression Omnibus:

- GSE142025: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142025
- GSE96804: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96804
- GSE30528: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30528
- GSE30529: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30529
- GSE166239: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166239
- GSE131882: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131882
- GSE111154: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111154

Processed effect-size tables, pathway-permutation results, and analysis manifests are available in the public versioned repository at https://github.com/denglizhen-113/dkd-focused-universe-reanalysis/tree/v1.0.2.

## Code availability

Analysis code, environment specifications, source-data tables, and machine-readable manifests are publicly available at https://github.com/denglizhen-113/dkd-focused-universe-reanalysis/tree/v1.0.2 (version tag v1.0.2). A DOI-minting archival record may be added in a subsequent release.

## References

1. Alicic, R. Z., Rooney, M. T. & Tuttle, K. R. Diabetic kidney disease: challenges, progress, and possibilities. *Clin. J. Am. Soc. Nephrol.* **12**, 2032–2045 (2017). https://doi.org/10.2215/CJN.11491116
2. Woroniecka, K. I. et al. Transcriptome analysis of human diabetic kidney disease. *Diabetes* **60**, 2354–2369 (2011). https://doi.org/10.2337/db10-1181
3. Pan, Y. et al. Dissection of glomerular transcriptional profile in patients with diabetic nephropathy: SRGAP2a protects podocyte structure and function. *Diabetes* **67**, 717–730 (2018). https://doi.org/10.2337/db17-0755
4. Fan, Y. et al. Comparison of kidney transcriptomic profiles of early and advanced diabetic nephropathy reveals potential new mechanisms for disease progression. *Diabetes* **68**, 2301–2314 (2019). https://doi.org/10.2337/db19-0204
5. Nordbø, O. P. et al. Transcriptomic analysis reveals partial epithelial-mesenchymal transition and inflammation as common pathogenic mechanisms in hypertensive nephrosclerosis and type 2 diabetic nephropathy. *Physiol. Rep.* **11**, e15825 (2023). https://doi.org/10.14814/phy2.15825
6. Ritchie, M. E. et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res.* **43**, e47 (2015). https://doi.org/10.1093/nar/gkv007
7. Benjamini, Y. & Hochberg, Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J. R. Stat. Soc. B* **57**, 289–300 (1995). https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
8. Hedges, L. V. Distribution theory for Glass's estimator of effect size and related estimators. *J. Educ. Stat.* **6**, 107–128 (1981). https://doi.org/10.2307/1164588
9. IntHout, J., Ioannidis, J. P. A. & Borm, G. F. The Hartung-Knapp-Sidik-Jonkman method for random effects meta-analysis is straightforward and considerably outperforms the standard DerSimonian-Laird method. *BMC Med. Res. Methodol.* **14**, 25 (2014). https://doi.org/10.1186/1471-2288-14-25
10. Röver, C., Knapp, G. & Friede, T. Hartung-Knapp-Sidik-Jonkman approach and its modification for random-effects meta-analysis with few studies. *BMC Med. Res. Methodol.* **15**, 99 (2015). https://doi.org/10.1186/s12874-015-0091-1
11. Langfelder, P. & Horvath, S. WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics* **9**, 559 (2008). https://doi.org/10.1186/1471-2105-9-559
12. Squair, J. W. et al. Confronting false discoveries in single-cell differential expression. *Nat. Commun.* **12**, 5692 (2021). https://doi.org/10.1038/s41467-021-25960-2
13. Lähnemann, D. et al. Eleven grand challenges in single-cell data science. *Genome Biol.* **21**, 31 (2020). https://doi.org/10.1186/s13059-020-1926-6

## Acknowledgements

The author acknowledges the investigators who generated and publicly shared the datasets reanalyzed in this study.

## Funding

The author received no specific funding for this work.

## Ethics statement

This study used publicly available, de-identified human transcriptomic data. No new human or animal samples were collected.

## Author contributions

Lizhen Deng: Conceptualization, methodology, software, formal analysis, investigation, data curation, visualization, writing—original draft, writing—review and editing, and project administration.

## Competing interests

The author declares no competing interests.

## AI disclosure

OpenAI Codex was used under author supervision for code generation, comparison of statistical options, statistical remediation, consistency checks, figure and table preparation, and language drafting. Figures were generated by AI-assisted code from the archived numerical tables; no generative image model was used. The author retains scholarly judgement and responsibility and must independently verify the code, numerical outputs, citations, interpretation, and final submitted text. AI was not used to generate or alter primary data.

## Figure legends

### Figure 1. Analysis boundaries used for statistical remediation.

The 145-gene analysis-defined universe was carried through discovery-direction estimation, the three-study synthesis, external gene testing, and pathway permutation. The outcome-dependent 61-gene subset was not used as the primary inferential family. GSE30529 was retained only as a secondary same-study compartment.

### Figure 2. Selection-independent three-study synthesis across the focused universe.

Panel a shows pooled Hedges' g against the modified Hartung–Knapp P value for genes with complete estimates in GSE142025, GSE96804, and GSE30528. Labels identify the ten smallest unadjusted P values for display only. Panel b shows every study effect, normal 95% confidence interval, random-effects weight, and the modified Hartung–Knapp pooled interval for eight fixed legacy display genes. These genes are not a separate inferential family. No gene met BH FDR<0.05 across the 145-gene family.

### Figure 3. Gene-level evaluation in GSE30528 glomeruli.

Discovery effects combine GSE142025 and GSE96804 without GSE30528. External effects are Hedges' g values from GSE30528. Blue points are direction-concordant genes; dark blue additionally meet external FDR<0.05; red points are direction-discordant. GSE30528 is one external source study.

### Figure 4. Pathway-level sample-label permutation in GSE30528.

Bars show the mean discovery-direction-aligned Hedges' g for each pathway. Labels report BH FDR over nine pathways. Disease labels were permuted 10,000 times with the 9-versus-13 allocation preserved and all genes recomputed jointly.

## Supplementary legend

### Supplementary Figure S1. Complete three-study estimates for 125 mapped genes.

Points show REML pooled Hedges' g and modified Hartung–Knapp 95% confidence intervals. The figure is descriptive; multiplicity-adjusted results are supplied in Supplementary Table S1.

### Supplementary Figure S2. Unsupervised PCA diagnostics.

PCA used the 1,000 most variable genes separately within each primary dataset. Group separation can reflect disease biology, composition, confounding, or unrecorded technical structure.

### Supplementary Figure S3. Within-dataset sample-correlation diagnostics.

Heatmaps show Pearson correlations computed from the 1,000 most variable genes within each primary dataset. The predefined robust-z outlier audit is reported in Supplementary Table S15.
