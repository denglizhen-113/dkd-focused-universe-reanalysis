from __future__ import annotations

import hashlib
import json
import math
import shutil
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch


ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = ROOT / "tables" / "stage21_m18_remediation"
FIGURE_DIR = ROOT / "figures" / "stage21_m18_remediation"
MANUSCRIPT_DIR = ROOT / "manuscript_ready" / "stage21_m18_remediation"
PACKAGE_DIR = ROOT / "submission_package" / "stage21_m18_remediation"
MAIN_FIGURE_DIR = PACKAGE_DIR / "main_figures"
SUPPLEMENT_DIR = PACKAGE_DIR / "supplementary"
SCRIPT_PATH = Path(__file__).resolve()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def save_figure(fig: plt.Figure, stem: str, *, dpi: int = 400) -> list[Path]:
    outputs = []
    for suffix in (".png", ".pdf"):
        path = FIGURE_DIR / f"{stem}{suffix}"
        fig.savefig(path, dpi=dpi if suffix == ".png" else None, bbox_inches="tight")
        outputs.append(path)
    plt.close(fig)
    return outputs


def figure_1_design() -> list[Path]:
    fig, ax = plt.subplots(figsize=(11.5, 4.8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    boxes = [
        (0.02, 0.60, 0.20, 0.25, "Analysis-defined universe\n145 genes in 9 pathways", "#dceaf7"),
        (0.28, 0.60, 0.20, 0.25, "Discovery direction\nGSE142025 + GSE96804", "#e8e4f4"),
        (0.54, 0.60, 0.20, 0.25, "Selection-independent synthesis\n+ GSE30528 glomeruli", "#e3f2df"),
        (0.80, 0.60, 0.18, 0.25, "REML + modified HK\nBH family = 145", "#f9ead7"),
        (0.28, 0.13, 0.20, 0.22, "External gene test\nGSE30528 only", "#fff4cf"),
        (0.54, 0.13, 0.20, 0.22, "Pathway label permutation\n10,000 replicates", "#fff4cf"),
        (0.80, 0.13, 0.18, 0.22, "GSE30529\nsecondary same-study context", "#eeeeee"),
    ]
    for x, y, width, height, label, color in boxes:
        patch = FancyBboxPatch(
            (x, y), width, height, boxstyle="round,pad=0.015", facecolor=color,
            edgecolor="#34495e", linewidth=1.2
        )
        ax.add_patch(patch)
        ax.text(x + width / 2, y + height / 2, label, ha="center", va="center", fontsize=10)
    arrows = [
        ((0.22, 0.725), (0.28, 0.725)),
        ((0.48, 0.725), (0.54, 0.725)),
        ((0.74, 0.725), (0.80, 0.725)),
        ((0.38, 0.60), (0.38, 0.35)),
        ((0.64, 0.60), (0.64, 0.35)),
        ((0.89, 0.60), (0.89, 0.35)),
    ]
    for start, end in arrows:
        ax.annotate("", xy=end, xytext=start, arrowprops={"arrowstyle": "->", "lw": 1.5, "color": "#34495e"})
    ax.text(
        0.5, 0.95, "Predefined analysis boundaries used for statistical remediation",
        ha="center", va="center", fontsize=13, fontweight="bold"
    )
    ax.text(
        0.02, 0.03,
        "The outcome-dependent 61-gene subset is not used as the primary inferential family. "
        "GSE30528 and GSE30529 are compartments of one source study.",
        fontsize=9, color="#333333"
    )
    return save_figure(fig, "Figure_1_analysis_design")


def figure_2_primary(primary: pd.DataFrame, per_cohort: pd.DataFrame) -> list[Path]:
    frame = primary.loc[primary["complete_three_study_mapping"]].copy()
    frame["minus_log10_p"] = -np.log10(frame["p_value_modified_hk"].clip(lower=1e-300))
    fig, (ax, forest_ax) = plt.subplots(
        1, 2, figsize=(14.5, 9.2), gridspec_kw={"width_ratios": [1.0, 1.38]}
    )
    ax.scatter(
        frame["pooled_effect"], frame["minus_log10_p"], s=30,
        c="#6c7a89", alpha=0.78, edgecolor="white", linewidth=0.35
    )
    ax.axvline(0, color="#444444", linewidth=0.9, linestyle="--")
    top = frame.nsmallest(10, "p_value_modified_hk")
    for row in top.itertuples(index=False):
        ax.annotate(
            row.gene_symbol,
            (row.pooled_effect, -math.log10(max(row.p_value_modified_hk, 1e-300))),
            xytext=(4, 4), textcoords="offset points", fontsize=8
        )
    ax.set_xlabel("Pooled Hedges' g (REML)")
    ax.set_ylabel("−log10 modified Hartung–Knapp P")
    ax.set_title("Selection-independent three-study synthesis across the 145-gene family")
    ax.text(
        0.02, 0.98,
        f"Complete mapping: {len(frame)}/145\nBH FDR<0.05: {int(primary['fdr_lt_0_05'].sum())}/145",
        transform=ax.transAxes, va="top", ha="left",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.9, "edgecolor": "#999999"},
        fontsize=10
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(-0.12, 1.04, "a", transform=ax.transAxes, fontsize=14, fontweight="bold")

    display_genes = ["C3", "COL1A1", "COL1A2", "FN1", "PECAM1", "PLG", "VCAM1", "VWF"]
    colors = {"GSE142025_C1": "#1b9e77", "GSE96804": "#7570b3", "GSE30528": "#d95f02"}
    labels = {"GSE142025_C1": "GSE142025", "GSE96804": "GSE96804", "GSE30528": "GSE30528"}
    forest_rows = []
    y = 0
    y_positions = []
    y_labels = []
    for gene in display_genes:
        gene_meta = primary.loc[primary["gene_symbol"].eq(gene)].iloc[0]
        gene_data = per_cohort.loc[
            per_cohort["gene_symbol"].eq(gene) & per_cohort["cohort"].isin(colors)
        ].copy()
        tau2 = float(gene_meta["tau2_reml"])
        valid = gene_data["variance"].notna() & gene_data["variance"].gt(0)
        total_weight = float((1.0 / (gene_data.loc[valid, "variance"] + tau2)).sum())
        for cohort in ("GSE142025_C1", "GSE96804", "GSE30528"):
            row = gene_data.loc[gene_data["cohort"].eq(cohort)].iloc[0]
            effect = float(row["hedges_g"])
            se = math.sqrt(float(row["variance"]))
            low, high = effect - 1.96 * se, effect + 1.96 * se
            weight = (1.0 / (float(row["variance"]) + tau2)) / total_weight if total_weight > 0 else math.nan
            forest_ax.errorbar(
                effect, y, xerr=[[effect - low], [high - effect]], fmt="o",
                color=colors[cohort], ecolor=colors[cohort],
                markersize=3.0 + 8.0 * math.sqrt(max(weight, 0.0)),
                elinewidth=1.0, capsize=1.8,
            )
            y_positions.append(y)
            y_labels.append(f"{gene}  {labels[cohort]}")
            forest_rows.append(
                {
                    "gene_symbol": gene, "estimate_type": "study", "study": cohort,
                    "hedges_g": effect, "ci_95_low": low, "ci_95_high": high,
                    "random_effect_weight": weight, "tau2_reml": tau2,
                }
            )
            y += 1
        pooled = float(gene_meta["pooled_effect"])
        low = float(gene_meta["ci_95_low_modified_hk"])
        high = float(gene_meta["ci_95_high_modified_hk"])
        forest_ax.errorbar(
            pooled, y, xerr=[[pooled - low], [high - pooled]], fmt="D",
            color="#111111", ecolor="#111111", markersize=4.4, elinewidth=1.2, capsize=2,
        )
        y_positions.append(y)
        y_labels.append(f"{gene}  pooled")
        forest_rows.append(
            {
                "gene_symbol": gene, "estimate_type": "modified_HK_pooled", "study": "pooled",
                "hedges_g": pooled, "ci_95_low": low, "ci_95_high": high,
                "random_effect_weight": 1.0, "tau2_reml": tau2,
            }
        )
        y += 1.0
        forest_ax.axhline(y - 0.45, color="#eeeeee", linewidth=0.7)
    pd.DataFrame(forest_rows).to_csv(TABLE_DIR / "forest_plot_display_genes.csv", index=False)
    forest_ax.axvline(0, color="#555555", linewidth=0.8, linestyle="--")
    forest_ax.set_yticks(y_positions, y_labels, fontsize=7.1)
    forest_ax.invert_yaxis()
    forest_ax.set_xlabel("Hedges' g with 95% CI")
    forest_ax.set_title("Per-study and pooled estimates for eight fixed display genes")
    forest_ax.spines[["top", "right"]].set_visible(False)
    forest_ax.text(-0.11, 1.04, "b", transform=forest_ax.transAxes, fontsize=14, fontweight="bold")
    forest_ax.text(
        0.01, -0.075,
        "Study intervals: normal 95% CI; marker area encodes relative REML weight (exact values: Table S17).\n"
        "Pooled intervals: modified Hartung–Knapp 95% CI. Fixed display set; not a testing family.",
        transform=forest_ax.transAxes, fontsize=8, va="top",
    )
    fig.tight_layout()
    return save_figure(fig, "Figure_2_primary_reml_hk")


def figure_3_external(external: pd.DataFrame) -> list[Path]:
    frame = external.loc[external["external_welch_p_value"].notna()].copy()
    colors = np.where(
        frame["external_fdr_lt_0_05"] & frame["direction_concordant"],
        "#2b8cbe",
        np.where(frame["direction_concordant"], "#a6bddb", "#e34a33"),
    )
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    ax.scatter(
        frame["discovery_effect_fixed"], frame["hedges_g"], c=colors,
        s=35, alpha=0.82, edgecolor="white", linewidth=0.35
    )
    ax.axhline(0, color="#555555", linestyle="--", linewidth=0.8)
    ax.axvline(0, color="#555555", linestyle="--", linewidth=0.8)
    limit = float(np.nanmax(np.abs(np.r_[frame["discovery_effect_fixed"], frame["hedges_g"]])))
    ax.plot([-limit, limit], [-limit, limit], color="#888888", linewidth=0.7, linestyle=":")
    ax.set_xlabel("Discovery fixed-effect Hedges' g\n(GSE142025 + GSE96804)")
    ax.set_ylabel("External glomerular Hedges' g\n(GSE30528)")
    ax.set_title("Gene-level evaluation in one external source study")
    ax.text(
        0.02, 0.98,
        f"Mapped: {len(frame)}/145\nDirection concordant: {int(frame['direction_concordant'].sum())}/{len(frame)}\nExternal BH FDR<0.05: {int(frame['external_fdr_lt_0_05'].sum())}/145",
        transform=ax.transAxes, va="top", ha="left",
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.9, "edgecolor": "#999999"},
        fontsize=9.5
    )
    ax.spines[["top", "right"]].set_visible(False)
    return save_figure(fig, "Figure_3_external_gene_validation")


def figure_4_pathways(pathways: pd.DataFrame) -> list[Path]:
    frame = pathways.sort_values("observed_mean_direction_aligned_hedges_g").copy()
    colors = np.where(frame["fdr_bh_9_pathways"] < 0.05, "#3182bd", "#bdbdbd")
    fig, ax = plt.subplots(figsize=(8.4, 5.8))
    y = np.arange(len(frame))
    ax.barh(y, frame["observed_mean_direction_aligned_hedges_g"], color=colors, alpha=0.9)
    ax.set_yticks(y, [name.replace("_", " ").title() for name in frame["pathway"]])
    ax.axvline(0, color="#555555", linewidth=0.8)
    ax.set_xlabel("Mean direction-aligned Hedges' g in GSE30528")
    ax.set_title("Sample-label permutation validation of nine focused pathways")
    for index, row in enumerate(frame.itertuples(index=False)):
        ax.text(
            row.observed_mean_direction_aligned_hedges_g + 0.02,
            index,
            f"FDR={row.fdr_bh_9_pathways:.3g}",
            va="center", fontsize=8.5
        )
    ax.spines[["top", "right"]].set_visible(False)
    return save_figure(fig, "Figure_4_pathway_permutation")


def supplementary_figure_1(primary: pd.DataFrame) -> list[Path]:
    frame = primary.loc[primary["complete_three_study_mapping"]].sort_values("pooled_effect").copy()
    y = np.arange(len(frame))
    low = frame["pooled_effect"] - frame["ci_95_low_modified_hk"]
    high = frame["ci_95_high_modified_hk"] - frame["pooled_effect"]
    fig, ax = plt.subplots(figsize=(10, 25))
    ax.errorbar(
        frame["pooled_effect"], y, xerr=np.vstack([low, high]), fmt="o",
        markersize=2.8, color="#4d4d4d", ecolor="#bdbdbd", elinewidth=0.65, capsize=1.2
    )
    ax.axvline(0, color="#b2182b", linestyle="--", linewidth=0.8)
    ax.set_yticks(y, frame["gene_symbol"], fontsize=6.7)
    ax.set_xlabel("Pooled Hedges' g with modified Hartung–Knapp 95% CI")
    ax.set_title("Complete three-study estimates for 125 mapped genes")
    ax.spines[["top", "right"]].set_visible(False)
    return save_figure(fig, "Supplementary_Figure_S1_all_gene_intervals", dpi=300)


def markdown_table(frame: pd.DataFrame, columns: list[str], labels: list[str]) -> str:
    lines = ["| " + " | ".join(labels) + " |", "| " + " | ".join(["---"] * len(labels)) + " |"]
    for row in frame[columns].itertuples(index=False, name=None):
        rendered = []
        for value in row:
            if isinstance(value, (float, np.floating)):
                rendered.append(f"{value:.4g}")
            else:
                rendered.append(str(value))
        lines.append("| " + " | ".join(rendered) + " |")
    return "\n".join(lines)


def build_manuscript(primary: pd.DataFrame, external: pd.DataFrame, pathways: pd.DataFrame) -> str:
    tub = pd.read_csv(TABLE_DIR / "gse30529_donor_disjoint_gene_validation_145.csv")
    tub_pathways = pd.read_csv(TABLE_DIR / "gse30529_donor_disjoint_pathway_permutation.csv")
    paired = pd.read_csv(TABLE_DIR / "paired_compartment_expression_correlations_145.csv")
    wgcna = pd.read_csv(TABLE_DIR / "wgcna_descriptive_integration_145.csv")
    outlier_sensitivity = pd.read_csv(TABLE_DIR / "primary_outlier_exclusion_sensitivity_145.csv")
    sample_outliers = pd.read_csv(TABLE_DIR / "sample_correlation_outlier_diagnostics.csv")
    pathway_table = markdown_table(
        pathways.sort_values("fdr_bh_9_pathways"),
        [
            "pathway", "genes_mapped_with_discovery_direction",
            "observed_mean_direction_aligned_hedges_g", "one_sided_permutation_p",
            "fdr_bh_9_pathways",
        ],
        ["Pathway", "Mapped genes", "Aligned Hedges' g", "Permutation P", "BH FDR"],
    )
    complete = int(primary["complete_three_study_mapping"].sum())
    external_mapped = int(external["external_welch_p_value"].notna().sum())
    external_concordant = int(
        external.loc[external["external_welch_p_value"].notna(), "direction_concordant"].sum()
    )
    external_significant = int(external["external_fdr_lt_0_05"].sum())
    tub_mapped = int(tub["welch_p_value"].notna().sum())
    tub_concordant = int(tub.loc[tub["welch_p_value"].notna(), "direction_concordant"].sum())
    tub_significant = int(tub["fdr_lt_0_05"].sum())
    tub_pathway_significant = int(tub_pathways["fdr_bh_9_pathways"].lt(0.05).sum())
    paired_median_r = float(paired["pearson_r"].median())
    black_total = int(wgcna["in_black_module"].sum())
    black_external = int((wgcna["in_black_module"] & wgcna["external_fdr_lt_0_05"].fillna(False)).sum())
    external_outside_black = int((~wgcna["in_black_module"] & wgcna["external_fdr_lt_0_05"].fillna(False)).sum())
    outlier_count = int(sample_outliers["predefined_outlier_flag"].sum())
    outlier_sensitivity_significant = int(outlier_sensitivity["fdr_lt_0_05"].sum())
    text = f"""# Focused-Universe Reanalysis of Complement, Coagulation, and Extracellular-Matrix Transcriptional Programs in Diabetic Kidney Disease

Lizhen Deng

College of Life Science and Technology, Huazhong University of Science and Technology, Wuhan, Hubei, China

Corresponding author: Lizhen Deng (3070116993@qq.com)

ORCID: 0009-0003-2428-8176

## Abstract

Outcome-dependent gene filtering can invalidate conventional confirmation claims when discovery datasets are reused in meta-analysis. We reanalyzed nine biologically defined gene sets comprising 145 unique genes in human diabetic kidney disease transcriptomes. Discovery direction was estimated from GSE142025 and GSE96804. A selection-independent three-study synthesis added the prespecified glomerular GSE30528 compartment of source study GSE30122 and used restricted maximum likelihood with modified Hartung–Knapp inference, prediction intervals, and Benjamini–Hochberg correction across all 145 genes. Complete three-study estimates were available for {complete}/145 genes; none met meta-analysis FDR<0.05 (0/145). In the single external source study, {external_mapped}/145 genes mapped, {external_concordant}/{external_mapped} followed the discovery direction, and {external_significant}/145 met external Welch-test FDR<0.05. Sample-label permutation, which preserved gene correlation, supported 8/9 predefined pathways after pathway-level FDR correction; leukocyte adhesion did not meet the threshold. These results support pathway-level transcriptional concordance within one external glomerular study but do not establish independent multi-study confirmation, causality, pathway activity, or cell-type localization. The analysis illustrates how broader, selection-independent families and small-sample meta-analytic inference materially alter conclusions drawn from a filtered gene subset.

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

Complete Hedges' g estimates from GSE142025, GSE96804, and GSE30528 were available for {complete}/145 genes. Random-effects models used restricted maximum likelihood and modified Hartung–Knapp t inference with two degrees of freedom. Non-estimable members were retained as P=1 for correction over the complete 145-gene family. No gene met meta-analysis FDR<0.05 (0/145; Fig. 2 and Supplementary Table S1). The unmodified HKSJ intervals and both unmodified and modified prediction intervals are supplied as sensitivities. Prediction intervals were wide because only three source studies contributed, underscoring that a stable mean effect across future settings could not be established for individual genes.

### Gene-level evaluation in one external source study

In GSE30528 glomeruli, {external_mapped}/145 focused genes mapped to the platform and had estimable Welch tests. Of these, {external_concordant}/{external_mapped} followed the discovery direction, and {external_significant}/145 met FDR<0.05 after correction over the complete focused family (Fig. 3 and Supplementary Table S2). These are external results relative to the GSE142025/GSE96804 discovery direction, but they derive from one source study and therefore do not establish replication across independent source studies.

### Donor-disjoint tubulointerstitial sensitivity

Nine cross-accession donor pairs were inferred from sample titles and the shared GSE30122 parent study. Removing the corresponding GSE30529 samples left 5 DKD and 8 control tubulointerstitial samples. In this donor-disjoint sensitivity, {tub_mapped}/145 genes mapped, {tub_concordant}/{tub_mapped} followed the discovery direction, and {tub_significant}/145 met within-compartment FDR<0.05. Sample-label permutation supported {tub_pathway_significant}/9 pathways after pathway-level FDR correction (Supplementary Tables S9–S11). These results demonstrate non-uniform gene-level behavior across compartments: 7 genes met FDR<0.05 in both compartments, 32 only in GSE30528, and 13 only in donor-disjoint GSE30529. They do not convert the two compartments into independent studies.

### Pathway-level sample-label permutation

To retain gene-gene correlation, disease labels were permuted among the 22 GSE30528 samples while preserving the 9-DKD versus 13-control allocation, and all genes were recomputed jointly in each of 10,000 permutations. Eight of nine pathways met pathway-level FDR<0.05 (Fig. 4; Table 1). Complement activation, endothelial activation, macrophage recruitment, fibrosis, extracellular-matrix remodeling, coagulation, hypoxia, and transforming-growth-factor-beta signaling were supported. Leukocyte adhesion was directionally positive but did not meet the adjusted threshold (FDR=0.0998).

### Unsupervised structure, outlier, and WGCNA audits

PCA and sample-correlation diagnostics were run separately within each primary dataset. The predefined median-correlation rule flagged {outlier_count} samples; because all flagged samples were DKD samples and their biological versus technical origin could not be distinguished, they were retained in the primary analysis. Excluding them in a sensitivity analysis did not change the individual-gene conclusion ({outlier_sensitivity_significant}/145 meta-analysis FDR<0.05; Supplementary Figs. S2–S3 and Supplementary Tables S14–S16).

WGCNA was retained only as an internal GSE96804 description. Of the 145 focused genes, {black_total} mapped to the DKD-associated black module. Nineteen of these also met GSE30528 FDR<0.05, whereas 39 black-module genes did not; conversely, {external_outside_black} GSE30528 FDR-positive genes were outside the black module or outside the WGCNA top-variance universe (Supplementary Table S12). This partial overlap is compatible with compartment, cohort, variance-filter, and network-construction differences and is not independent validation.

### Table 1. External pathway-level permutation results

{pathway_table}

## Discussion

The principal finding is narrower than the earlier filtered-subset interpretation. When the full 145-gene family was analyzed with REML and modified Hartung–Knapp inference, no individual gene met meta-analysis FDR<0.05. The result does not imply that all underlying expression differences are absent. It shows that the available three-study evidence is insufficient for multiplicity-adjusted, small-sample random-effects claims at individual-gene resolution.

The independent evidence is stronger at pathway level within GSE30528. Eight predefined pathways were concordant under a sample-label permutation that retained the observed cross-gene structure. This supports an exploratory interpretation in which complement, coagulation, vascular, inflammatory, hypoxia, transforming-growth-factor-beta, and extracellular-matrix programs recur in the glomerular compartment of one external DKD study. Prior human kidney studies provide biological context for these themes [2–5], but the present reanalysis does not establish pathway activation, protein abundance, upstream regulation, or causality.

The contrast between 0/145 meta-analytic gene discoveries and 39/145 external Welch-test discoveries is not contradictory. The former uses random-effects inference over only three heterogeneous source studies with a t reference distribution and full-family multiplicity control. The latter tests one external compartment and does not estimate between-study heterogeneity. The pathway permutation then aggregates direction-aligned effects and gains stability while preserving gene correlation. These estimands answer different questions and should not be combined into a single “confirmed gene” count.

Model multiplicity was addressed by fixing one primary model—REML with modified Hartung–Knapp inference—and one 145-gene correction family. Unmodified HKSJ intervals, prediction intervals, sex adjustment, and correlation-outlier exclusion are labeled sensitivities and were not searched to choose a favorable result. A Bonferroni division by the six historical meta-analysis variants is therefore not applied to the new primary model; the historical variants no longer constitute six co-primary opportunities for discovery.

GSE30528 and GSE30529 should not be presented as separate sources of validation. Both are compartments of GSE30122, and nine cross-accession donor pairs can be inferred from sample titles and the shared parent study. Restricting the primary external estimand to GSE30528 avoids introducing an unestimated cross-compartment covariance. The earlier filtered 52-gene analysis reported 39/52 direction concordance in donor-disjoint GSE30529 but did not meet its matched-background criterion. That statistic is retained only as historical context because the candidate set and null construction were outcome-dependent; it is not combined with the new 145-gene tests.

The new full-family sensitivity clarifies, but does not erase, the compartment difference. Donor removal reduced GSE30529 from 10 DKD and 12 control samples to 5 and 8, respectively, making single-gene estimates less precise and more sensitive to individual samples. The two microdissected compartments also contain different cell mixtures and may carry genuinely different DKD responses: glomerular endothelial, mesangial, and podocyte programs need not track tubulointerstitial epithelial, stromal, vascular, and immune signals gene by gene. Manual microdissection, residual vascular content, probe behavior, and incomplete clinical metadata provide additional technical or compositional explanations. Consistent with this heterogeneity, only 7 genes met full-family FDR<0.05 in both compartments, whereas 32 were GSE30528-only and 13 were donor-disjoint-GSE30529-only. At pathway level, however, 8/9 sets were supported in each compartment, although the unsupported set differed. This pattern suggests that broad biological programs may recur while their contributing genes and effect magnitudes vary by compartment. It does not prove a compartment-specific mechanism, because donor pairing was inferred from titles, covariates were unavailable, and the retained GSE30529 sample was small. The correct conclusion is therefore same-study, compartment-sensitive concordance—not universal DKD replication and not evidence from two independent validation cohorts.

Across the nine paired donors, the median gene-wise cross-compartment expression correlation was only {paired_median_r:.3f}. This descriptive quantity is not the correlation between cohort-level disease-effect estimates and therefore cannot identify the rho parameter used in the retired correlated-compartment meta-analysis. Rather than choosing rho from these nine inferred pairs, the primary analysis removed rho dependence by representing GSE30122 with one prespecified compartment. The former rho-grid models are retained only as historical sensitivity analyses and are not used to select genes or conclusions.

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

The analysis archived the focused universe, per-study effects, discovery directions, complete random-effects output, external gene tests, pathway permutation results, donor-disjoint sensitivities, sample-level metadata, diagnostics, analysis settings, input hashes, and an output manifest. The execution environment used Python 3.13.9, NumPy 2.3.5, pandas 2.3.3, SciPy 1.16.3, statsmodels 0.14.5, and Matplotlib 3.10.6. M18 recomputed statistical results from frozen processed matrices and cached normalized GEO SOFT values; it did not independently repeat every upstream raw-read or CEL-file preprocessing step. A public immutable repository release and archival DOI remain required before submission.

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

Processed effect-size tables, pathway-permutation results, and analysis manifests require a public archival record with an immutable DOI before submission.

## Code availability

Analysis code, environment specifications, source-data tables, and machine-readable manifests require a public version-controlled release and a permanent archival DOI before submission.

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
"""
    return textwrap.dedent(text).strip() + "\n"


def build_cover_letter() -> str:
    return textwrap.dedent(
        """
        Dear Editors,

        Please consider our manuscript, “Focused-Universe Reanalysis of Complement, Coagulation, and Extracellular-Matrix Transcriptional Programs in Diabetic Kidney Disease,” for publication as an exploratory secondary analysis.

        The manuscript directly addresses two inferential risks identified during internal review. First, we replaced conventional inference over an outcome-filtered 61-gene subset with a selection-independent analysis of the complete 145-gene focused universe. Second, we treat GSE30528 and GSE30529 as two compartments of one source study rather than as independent validation cohorts. The revised primary synthesis uses REML with modified Hartung–Knapp inference, prediction intervals, and FDR control across all 145 genes.

        The revised result is deliberately bounded. No individual gene met meta-analysis FDR<0.05 across the complete family. In one external glomerular study, 94/141 mapped genes followed the discovery direction and 39/145 met external FDR<0.05. Sample-label permutation preserving gene correlation supported eight of nine predefined pathways. We interpret these findings as pathway-level concordance within one external study, not multi-study confirmation. Cell-type localization claims have been removed from the title, abstract, and principal conclusions.

        All datasets are public. The submission will include complete source tables, a machine-readable figure/table manifest, input hashes, environment specifications, and public code and processed-data links. Author details, the final repository URL, and archival DOI must be inserted before this letter is submitted.

        Sincerely,

        CORRESPONDING AUTHOR NAME REQUIRED
        AFFILIATION REQUIRED
        EMAIL REQUIRED
        """
    ).strip() + "\n"


def build_supplement(primary: pd.DataFrame, external: pd.DataFrame, pathways: pd.DataFrame) -> str:
    return textwrap.dedent(
        f"""
        # Supplementary Information

        **Manuscript title:** Focused-Universe Reanalysis of Complement, Coagulation, and Extracellular-Matrix Transcriptional Programs in Diabetic Kidney Disease

        **Authors:** AUTHOR LIST REQUIRED BEFORE SUBMISSION

        ## Supplementary Figure S1

        Complete three-study REML estimates with modified Hartung–Knapp 95% confidence intervals for the {int(primary['complete_three_study_mapping'].sum())} genes with complete mapping.

        ## Supplementary Table S1

        Complete 145-gene random-effects results, including REML tau-squared, unmodified and modified Hartung–Knapp uncertainty, prediction intervals, and 145-family FDR.

        ## Supplementary Table S2

        External GSE30528 gene-level results for all 145 genes, including mapping status, discovery direction, external effect, Welch P value, and 145-family FDR.

        ## Supplementary Table S3

        Nine-pathway sample-label permutation results. All genes were recomputed jointly in each of 10,000 label permutations.

        ## Supplementary Table S4

        Per-study Hedges' g estimates for GSE142025 C1, GSE96804, GSE30528, and secondary same-study compartment GSE30529.

        ## Supplementary Table S5

        Complete 145-row focused-universe lineage table, including the prior outcome-dependent prioritization status. This table is supplied for transparency and is not an inferential result.

        ## Supplementary Figure S2

        Unsupervised PCA diagnostics for GSE142025 C1, GSE96804, and GSE30528 using the 1,000 most variable genes within each dataset. These plots do not identify batch effects in the absence of batch labels.

        ## Supplementary Table S6

        Sample-metadata availability audit for batch, age, sex, renal-function, proteinuria, albuminuria, and medication fields.

        ## Supplementary Table S7

        GSE96804 sex-adjustment sensitivity for the complete 145-gene family.

        ## Supplementary Table S8

        PCA variance and disease-group diagnostic summaries for the three primary datasets.

        ## Supplementary Table S9

        Donor-disjoint GSE30529 gene-level sensitivity for the complete 145-gene family after removal of nine inferred shared-donor samples.

        ## Supplementary Table S10

        Donor-disjoint GSE30529 pathway sample-label permutation results using a fixed 5-DKD versus 8-control allocation.

        ## Supplementary Table S11

        Gene-level significance-switch summary comparing GSE30528 with donor-disjoint GSE30529. The `display_color_code` field provides the requested red/blue/purple/gray quick index without embedding spreadsheet-only formatting.

        ## Supplementary Table S12

        Descriptive WGCNA integration for all 145 focused genes. WGCNA is internal to GSE96804 and is not external validation.

        ## Supplementary Table S13

        Gene-wise expression correlations across nine title-inferred cross-compartment donor pairs. These correlations do not estimate covariance between cohort-level disease effects.

        ## Supplementary Table S14

        Complete sample-level metadata and analysis-role table for the primary datasets and the GSE30529 sensitivity.

        ## Supplementary Figure S3

        Within-dataset sample-correlation heatmaps based on the 1,000 most variable genes.

        ## Supplementary Table S15

        Sample-level median-correlation diagnostics and predefined robust-z outlier flags.

        ## Supplementary Table S16

        Complete 145-gene primary meta-analysis sensitivity after excluding the flagged samples.

        ## Supplementary Table S17

        Source data for the Figure 2 fixed-gene forest display, including each study estimate, confidence interval, random-effects weight, and pooled modified Hartung–Knapp interval.

        ## Supplementary Table S18

        Pairwise within-dataset sample correlations underlying Supplementary Figure S3.
        """
    ).strip() + "\n"


def write_manifest(records: list[dict[str, object]]) -> None:
    frame = pd.DataFrame(records)
    frame.to_csv(PACKAGE_DIR / "figure_table_manifest.csv", index=False)
    (PACKAGE_DIR / "figure_table_manifest.json").write_text(
        json.dumps(records, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def main() -> None:
    for directory in (FIGURE_DIR, MANUSCRIPT_DIR, PACKAGE_DIR, MAIN_FIGURE_DIR, SUPPLEMENT_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    primary = pd.read_csv(TABLE_DIR / "primary_reml_hk_145.csv")
    per_cohort = pd.read_csv(TABLE_DIR / "per_cohort_effect_sizes_145.csv")
    external = pd.read_csv(TABLE_DIR / "external_gse30528_gene_validation_145.csv")
    pathways = pd.read_csv(TABLE_DIR / "external_gse30528_pathway_permutation.csv")

    figure_outputs = {
        "Figure_1": figure_1_design(),
        "Figure_2": figure_2_primary(primary, per_cohort),
        "Figure_3": figure_3_external(external),
        "Figure_4": figure_4_pathways(pathways),
        "Supplementary_Figure_S1": supplementary_figure_1(primary),
    }
    pca_png = FIGURE_DIR / "Supplementary_Figure_S2_pca_diagnostics.png"
    pca_pdf = FIGURE_DIR / "Supplementary_Figure_S2_pca_diagnostics.pdf"
    if pca_png.exists() and pca_pdf.exists():
        figure_outputs["Supplementary_Figure_S2"] = [pca_png, pca_pdf]
    correlation_png = FIGURE_DIR / "Supplementary_Figure_S3_sample_correlations.png"
    correlation_pdf = FIGURE_DIR / "Supplementary_Figure_S3_sample_correlations.pdf"
    if correlation_png.exists() and correlation_pdf.exists():
        figure_outputs["Supplementary_Figure_S3"] = [correlation_png, correlation_pdf]
    for item_id, outputs in figure_outputs.items():
        target_dir = SUPPLEMENT_DIR if item_id.startswith("Supplementary") else MAIN_FIGURE_DIR
        for output in outputs:
            suffix = output.suffix
            target = target_dir / f"{item_id}{suffix}"
            shutil.copy2(output, target)

    manuscript = build_manuscript(primary, external, pathways)
    manuscript_path = MANUSCRIPT_DIR / "scientific_reports_remediated_manuscript.md"
    manuscript_path.write_text(manuscript, encoding="utf-8")
    shutil.copy2(manuscript_path, PACKAGE_DIR / manuscript_path.name)

    cover = build_cover_letter()
    (PACKAGE_DIR / "cover_letter_required_fields_template.md").write_text(cover, encoding="utf-8")
    supplement = build_supplement(primary, external, pathways)
    (SUPPLEMENT_DIR / "supplementary_information.md").write_text(supplement, encoding="utf-8")

    table_sources = {
        "Supplementary_Table_S1.csv": TABLE_DIR / "primary_reml_hk_145.csv",
        "Supplementary_Table_S2.csv": TABLE_DIR / "external_gse30528_gene_validation_145.csv",
        "Supplementary_Table_S3.csv": TABLE_DIR / "external_gse30528_pathway_permutation.csv",
        "Supplementary_Table_S4.csv": TABLE_DIR / "per_cohort_effect_sizes_145.csv",
        "Supplementary_Table_S5.csv": ROOT / "tables" / "phase3" / "core_candidate_gene_table.csv",
        "Supplementary_Table_S6.csv": TABLE_DIR / "sample_metadata_availability.csv",
        "Supplementary_Table_S7.csv": TABLE_DIR / "gse96804_sex_adjustment_sensitivity.csv",
        "Supplementary_Table_S8.csv": TABLE_DIR / "pca_diagnostic_summary.csv",
        "Supplementary_Table_S9.csv": TABLE_DIR / "gse30529_donor_disjoint_gene_validation_145.csv",
        "Supplementary_Table_S10.csv": TABLE_DIR / "gse30529_donor_disjoint_pathway_permutation.csv",
        "Supplementary_Table_S11.csv": TABLE_DIR / "compartment_significance_switch_summary_145.csv",
        "Supplementary_Table_S12.csv": TABLE_DIR / "wgcna_descriptive_integration_145.csv",
        "Supplementary_Table_S13.csv": TABLE_DIR / "paired_compartment_expression_correlations_145.csv",
        "Supplementary_Table_S14.csv": TABLE_DIR / "complete_sample_metadata.csv",
        "Supplementary_Table_S15.csv": TABLE_DIR / "sample_correlation_outlier_diagnostics.csv",
        "Supplementary_Table_S16.csv": TABLE_DIR / "primary_outlier_exclusion_sensitivity_145.csv",
        "Supplementary_Table_S17.csv": TABLE_DIR / "forest_plot_display_genes.csv",
        "Supplementary_Table_S18.csv": TABLE_DIR / "sample_pairwise_correlations.csv",
    }
    for name, source in table_sources.items():
        shutil.copy2(source, SUPPLEMENT_DIR / name)

    author_required = textwrap.dedent(
        """
        # Required author-supplied fields before submission

        - Full author names in publication order
        - Institutional affiliations and postal addresses
        - Corresponding author name and email
        - ORCID identifiers
        - CRediT author-contribution statement
        - Public code repository URL
        - Immutable release URL and archival DOI
        - Public processed-data/source-data archive URL and DOI
        - Author confirmation of funding and competing-interest statements
        - Final human visual approval of every manuscript and supplementary PDF page
        """
    ).strip() + "\n"
    (PACKAGE_DIR / "AUTHOR_AND_RELEASE_FIELDS_REQUIRED.md").write_text(author_required, encoding="utf-8")

    records = [
        {
            "item_id": "Figure_1", "title": "Analysis boundaries used for statistical remediation",
            "file_png": "main_figures/Figure_1.png", "file_pdf": "main_figures/Figure_1.pdf",
            "source_data": "config/phase3_focus_gene_sets.gmt; tables/stage21_m18_remediation/analysis_settings.json",
            "generator": SCRIPT_PATH.relative_to(ROOT).as_posix() + ":figure_1_design",
            "expected_panels": 1, "semantic_keywords": "145 genes; discovery; GSE30528; REML; modified Hartung-Knapp",
            "cited_in": "Results: Analysis family and study roles; Figure legend 1",
        },
        {
            "item_id": "Figure_2", "title": "Selection-independent three-study synthesis",
            "file_png": "main_figures/Figure_2.png", "file_pdf": "main_figures/Figure_2.pdf",
            "source_data": "tables/stage21_m18_remediation/primary_reml_hk_145.csv; tables/stage21_m18_remediation/per_cohort_effect_sizes_145.csv; tables/stage21_m18_remediation/forest_plot_display_genes.csv",
            "generator": SCRIPT_PATH.relative_to(ROOT).as_posix() + ":figure_2_primary",
            "expected_panels": 2, "semantic_keywords": "REML; modified Hartung-Knapp; 0/145; study effects; forest plot; weights",
            "cited_in": "Results: Selection-independent three-study synthesis; Figure legend 2",
        },
        {
            "item_id": "Figure_3", "title": "Gene-level evaluation in one external source study",
            "file_png": "main_figures/Figure_3.png", "file_pdf": "main_figures/Figure_3.pdf",
            "source_data": "tables/stage21_m18_remediation/external_gse30528_gene_validation_145.csv",
            "generator": SCRIPT_PATH.relative_to(ROOT).as_posix() + ":figure_3_external",
            "expected_panels": 1, "semantic_keywords": "GSE30528; 141 mapped; 94 concordant; 39 FDR",
            "cited_in": "Results: Gene-level evaluation; Figure legend 3",
        },
        {
            "item_id": "Figure_4", "title": "Pathway-level sample-label permutation",
            "file_png": "main_figures/Figure_4.png", "file_pdf": "main_figures/Figure_4.pdf",
            "source_data": "tables/stage21_m18_remediation/external_gse30528_pathway_permutation.csv",
            "generator": SCRIPT_PATH.relative_to(ROOT).as_posix() + ":figure_4_pathways",
            "expected_panels": 1, "semantic_keywords": "9 pathways; 10000 permutations; 8 FDR significant",
            "cited_in": "Results: Pathway-level sample-label permutation; Figure legend 4",
        },
        {
            "item_id": "Supplementary_Figure_S1", "title": "Complete three-study intervals",
            "file_png": "supplementary/Supplementary_Figure_S1.png",
            "file_pdf": "supplementary/Supplementary_Figure_S1.pdf",
            "source_data": "tables/stage21_m18_remediation/primary_reml_hk_145.csv",
            "generator": SCRIPT_PATH.relative_to(ROOT).as_posix() + ":supplementary_figure_1",
            "expected_panels": 1, "semantic_keywords": "125 genes; modified Hartung-Knapp confidence intervals",
            "cited_in": "Supplementary legend S1",
        },
    ]
    if "Supplementary_Figure_S2" in figure_outputs:
        records.append(
            {
                "item_id": "Supplementary_Figure_S2", "title": "Unsupervised PCA diagnostics",
                "file_png": "supplementary/Supplementary_Figure_S2.png",
                "file_pdf": "supplementary/Supplementary_Figure_S2.pdf",
                "source_data": "tables/stage21_m18_remediation/pca_diagnostic_summary.csv",
                "generator": "scripts/stage21_m18_remediation/run_m18_qc.py:pca_audit",
                "expected_panels": 3, "semantic_keywords": "GSE142025; GSE96804; GSE30528; PCA; 1000 variable genes",
                "cited_in": "Supplementary Information: Figure S2",
            }
        )
    if "Supplementary_Figure_S3" in figure_outputs:
        records.append(
            {
                "item_id": "Supplementary_Figure_S3", "title": "Within-dataset sample correlations",
                "file_png": "supplementary/Supplementary_Figure_S3.png",
                "file_pdf": "supplementary/Supplementary_Figure_S3.pdf",
                "source_data": "tables/stage21_m18_remediation/sample_pairwise_correlations.csv; tables/stage21_m18_remediation/sample_correlation_outlier_diagnostics.csv",
                "generator": "scripts/stage21_m18_remediation/run_m18_qc.py:sample_correlation_audit",
                "expected_panels": 3, "semantic_keywords": "sample correlations; outlier audit; robust z",
                "cited_in": "Supplementary Information: Figure S3",
            }
        )
    for index, (name, source) in enumerate(table_sources.items(), 1):
        records.append(
            {
                "item_id": f"Supplementary_Table_S{index}", "title": name.removesuffix(".csv").replace("_", " "),
                "file_png": "", "file_pdf": f"supplementary/{name}",
                "source_data": source.relative_to(ROOT).as_posix(),
                "generator": "verbatim source-table copy",
                "expected_panels": 0, "semantic_keywords": name.replace("_", " "),
                "cited_in": f"Supplementary Information: Table S{index}",
            }
        )
    write_manifest(records)

    package_files = [path for path in PACKAGE_DIR.rglob("*") if path.is_file()]
    package_manifest = pd.DataFrame(
        [
            {
                "path": path.relative_to(PACKAGE_DIR).as_posix(),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in sorted(package_files)
            if path.name != "package_file_manifest.csv"
        ]
    )
    package_manifest.to_csv(PACKAGE_DIR / "package_file_manifest.csv", index=False)
    print(f"manuscript={manuscript_path}")
    print(f"package_files={len(package_manifest)}")


if __name__ == "__main__":
    main()
