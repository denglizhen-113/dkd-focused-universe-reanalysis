from __future__ import annotations

import importlib.util
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import t as student_t
from scipy.stats import ttest_ind


ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = ROOT / "tables" / "stage21_m18_remediation"
DOC_DIR = ROOT / "docs" / "stage21_m18_remediation"
FIGURE_DIR = ROOT / "figures" / "stage21_m18_remediation"
M8_PATH = ROOT / "scripts" / "stage21_m8_validation" / "run_m8_analysis.py"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def bh(p: pd.Series) -> pd.Series:
    values = pd.to_numeric(p, errors="coerce").fillna(1.0).to_numpy(float)
    n = len(values)
    order = np.argsort(values, kind="mergesort")
    ranked = values[order]
    adjusted = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adjusted, 1.0)
    return pd.Series(out, index=p.index)


def pca_coordinates(expression: pd.DataFrame, top_variable: int = 1000) -> tuple[pd.DataFrame, list[float]]:
    matrix = expression.apply(pd.to_numeric, errors="coerce")
    variance = matrix.var(axis=1, ddof=1).sort_values(ascending=False)
    genes = variance.head(min(top_variable, len(variance))).index
    sample_by_gene = matrix.loc[genes].T
    sample_by_gene = sample_by_gene.fillna(sample_by_gene.mean(axis=0))
    centered = sample_by_gene - sample_by_gene.mean(axis=0)
    scales = sample_by_gene.std(axis=0, ddof=1).replace(0, 1)
    standardized = centered / scales
    u, singular, _ = np.linalg.svd(standardized.to_numpy(float), full_matrices=False)
    scores = u[:, :2] * singular[:2]
    explained = np.square(singular) / np.square(singular).sum()
    frame = pd.DataFrame(scores, index=sample_by_gene.index, columns=["PC1", "PC2"])
    frame.index.name = "sample_id"
    return frame.reset_index(), [float(explained[0]), float(explained[1])]


def prepare_gse30528(m8) -> tuple[pd.DataFrame, pd.DataFrame]:
    soft = ROOT / "data_raw" / "geo_soft" / "GSE30528" / "GSE30528_family.soft.gz"
    probe_matrix, probe_to_gene, _ = m8.parse_geo_soft(soft)
    selected = m8.select_highest_mean_probes(probe_matrix, probe_to_gene)
    expression = selected.drop(columns=["selected_probe_id", "selected_probe_mean"])
    metadata = pd.read_csv(ROOT / "data_raw" / "geo_soft" / "GSE30528" / "samples.csv")
    metadata["sample_id"] = metadata["gsm_accession"].astype(str)
    metadata["disease_group_clean"] = np.where(
        metadata["characteristics"].str.contains("diabetic kidney disease", case=False, na=False),
        "DKD", "Control"
    )
    return expression, metadata


def primary_dataset_inputs(m8) -> list[tuple[str, pd.DataFrame, pd.DataFrame]]:
    datasets: list[tuple[str, pd.DataFrame, pd.DataFrame]] = []
    expression_142 = m8.load_clean_expression(
        ROOT / "data_processed" / "GSE142025" / "expression_matrix_clean.csv"
    )
    metadata_142 = pd.read_csv(
        ROOT / "data_processed" / "GSE142025" / "sample_annotation_clean.csv"
    )
    metadata_142 = metadata_142.loc[
        metadata_142["disease_group_clean"].isin(["Advanced_DKD", "Control"])
    ].copy()
    expression_142 = expression_142.loc[:, metadata_142["sample_id"].tolist()]
    datasets.append(("GSE142025 C1", expression_142, metadata_142))

    expression_96804 = m8.load_clean_expression(
        ROOT / "data_processed" / "GSE96804" / "expression_matrix_clean.csv"
    )
    metadata_96804 = pd.read_csv(
        ROOT / "data_processed" / "GSE96804" / "sample_annotation_clean.csv"
    )
    datasets.append(("GSE96804", expression_96804, metadata_96804))

    expression_30528, metadata_30528 = prepare_gse30528(m8)
    datasets.append(("GSE30528", expression_30528, metadata_30528))
    return datasets


def complete_sample_metadata() -> pd.DataFrame:
    columns = [
        "dataset", "sample_id", "title", "platform", "tissue_region",
        "disease_group_clean", "sex", "age", "batch", "eGFR", "proteinuria",
        "albuminuria", "medication", "donor_id_inferred", "primary_analysis_role",
        "expression_input_type", "normalization_or_processing", "metadata_limit",
    ]
    rows = []
    for dataset in ("GSE142025", "GSE96804"):
        metadata = pd.read_csv(ROOT / "data_processed" / dataset / "sample_annotation_clean.csv")
        if dataset == "GSE142025":
            metadata = metadata.loc[
                metadata["disease_group_clean"].isin(["Advanced_DKD", "Control"])
            ].copy()
            role = "primary synthesis: C1 advanced DKD versus control"
            processing = "log2-transformed quantile-normalized expression"
        else:
            role = "primary synthesis and discovery direction"
            processing = "log2-scale normalized microarray expression"
        for item in metadata.itertuples(index=False):
            rows.append(
                {
                    "dataset": dataset,
                    "sample_id": str(item.sample_id),
                    "title": str(item.title),
                    "platform": str(item.platform),
                    "tissue_region": str(item.tissue_region),
                    "disease_group_clean": str(item.disease_group_clean),
                    "sex": getattr(item, "sex", np.nan),
                    "age": getattr(item, "age", np.nan),
                    "batch": getattr(item, "batch", np.nan),
                    "eGFR": getattr(item, "eGFR", np.nan),
                    "proteinuria": getattr(item, "proteinuria", np.nan),
                    "albuminuria": getattr(item, "albuminuria", np.nan),
                    "medication": getattr(item, "medication", np.nan),
                    "donor_id_inferred": "",
                    "primary_analysis_role": role,
                    "expression_input_type": "project frozen processed matrix derived from public expression data",
                    "normalization_or_processing": processing,
                    "metadata_limit": "Only fields exposed in archived public annotation are reported",
                }
            )
    shared = {"62", "67", "164", "168", "178", "76", "77", "81", "82"}
    for dataset, tissue, role in (
        ("GSE30528", "glomeruli", "primary external compartment and primary synthesis"),
        ("GSE30529", "tubulointerstitium", "secondary same-source compartment sensitivity"),
    ):
        metadata = pd.read_csv(ROOT / "data_raw" / "geo_soft" / dataset / "samples.csv")
        for item in metadata.itertuples(index=False):
            donor_match = re.search(r"Kidney\s+(\d+)", str(item.title), flags=re.IGNORECASE)
            donor = donor_match.group(1) if donor_match else ""
            group = "DKD" if "diabetic kidney disease" in str(item.characteristics).lower() else "Control"
            sample_role = role
            if dataset == "GSE30529" and donor in shared:
                sample_role += "; excluded from donor-disjoint sensitivity"
            elif dataset == "GSE30529":
                sample_role += "; retained in donor-disjoint sensitivity"
            rows.append(
                {
                    "dataset": dataset,
                    "sample_id": str(item.gsm_accession),
                    "title": str(item.title),
                    "platform": "GPL571",
                    "tissue_region": tissue,
                    "disease_group_clean": group,
                    "sex": np.nan,
                    "age": np.nan,
                    "batch": np.nan,
                    "eGFR": np.nan,
                    "proteinuria": np.nan,
                    "albuminuria": np.nan,
                    "medication": np.nan,
                    "donor_id_inferred": donor,
                    "primary_analysis_role": sample_role,
                    "expression_input_type": "cached processed GEO SOFT series-matrix values",
                    "normalization_or_processing": "GEO-supplied normalized values; label-independent highest-all-sample-mean probe",
                    "metadata_limit": "Cross-accession donor identity inferred from titles; not an explicit shared identifier",
                }
            )
    return pd.DataFrame(rows, columns=columns)


def sample_correlation_audit(m8):
    datasets = primary_dataset_inputs(m8)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    sample_rows = []
    pair_rows = []
    for ax, (dataset, expression, metadata) in zip(axes, datasets):
        numeric = expression.apply(pd.to_numeric, errors="coerce")
        variance = numeric.var(axis=1, ddof=1).sort_values(ascending=False)
        genes = variance.head(min(1000, len(variance))).index
        matrix = numeric.loc[genes].T
        matrix = matrix.fillna(matrix.mean(axis=0))
        correlation = matrix.T.corr(method="pearson")
        off_diagonal = correlation.copy()
        np.fill_diagonal(off_diagonal.values, np.nan)
        median_correlation = off_diagonal.median(axis=1, skipna=True)
        center = float(median_correlation.median())
        mad = float(np.median(np.abs(median_correlation.to_numpy(float) - center)))
        robust_z = (
            (median_correlation - center) / (1.4826 * mad)
            if mad > 0 else pd.Series(0.0, index=median_correlation.index)
        )
        group_map = metadata.set_index("sample_id")["disease_group_clean"]
        for sample in correlation.index:
            sample_rows.append(
                {
                    "dataset": dataset,
                    "sample_id": sample,
                    "disease_group": group_map.get(sample, "UNKNOWN"),
                    "genes_used": len(genes),
                    "median_sample_correlation": float(median_correlation.loc[sample]),
                    "robust_z_median_correlation": float(robust_z.loc[sample]),
                    "predefined_outlier_flag": bool(robust_z.loc[sample] < -3.0),
                    "outlier_rule": "robust z of median Pearson correlation below -3",
                }
            )
        for index, first in enumerate(correlation.index):
            for second in correlation.index[index + 1:]:
                pair_rows.append(
                    {
                        "dataset": dataset,
                        "sample_1": first,
                        "sample_2": second,
                        "pearson_correlation": float(correlation.loc[first, second]),
                    }
                )
        image = ax.imshow(
            correlation.to_numpy(float),
            norm=TwoSlopeNorm(vmin=-0.6, vcenter=0.0, vmax=1.0),
            cmap="coolwarm",
            aspect="auto",
        )
        ax.set_title(f"{dataset}\n{len(correlation)} samples")
        ax.set_xlabel("Samples")
        ax.set_ylabel("Samples")
        ax.set_xticks([])
        ax.set_yticks([])
    fig.colorbar(image, ax=axes.ravel().tolist(), shrink=0.78, label="Pearson correlation")
    fig.suptitle("Within-dataset sample correlations using the 1,000 most variable genes", y=0.98)
    fig.subplots_adjust(wspace=0.18, top=0.84, right=0.93)
    return pd.DataFrame(sample_rows), pd.DataFrame(pair_rows), fig


def metadata_audit() -> pd.DataFrame:
    rows = []
    for dataset in ("GSE142025", "GSE96804"):
        metadata = pd.read_csv(ROOT / "data_processed" / dataset / "sample_annotation_clean.csv")
        row = {"dataset": dataset, "n_samples": len(metadata)}
        for field in ("batch", "age", "sex", "eGFR", "proteinuria", "albuminuria", "medication"):
            series = metadata[field].replace("", np.nan) if field in metadata else pd.Series(dtype=object)
            row[f"{field}_nonmissing"] = int(series.notna().sum())
            row[f"{field}_unique_nonmissing"] = int(series.dropna().astype(str).nunique())
        rows.append(row)
    for dataset in ("GSE30528", "GSE30529"):
        metadata = pd.read_csv(ROOT / "data_raw" / "geo_soft" / dataset / "samples.csv")
        row = {"dataset": dataset, "n_samples": len(metadata)}
        for field in ("batch", "age", "sex", "eGFR", "proteinuria", "albuminuria", "medication"):
            row[f"{field}_nonmissing"] = 0
            row[f"{field}_unique_nonmissing"] = 0
        rows.append(row)
    return pd.DataFrame(rows)


def sex_adjustment_gse96804() -> pd.DataFrame:
    universe = pd.read_csv(TABLE_DIR / "focused_universe_145.csv")
    genes = universe["gene_symbol"].astype(str).str.upper().tolist()
    expression = pd.read_csv(
        ROOT / "data_processed" / "GSE96804" / "expression_matrix_clean.csv"
    ).set_index("gene_symbol")
    expression.index = expression.index.astype(str).str.upper()
    metadata = pd.read_csv(
        ROOT / "data_processed" / "GSE96804" / "sample_annotation_clean.csv"
    ).set_index("sample_id")
    samples = [sample for sample in expression.columns if sample in metadata.index]
    metadata = metadata.loc[samples]
    disease = metadata["disease_group_clean"].eq("DKD").astype(float).to_numpy()
    sex_text = metadata["sex"].astype(str).str.lower()
    valid_samples = sex_text.isin(["male", "female"]).to_numpy()
    samples = list(np.asarray(samples)[valid_samples])
    disease = disease[valid_samples]
    male = sex_text.loc[samples].eq("male").astype(float).to_numpy()
    x = np.column_stack([np.ones(len(samples)), disease, male])
    inverse = np.linalg.inv(x.T @ x)
    rows = []
    for gene in genes:
        if gene not in expression.index:
            rows.append({"gene_symbol": gene, "mapping_status": "NOT_MAPPED"})
            continue
        y = expression.loc[gene, samples].to_numpy(float)
        beta = inverse @ x.T @ y
        residual = y - x @ beta
        df = len(y) - x.shape[1]
        sigma2 = float(residual @ residual / df)
        se = math.sqrt(sigma2 * inverse[1, 1])
        p = float(2 * student_t.sf(abs(beta[1] / se), df)) if se > 0 else math.nan
        raw_dkd = y[disease == 1]
        raw_control = y[disease == 0]
        raw_diff = float(np.mean(raw_dkd) - np.mean(raw_control))
        raw_p = float(ttest_ind(raw_dkd, raw_control, equal_var=False).pvalue)
        rows.append(
            {
                "gene_symbol": gene,
                "mapping_status": "MAPPED",
                "n_samples": len(y),
                "n_dkd": int(disease.sum()),
                "n_control": int((1 - disease).sum()),
                "n_male": int(male.sum()),
                "n_female": int((1 - male).sum()),
                "unadjusted_mean_difference": raw_diff,
                "unadjusted_welch_p": raw_p,
                "sex_adjusted_disease_beta": float(beta[1]),
                "sex_adjusted_disease_se": se,
                "sex_adjusted_disease_p": p,
                "direction_changed_after_sex_adjustment": np.sign(raw_diff) != np.sign(beta[1]),
            }
        )
    frame = pd.DataFrame(rows)
    frame["sex_adjusted_fdr_bh_145"] = bh(frame["sex_adjusted_disease_p"])
    return frame


def pca_audit(m8) -> tuple[pd.DataFrame, plt.Figure]:
    datasets = primary_dataset_inputs(m8)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6))
    rows = []
    colors = {"Advanced_DKD": "#d7301f", "DKD": "#d7301f", "Control": "#3182bd"}
    for ax, (dataset, expression, metadata) in zip(axes, datasets):
        coordinates, explained = pca_coordinates(expression)
        group_map = metadata.set_index("sample_id")["disease_group_clean"]
        coordinates["group"] = coordinates["sample_id"].map(group_map)
        for group, group_frame in coordinates.groupby("group"):
            ax.scatter(
                group_frame["PC1"], group_frame["PC2"], label=group,
                c=colors.get(group, "#777777"), s=32, alpha=0.82,
                edgecolor="white", linewidth=0.35
            )
        groups = list(coordinates["group"].dropna().unique())
        pc1_p = math.nan
        pc2_p = math.nan
        if len(groups) == 2:
            first = coordinates.loc[coordinates["group"].eq(groups[0])]
            second = coordinates.loc[coordinates["group"].eq(groups[1])]
            pc1_p = float(ttest_ind(first["PC1"], second["PC1"], equal_var=False).pvalue)
            pc2_p = float(ttest_ind(first["PC2"], second["PC2"], equal_var=False).pvalue)
        rows.append(
            {
                "dataset": dataset,
                "n_samples": len(coordinates),
                "genes_used": min(1000, len(expression)),
                "pc1_variance_percent": explained[0] * 100,
                "pc2_variance_percent": explained[1] * 100,
                "pc1_group_welch_p": pc1_p,
                "pc2_group_welch_p": pc2_p,
                "interpretation": "Unsupervised diagnostic; group separation may reflect biology, confounding, or both",
            }
        )
        ax.set_title(dataset)
        ax.set_xlabel(f"PC1 ({explained[0] * 100:.1f}%)")
        ax.set_ylabel(f"PC2 ({explained[1] * 100:.1f}%)")
        ax.legend(frameon=False, fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Unsupervised PCA of the 1,000 most variable genes within each primary dataset", y=1.02)
    fig.tight_layout()
    return pd.DataFrame(rows), fig


def main() -> None:
    for directory in (TABLE_DIR, DOC_DIR, FIGURE_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    m8 = load_module("stage21_m8_for_qc", M8_PATH)
    availability = metadata_audit()
    sex_adjustment = sex_adjustment_gse96804()
    pca, figure = pca_audit(m8)
    complete_metadata = complete_sample_metadata()
    sample_correlation, pairwise_correlation, correlation_figure = sample_correlation_audit(m8)
    availability.to_csv(TABLE_DIR / "sample_metadata_availability.csv", index=False)
    sex_adjustment.to_csv(TABLE_DIR / "gse96804_sex_adjustment_sensitivity.csv", index=False)
    pca.to_csv(TABLE_DIR / "pca_diagnostic_summary.csv", index=False)
    complete_metadata.to_csv(TABLE_DIR / "complete_sample_metadata.csv", index=False)
    sample_correlation.to_csv(TABLE_DIR / "sample_correlation_outlier_diagnostics.csv", index=False)
    pairwise_correlation.to_csv(TABLE_DIR / "sample_pairwise_correlations.csv", index=False)
    figure.savefig(FIGURE_DIR / "Supplementary_Figure_S2_pca_diagnostics.png", dpi=400, bbox_inches="tight")
    figure.savefig(FIGURE_DIR / "Supplementary_Figure_S2_pca_diagnostics.pdf", bbox_inches="tight")
    plt.close(figure)
    correlation_figure.savefig(
        FIGURE_DIR / "Supplementary_Figure_S3_sample_correlations.png", dpi=400, bbox_inches="tight"
    )
    correlation_figure.savefig(
        FIGURE_DIR / "Supplementary_Figure_S3_sample_correlations.pdf", bbox_inches="tight"
    )
    plt.close(correlation_figure)

    mapped = sex_adjustment["mapping_status"].eq("MAPPED")
    changes = int(sex_adjustment.loc[mapped, "direction_changed_after_sex_adjustment"].sum())
    adjusted_significant = int(sex_adjustment["sex_adjusted_fdr_bh_145"].lt(0.05).sum())
    correlation_outliers = int(sample_correlation["predefined_outlier_flag"].sum())
    report = f"""# Covariate and unsupervised-structure audit

Public sample annotations do not supply a common adjustment set across the three primary studies. Batch and age are unavailable for every primary dataset. Sex is complete in GSE96804 but unavailable in GSE142025 and GSE30528. Renal-function, proteinuria, albuminuria, and medication fields are absent or unusable for the primary contrasts. A common multivariable sensitivity analysis is therefore not estimable from the archived public metadata.

Within GSE96804, sex adjustment was feasible for the complete focused family. Among mapped genes, {changes} changed disease-effect direction after adding sex to an intercept-plus-disease linear model; {adjusted_significant}/145 met sex-adjusted BH FDR<0.05. This is a dataset-specific sensitivity and does not repair missing covariates in the other studies.

Unsupervised PCA was run separately within GSE142025 C1, GSE96804, and GSE30528 using the 1,000 most variable genes. The plots are diagnostics, not batch corrections. Any disease-group separation may reflect disease biology, composition, unrecorded technical structure, or confounding. Because no batch labels were available, PCA cannot identify or quantify batch-specific effects.

Within-dataset sample-correlation heatmaps used the same 1,000-variable-gene rule. A sample was predefined as a potential correlation outlier when the robust z score of its median Pearson correlation was below -3; {correlation_outliers} samples met that rule. This audit can identify unusual global profiles, but it cannot label their cause as biological or technical without batch and processing metadata.
"""
    (DOC_DIR / "M18_COVARIATE_AND_PCA_AUDIT.md").write_text(report, encoding="utf-8")
    print(
        f"sex_direction_changes={changes} adjusted_fdr={adjusted_significant} "
        f"sample_correlation_outliers={correlation_outliers} metadata_rows={len(complete_metadata)}"
    )


if __name__ == "__main__":
    main()
