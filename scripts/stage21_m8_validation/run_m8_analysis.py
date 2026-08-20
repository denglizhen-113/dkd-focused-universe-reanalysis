#!/usr/bin/env python3
"""Stage 21-M8 method correction and 61-gene validation.

Execution date: 2026-07-29
Random seed: 20260729 (no inferential result uses random sampling)

Absolute inputs:
- C:/Users/denglizhen/Documents/凝血通讯轴/tables/phase3_5_limma_audit/core_candidate_gene_table_limma_audited.csv
- C:/Users/denglizhen/Documents/凝血通讯轴/data_processed/GSE142025/expression_matrix_clean.csv
- C:/Users/denglizhen/Documents/凝血通讯轴/data_processed/GSE142025/sample_annotation_clean.csv
- C:/Users/denglizhen/Documents/凝血通讯轴/data_processed/GSE96804/expression_matrix_clean.csv
- C:/Users/denglizhen/Documents/凝血通讯轴/data_processed/GSE96804/sample_annotation_clean.csv
- C:/Users/denglizhen/Documents/凝血通讯轴/data_raw/geo_soft/GSE30528/GSE30528_family.soft.gz
- C:/Users/denglizhen/Documents/凝血通讯轴/data_raw/geo_soft/GSE30528/samples.csv
- C:/Users/denglizhen/Documents/凝血通讯轴/data_raw/geo_soft/GSE30529/GSE30529_family.soft.gz
- C:/Users/denglizhen/Documents/凝血通讯轴/data_raw/geo_soft/GSE30529/samples.csv
- C:/Users/denglizhen/Documents/凝血通讯轴/data_raw/geo_soft/GSE111154/GSE111154_family.soft.gz
- C:/Users/denglizhen/Documents/凝血通讯轴/data_raw/geo_soft/GSE111154/samples.csv
- C:/Users/denglizhen/Documents/凝血通讯轴/tables/stage10_scireports_evidence_expansion/external_validation_gene_level_results.csv
- C:/Users/denglizhen/Documents/凝血通讯轴/tables/stage12_high_confidence_upgrade/gene_level_random_effects_meta_analysis.csv

Absolute output directories:
- C:/Users/denglizhen/Documents/凝血通讯轴/tables/stage21_m8_validation/
- C:/Users/denglizhen/Documents/凝血通讯轴/docs/stage21_m8_validation/
- C:/Users/denglizhen/Documents/凝血通讯轴/logs/stage21_m8_validation/
- C:/Users/denglizhen/Documents/凝血通讯轴/scripts/stage21_m8_validation/

Output files: contrast_C1_results.csv, contrast_C2_results.csv,
C2_vs_stage10_reproduction_check.csv, GSE30528_61gene_validation.csv,
GSE30529_61gene_validation.csv, direction_concordance_test.csv,
effect_size_correlation.csv, per_cohort_effect_sizes.csv,
meta_analysis_unified_scale_61genes.csv, stage12_vs_m8_comparison.csv,
input_file_manifest.csv, analysis_settings.json, five requested Markdown reports,
requirements.txt, and m8_run.log.

No Stage 10, Stage 12, Stage 14, manuscript, WGCNA, single-nucleus, or
CellPhoneDB artifact is modified.
"""

from __future__ import annotations

import gzip
import hashlib
import inspect
import json
import logging
import math
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import binomtest, chi2, norm, spearmanr, ttest_ind
from statsmodels.stats.meta_analysis import effectsize_smd
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = Path(__file__).resolve()
TABLE_DIR = ROOT / "tables" / "stage21_m8_validation"
DOC_DIR = ROOT / "docs" / "stage21_m8_validation"
LOG_DIR = ROOT / "logs" / "stage21_m8_validation"
RANDOM_SEED = 20260729
MAIN_GENES = ["C3", "COL1A1", "COL1A2", "FN1", "PECAM1", "PLG", "VCAM1", "VWF"]
SIGNATURES = {
    "core_remodeling": ["C3", "COL1A1", "COL1A2", "VCAM1"],
    "ECM": ["COL1A1", "COL1A2", "FN1"],
    "complement_coagulation": ["C3", "VWF"],
    "vascular_tubular_context": ["VCAM1", "VWF", "PECAM1"],
}


class UTCFormatter(logging.Formatter):
    """Format log record timestamps as UTC when the rendered suffix is Z."""

    converter = time.gmtime


def source_line(function: Any) -> int:
    """Return the current source line for an analysis function."""

    return int(inspect.getsourcelines(function)[1])


def bh_adjust(p_values: np.ndarray | pd.Series | list[float]) -> np.ndarray:
    """Apply BH to finite values while retaining explicit missing rows."""

    values = np.asarray(p_values, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    if valid.any():
        adjusted[valid] = multipletests(values[valid], method="fdr_bh")[1]
    return adjusted


def select_highest_mean_probes(
    probe_matrix: pd.DataFrame, probe_to_gene: pd.Series
) -> pd.DataFrame:
    """Select one probe per gene by the Stage 10 all-sample mean rule."""

    common = probe_matrix.index.intersection(probe_to_gene.index)
    work = probe_matrix.loc[common].copy()
    work.insert(0, "gene_symbol", probe_to_gene.loc[common].astype(str).str.upper())
    work.insert(1, "selected_probe_id", work.index.astype(str))
    sample_columns = list(probe_matrix.columns)
    work["selected_probe_mean"] = work[sample_columns].mean(axis=1)
    work = work.sort_values(
        ["gene_symbol", "selected_probe_mean", "selected_probe_id"],
        ascending=[True, False, True],
        kind="mergesort",
    )
    selected = work.drop_duplicates("gene_symbol", keep="first").set_index("gene_symbol")
    return selected[["selected_probe_id", "selected_probe_mean", *sample_columns]]


def hedges_g_from_samples(case: np.ndarray, control: np.ndarray) -> dict[str, float]:
    """Return Hedges' g and sampling variance from two independent groups."""

    case_values = np.asarray(case, dtype=float)
    control_values = np.asarray(control, dtype=float)
    case_values = case_values[np.isfinite(case_values)]
    control_values = control_values[np.isfinite(control_values)]
    if len(case_values) < 2 or len(control_values) < 2:
        return {"hedges_g": np.nan, "variance": np.nan, "se": np.nan}
    case_sd = case_values.std(ddof=1)
    control_sd = control_values.std(ddof=1)
    if not np.isfinite(case_sd) or not np.isfinite(control_sd):
        return {"hedges_g": np.nan, "variance": np.nan, "se": np.nan}
    pooled_numerator = (len(case_values) - 1) * case_sd**2 + (
        len(control_values) - 1
    ) * control_sd**2
    if pooled_numerator <= 0:
        return {"hedges_g": np.nan, "variance": np.nan, "se": np.nan}
    effect, variance = effectsize_smd(
        case_values.mean(),
        case_sd,
        len(case_values),
        control_values.mean(),
        control_sd,
        len(control_values),
    )
    return {
        "hedges_g": float(effect),
        "variance": float(variance),
        "se": float(math.sqrt(variance)),
    }


def dersimonian_laird(
    effects: np.ndarray | pd.Series | list[float],
    variances: np.ndarray | pd.Series | list[float],
) -> dict[str, float]:
    """DerSimonian-Laird random-effects synthesis with heterogeneity tests."""

    yi = np.asarray(effects, dtype=float)
    vi = np.asarray(variances, dtype=float)
    valid = np.isfinite(yi) & np.isfinite(vi) & (vi > 0)
    yi = yi[valid]
    vi = vi[valid]
    if len(yi) < 2:
        return {
            key: np.nan
            for key in (
                "pooled_effect",
                "pooled_se",
                "ci_low",
                "ci_high",
                "p_value",
                "q",
                "q_p_value",
                "i2",
                "tau2",
                "k",
            )
        }
    fixed_weights = 1.0 / vi
    fixed_effect = float(np.sum(fixed_weights * yi) / np.sum(fixed_weights))
    q_value = float(np.sum(fixed_weights * (yi - fixed_effect) ** 2))
    df = len(yi) - 1
    c_value = float(
        np.sum(fixed_weights)
        - np.sum(fixed_weights**2) / np.sum(fixed_weights)
    )
    tau2 = max(0.0, (q_value - df) / c_value) if c_value > 0 else 0.0
    random_weights = 1.0 / (vi + tau2)
    pooled = float(np.sum(random_weights * yi) / np.sum(random_weights))
    pooled_se = float(math.sqrt(1.0 / np.sum(random_weights)))
    z_value = pooled / pooled_se
    p_value = float(2.0 * norm.sf(abs(z_value)))
    ci_low = float(pooled - 1.959963984540054 * pooled_se)
    ci_high = float(pooled + 1.959963984540054 * pooled_se)
    i2 = max(0.0, (q_value - df) / q_value) * 100.0 if q_value > 0 else 0.0
    return {
        "pooled_effect": pooled,
        "pooled_se": pooled_se,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "p_value": p_value,
        "q": q_value,
        "q_p_value": float(chi2.sf(q_value, df)),
        "i2": float(i2),
        "tau2": float(tau2),
        "k": float(len(yi)),
    }


def direction_binom_test(successes: int, total: int) -> dict[str, Any]:
    """One-sided exact binomial direction-concordance test."""

    if total <= 0 or successes < 0 or successes > total:
        return {"successes": successes, "total": total, "p_value": np.nan}
    test = binomtest(successes, total, p=0.5, alternative="greater")
    return {
        "successes": int(successes),
        "total": int(total),
        "p_value": float(test.pvalue),
    }


def build_gse142025_contrasts(
    metadata: pd.DataFrame,
) -> dict[str, dict[str, list[str]]]:
    """Build the two frozen GSE142025 contrast memberships."""

    required = {"sample_id", "disease_group_clean"}
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f"Missing metadata columns: {sorted(missing)}")
    if metadata["sample_id"].duplicated().any():
        raise ValueError("GSE142025 sample_id values must be unique")
    advanced = metadata.loc[
        metadata["disease_group_clean"].eq("Advanced_DKD"), "sample_id"
    ].astype(str).tolist()
    early = metadata.loc[
        metadata["disease_group_clean"].eq("Early_DKD"), "sample_id"
    ].astype(str).tolist()
    control = metadata.loc[
        metadata["disease_group_clean"].eq("Control"), "sample_id"
    ].astype(str).tolist()
    return {
        "C1_advanced_vs_control": {"case": advanced, "control": control},
        "C2_advanced_vs_early_plus_control": {
            "case": advanced,
            "control": [*early, *control],
        },
    }


def load_candidates(path: Path) -> pd.DataFrame:
    """Load and validate the frozen 61-gene denominator."""

    candidates = pd.read_csv(path)
    required = {
        "gene_symbol",
        "theme",
        "advanced_logFC_limma",
        "GSE96804_logFC_limma",
    }
    missing = required.difference(candidates.columns)
    if missing:
        raise ValueError(f"Candidate table is missing columns: {sorted(missing)}")
    if len(candidates) != 61 or candidates["gene_symbol"].nunique() != 61:
        raise ValueError("The frozen candidate denominator must contain 61 unique genes")
    candidates = candidates.copy()
    candidates["gene_symbol"] = candidates["gene_symbol"].astype(str).str.upper()
    candidates["candidate_source_row"] = np.arange(2, len(candidates) + 2)
    advanced = pd.to_numeric(candidates["advanced_logFC_limma"], errors="coerce")
    gse96804 = pd.to_numeric(candidates["GSE96804_logFC_limma"], errors="coerce")
    candidates["discovery_effect"] = advanced
    candidates["discovery_direction"] = np.where(
        advanced > 0, "UP", np.where(advanced < 0, "DOWN", "ZERO")
    )
    candidates["discovery_direction_matches_gse96804"] = np.sign(advanced).eq(
        np.sign(gse96804)
    )
    return candidates


def load_clean_expression(path: Path) -> pd.DataFrame:
    """Load a frozen clean gene-by-sample expression matrix."""

    matrix = pd.read_csv(path)
    if "gene_symbol" not in matrix.columns:
        raise ValueError(f"Missing gene_symbol column in {path}")
    matrix["gene_symbol"] = matrix["gene_symbol"].astype(str).str.upper()
    if matrix["gene_symbol"].duplicated().any():
        duplicated = matrix.loc[matrix["gene_symbol"].duplicated(), "gene_symbol"].tolist()
        raise ValueError(f"Duplicate gene symbols in {path}: {duplicated[:10]}")
    matrix = matrix.set_index("gene_symbol")
    return matrix.apply(pd.to_numeric, errors="coerce")


def parse_geo_soft(path: Path) -> tuple[pd.DataFrame, pd.Series, list[str]]:
    """Parse Stage 10-compatible platform symbols and sample expression tables."""

    probe_to_gene: dict[str, str] = {}
    sample_values: dict[str, pd.Series] = {}
    current_sample: str | None = None
    platform_headers: list[str] = []
    with gzip.open(path, "rt", errors="ignore") as handle:
        iterator = iter(handle)
        for raw_line in iterator:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE ="):
                current_sample = line.split("=", 1)[1].strip()
            elif line == "!platform_table_begin":
                platform_headers = next(iterator).rstrip("\n").split("\t")
                for platform_line in iterator:
                    platform_line = platform_line.rstrip("\n")
                    if platform_line == "!platform_table_end":
                        break
                    fields = platform_line.split("\t")
                    record = dict(zip(platform_headers, fields))
                    gene = record.get("Gene Symbol", "").split(" /// ")[0].strip().upper()
                    probe_id = record.get("ID", "").strip()
                    if probe_id and gene:
                        probe_to_gene[probe_id] = gene
            elif line == "!sample_table_begin":
                if current_sample is None:
                    raise ValueError(f"Sample table without sample identifier in {path}")
                sample_headers = next(iterator).rstrip("\n").split("\t")
                values: dict[str, float] = {}
                for sample_line in iterator:
                    sample_line = sample_line.rstrip("\n")
                    if sample_line == "!sample_table_end":
                        break
                    fields = sample_line.split("\t")
                    if len(fields) < 2:
                        continue
                    try:
                        values[fields[0]] = float(fields[1])
                    except ValueError:
                        continue
                sample_values[current_sample] = pd.Series(values, dtype=float)
    if not sample_values:
        raise ValueError(f"No sample expression tables parsed from {path}")
    probe_matrix = pd.DataFrame(sample_values)
    mapping = pd.Series(probe_to_gene, name="gene_symbol", dtype=str)
    mapping.index.name = "probe_id"
    return probe_matrix, mapping, platform_headers


def run_welch_contrast(
    expression: pd.DataFrame,
    case_ids: list[str],
    control_ids: list[str],
    candidates: pd.DataFrame,
    *,
    dataset: str,
    contrast_name: str,
    selected_probe_metadata: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Run a frozen 61-gene Welch contrast with explicit BH correction."""

    if set(case_ids).intersection(control_ids):
        raise ValueError(f"Overlapping contrast sides for {dataset} {contrast_name}")
    missing_samples = sorted(set([*case_ids, *control_ids]).difference(expression.columns))
    if missing_samples:
        raise ValueError(f"Expression matrix missing samples: {missing_samples}")
    matrix = expression.copy()
    matrix.index = matrix.index.astype(str).str.upper()
    source_rows = {gene: idx + 2 for idx, gene in enumerate(matrix.index)}
    rows: list[dict[str, Any]] = []
    analysis_line = source_line(run_welch_contrast)
    for candidate in candidates.itertuples(index=False):
        gene = str(candidate.gene_symbol).upper()
        base: dict[str, Any] = {
            "dataset": dataset,
            "contrast": contrast_name,
            "gene_symbol": gene,
            "theme": candidate.theme,
            "candidate_source_row": int(candidate.candidate_source_row),
            "discovery_effect": float(candidate.discovery_effect),
            "discovery_direction": candidate.discovery_direction,
            "n_case": len(case_ids),
            "n_control": len(control_ids),
            "case_ids": ";".join(case_ids),
            "control_ids": ";".join(control_ids),
            "test_method": "Welch two-sample t-test (equal_var=False; two-sided)",
            "multiple_testing": "Benjamini-Hochberg across frozen 61-gene family",
            "analysis_code_path": str(SCRIPT_PATH),
            "analysis_function": "run_welch_contrast",
            "analysis_code_line": analysis_line,
            "expression_source_row": source_rows.get(gene, np.nan),
            "selected_probe_id": "",
            "selected_probe_mean": np.nan,
        }
        if selected_probe_metadata is not None and gene in selected_probe_metadata.index:
            base["selected_probe_id"] = str(
                selected_probe_metadata.loc[gene, "selected_probe_id"]
            )
            base["selected_probe_mean"] = float(
                selected_probe_metadata.loc[gene, "selected_probe_mean"]
            )
        if gene not in matrix.index:
            base.update(
                {
                    "mapping_status": "NOT_MAPPED",
                    "case_mean": np.nan,
                    "control_mean": np.nan,
                    "case_sd": np.nan,
                    "control_sd": np.nan,
                    "log2FC": np.nan,
                    "p_value": np.nan,
                    "direction": "NOT_MAPPED",
                }
            )
            rows.append(base)
            continue
        case = pd.to_numeric(matrix.loc[gene, case_ids], errors="coerce").to_numpy(float)
        control = pd.to_numeric(
            matrix.loc[gene, control_ids], errors="coerce"
        ).to_numpy(float)
        test = ttest_ind(case, control, equal_var=False, nan_policy="omit")
        effect = float(np.nanmean(case) - np.nanmean(control))
        base.update(
            {
                "mapping_status": "MAPPED_TESTED",
                "case_mean": float(np.nanmean(case)),
                "control_mean": float(np.nanmean(control)),
                "case_sd": float(np.nanstd(case, ddof=1)),
                "control_sd": float(np.nanstd(control, ddof=1)),
                "log2FC": effect,
                "p_value": float(test.pvalue),
                "direction": "UP" if effect > 0 else "DOWN" if effect < 0 else "ZERO",
            }
        )
        rows.append(base)
    result = pd.DataFrame(rows)
    if len(result) != 61 or result["gene_symbol"].nunique() != 61:
        raise AssertionError("Welch output did not retain the frozen 61-gene denominator")
    result["adjusted_p_value"] = bh_adjust(result["p_value"].to_numpy(float))
    result["direction_matches_discovery"] = np.where(
        result["mapping_status"].eq("MAPPED_TESTED"),
        result["direction"].eq(result["discovery_direction"]),
        pd.NA,
    )
    return result


def fixed_effect_pool(
    effects: np.ndarray | pd.Series | list[float],
    variances: np.ndarray | pd.Series | list[float],
) -> dict[str, float]:
    """Inverse-variance fixed-effect pool used for combined validation direction."""

    yi = np.asarray(effects, dtype=float)
    vi = np.asarray(variances, dtype=float)
    valid = np.isfinite(yi) & np.isfinite(vi) & (vi > 0)
    yi = yi[valid]
    vi = vi[valid]
    if len(yi) == 0:
        return {"effect": np.nan, "se": np.nan, "k": 0.0}
    weights = 1.0 / vi
    return {
        "effect": float(np.sum(weights * yi) / np.sum(weights)),
        "se": float(math.sqrt(1.0 / np.sum(weights))),
        "k": float(len(yi)),
    }


def classify_conclusion_change(
    old_effect: float, old_fdr: float, new_effect: float, new_fdr: float
) -> str:
    """Classify direction and FDR<0.05 changes without minimizing either."""

    if not all(np.isfinite([old_effect, old_fdr, new_effect, new_fdr])):
        return "NOT_COMPARABLE"
    direction_reversal = np.sign(old_effect) != np.sign(new_effect)
    old_significant = old_fdr < 0.05
    new_significant = new_fdr < 0.05
    changes: list[str] = []
    if direction_reversal:
        changes.append("DIRECTION_REVERSAL")
    if old_significant and not new_significant:
        changes.append("SIGNIFICANCE_LOSS")
    elif not old_significant and new_significant:
        changes.append("SIGNIFICANCE_GAIN")
    return "_AND_".join(changes) if changes else "UNCHANGED"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def find_line(path: Path, needle: str) -> int:
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if needle in line:
            return line_number
    raise ValueError(f"Needle not found in {path}: {needle}")


def markdown_table(frame: pd.DataFrame, columns: list[str]) -> str:
    def clean(value: Any) -> str:
        if pd.isna(value):
            return "NA"
        return str(value).replace("|", "\\|").replace("\n", " ")

    header = "| " + " | ".join(columns) + " |"
    divider = "| " + " | ".join(["---"] * len(columns)) + " |"
    body = [
        "| " + " | ".join(clean(row[column]) for column in columns) + " |"
        for _, row in frame[columns].iterrows()
    ]
    return "\n".join([header, divider, *body])


def write_table(frame: pd.DataFrame, filename: str, logger: logging.Logger) -> Path:
    path = TABLE_DIR / filename
    frame.to_csv(path, index=False, encoding="utf-8", lineterminator="\n")
    logger.info("OUTPUT_TABLE path=%s rows=%d columns=%d", path, len(frame), len(frame.columns))
    return path


def write_report(text: str, filename: str, logger: logging.Logger) -> Path:
    path = DOC_DIR / filename
    path.write_text(text.rstrip() + "\n", encoding="utf-8")
    logger.info(
        "OUTPUT_REPORT path=%s lines=%d", path, len(path.read_text(encoding="utf-8").splitlines())
    )
    return path


def build_reproduction_check(
    c2_results: pd.DataFrame, historical: pd.DataFrame
) -> pd.DataFrame:
    stage10 = historical.loc[
        historical["dataset"].eq("GSE142025_advanced_vs_control")
    ][["gene_symbol", "logFC", "P.Value", "adj.P.Val"]].copy()
    stage10["stage10_result_source_row"] = stage10.index.to_numpy() + 2
    stage10 = stage10.rename(
        columns={
            "logFC": "stage10_log2FC",
            "P.Value": "stage10_p_value",
            "adj.P.Val": "stage10_adjusted_p_value_8gene_family",
        }
    )
    columns = [
        "gene_symbol",
        "m8_log2FC",
        "m8_p_value",
        "m8_adjusted_p_value_61gene_family",
        "candidate_source_row",
        "expression_source_row",
        "analysis_code_path",
        "analysis_function",
        "analysis_code_line",
    ]
    m8 = c2_results.rename(
        columns={
            "log2FC": "m8_log2FC",
            "p_value": "m8_p_value",
            "adjusted_p_value": "m8_adjusted_p_value_61gene_family",
        }
    )[columns]
    check = m8.merge(stage10, on="gene_symbol", how="left")
    check["historical_comparator_available"] = check["stage10_log2FC"].notna()
    check["absolute_log2FC_difference"] = (
        check["m8_log2FC"] - check["stage10_log2FC"]
    ).abs()
    check["absolute_raw_p_difference"] = (
        check["m8_p_value"] - check["stage10_p_value"]
    ).abs()
    check["absolute_adjusted_p_difference"] = (
        check["m8_adjusted_p_value_61gene_family"]
        - check["stage10_adjusted_p_value_8gene_family"]
    ).abs()
    check["effect_matches_1e_12"] = np.where(
        check["historical_comparator_available"],
        check["absolute_log2FC_difference"].le(1e-12),
        pd.NA,
    )
    check["raw_p_matches_1e_12"] = np.where(
        check["historical_comparator_available"],
        check["absolute_raw_p_difference"].le(1e-12),
        pd.NA,
    )
    check["adjusted_p_matches_1e_12"] = np.where(
        check["historical_comparator_available"],
        check["absolute_adjusted_p_difference"].le(1e-12),
        pd.NA,
    )
    check["reproduction_status"] = np.select(
        [
            ~check["historical_comparator_available"],
            check["effect_matches_1e_12"].eq(True)
            & check["raw_p_matches_1e_12"].eq(True)
            & check["adjusted_p_matches_1e_12"].eq(True),
            check["effect_matches_1e_12"].eq(True)
            & check["raw_p_matches_1e_12"].eq(True),
        ],
        [
            "NO_STAGE10_COMPARATOR",
            "EXACT_ALL_REPRODUCTION",
            "RAW_STATISTICS_REPRODUCED_ADJUSTED_P_DIFFERS_BY_TEST_FAMILY",
        ],
        default="RAW_STATISTIC_DIFFERENCE",
    )
    check["comparison_note"] = np.where(
        check["historical_comparator_available"],
        "Stage 10 BH family=8; M8 BH family=frozen 61",
        "Stage 10 did not calculate this candidate",
    )
    return check


def add_stage10_comparison(
    validation: pd.DataFrame, historical: pd.DataFrame, dataset: str
) -> pd.DataFrame:
    old = historical.loc[historical["dataset"].eq(dataset)][
        ["gene_symbol", "logFC", "P.Value", "adj.P.Val"]
    ].rename(
        columns={
            "logFC": "stage10_log2FC",
            "P.Value": "stage10_p_value",
            "adj.P.Val": "stage10_adjusted_p_value_8gene_family",
        }
    )
    result = validation.merge(old, on="gene_symbol", how="left")
    result["stage10_comparator_available"] = result["stage10_log2FC"].notna()
    result["stage10_effect_matches_1e_12"] = np.where(
        result["stage10_comparator_available"],
        (result["log2FC"] - result["stage10_log2FC"]).abs().le(1e-12),
        pd.NA,
    )
    result["stage10_raw_p_matches_1e_12"] = np.where(
        result["stage10_comparator_available"],
        (result["p_value"] - result["stage10_p_value"]).abs().le(1e-12),
        pd.NA,
    )
    finite_test_count = int(result["p_value"].notna().sum())
    result["stage10_adjusted_p_note"] = np.where(
        result["stage10_comparator_available"],
        f"Not expected to match: Stage 10 family=8; M8 finite-test family={finite_test_count} within frozen 61-row denominator",
        "No Stage 10 comparator",
    )
    return result


def build_per_cohort_effect_sizes(
    candidates: pd.DataFrame, cohort_definitions: list[dict[str, Any]]
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    analysis_line = source_line(build_per_cohort_effect_sizes)
    for cohort in cohort_definitions:
        matrix = cohort["expression"]
        matrix = matrix.copy()
        matrix.index = matrix.index.astype(str).str.upper()
        source_rows = {gene: index + 2 for index, gene in enumerate(matrix.index)}
        probe_metadata = cohort.get("probe_metadata")
        for candidate in candidates.itertuples(index=False):
            gene = candidate.gene_symbol
            base = {
                "cohort": cohort["cohort"],
                "gene_symbol": gene,
                "theme": candidate.theme,
                "candidate_source_row": int(candidate.candidate_source_row),
                "platform": cohort["platform"],
                "tissue": cohort["tissue"],
                "processing": cohort["processing"],
                "expression_scale": cohort["expression_scale"],
                "contrast": cohort["contrast"],
                "n_case": len(cohort["case_ids"]),
                "n_control": len(cohort["control_ids"]),
                "case_ids": ";".join(cohort["case_ids"]),
                "control_ids": ";".join(cohort["control_ids"]),
                "input_path": str(cohort["input_path"]),
                "expression_source_row": source_rows.get(gene, np.nan),
                "selected_probe_id": "",
                "effect_definition": "Hedges g standardized mean difference (case minus control)",
                "variance_definition": "statsmodels.effectsize_smd sampling variance",
                "analysis_code_path": str(SCRIPT_PATH),
                "analysis_function": "build_per_cohort_effect_sizes; hedges_g_from_samples",
                "analysis_code_line": analysis_line,
            }
            if probe_metadata is not None and gene in probe_metadata.index:
                base["selected_probe_id"] = str(
                    probe_metadata.loc[gene, "selected_probe_id"]
                )
            if gene not in matrix.index:
                base.update(
                    {
                        "mapping_status": "NOT_MAPPED",
                        "case_mean": np.nan,
                        "control_mean": np.nan,
                        "case_sd": np.nan,
                        "control_sd": np.nan,
                        "raw_mean_difference": np.nan,
                        "hedges_g": np.nan,
                        "variance": np.nan,
                        "se": np.nan,
                    }
                )
                rows.append(base)
                continue
            case = matrix.loc[gene, cohort["case_ids"]].to_numpy(float)
            control = matrix.loc[gene, cohort["control_ids"]].to_numpy(float)
            effect = hedges_g_from_samples(case, control)
            base.update(
                {
                    "mapping_status": (
                        "MAPPED_EFFECT_ESTIMATED"
                        if np.isfinite(effect["hedges_g"])
                        else "MAPPED_ZERO_OR_INVALID_VARIANCE"
                    ),
                    "case_mean": float(np.nanmean(case)),
                    "control_mean": float(np.nanmean(control)),
                    "case_sd": float(np.nanstd(case, ddof=1)),
                    "control_sd": float(np.nanstd(control, ddof=1)),
                    "raw_mean_difference": float(np.nanmean(case) - np.nanmean(control)),
                    **effect,
                }
            )
            rows.append(base)
    result = pd.DataFrame(rows)
    expected = len(candidates) * len(cohort_definitions)
    if len(result) != expected:
        raise AssertionError(f"Expected {expected} cohort-gene rows, found {len(result)}")
    return result


def build_meta_results(
    candidates: pd.DataFrame, per_cohort: pd.DataFrame, expected_cohorts: list[str]
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    analysis_line = source_line(build_meta_results)
    for candidate in candidates.itertuples(index=False):
        sub = per_cohort.loc[per_cohort["gene_symbol"].eq(candidate.gene_symbol)]
        valid = sub.loc[
            np.isfinite(sub["hedges_g"])
            & np.isfinite(sub["variance"])
            & sub["variance"].gt(0)
        ]
        missing = [cohort for cohort in expected_cohorts if cohort not in valid["cohort"].tolist()]
        meta = dersimonian_laird(valid["hedges_g"], valid["variance"])
        status = "META_ANALYZED" if len(valid) >= 2 else "NOT_META_ANALYZED"
        rows.append(
            {
                "gene_symbol": candidate.gene_symbol,
                "theme": candidate.theme,
                "candidate_source_row": int(candidate.candidate_source_row),
                "discovery_effect": float(candidate.discovery_effect),
                "discovery_direction": candidate.discovery_direction,
                "meta_status": status,
                "n_cohorts_expected": len(expected_cohorts),
                "n_cohorts_included": len(valid),
                "included_cohorts": ";".join(valid["cohort"].tolist()),
                "missing_or_nonestimable_cohorts": ";".join(missing),
                "pooled_hedges_g": meta["pooled_effect"],
                "pooled_se": meta["pooled_se"],
                "ci_95_low": meta["ci_low"],
                "ci_95_high": meta["ci_high"],
                "p_value": meta["p_value"],
                "cochran_q": meta["q"],
                "cochran_q_p_value": meta["q_p_value"],
                "i2_percent": meta["i2"],
                "tau2": meta["tau2"],
                "pooled_direction": (
                    "UP"
                    if meta["pooled_effect"] > 0
                    else "DOWN"
                    if meta["pooled_effect"] < 0
                    else "NOT_ESTIMATED"
                ),
                "model": "DerSimonian-Laird random effects",
                "effect_scale": "Hedges g standardized mean difference",
                "normal_ci_and_p": True,
                "analysis_code_path": str(SCRIPT_PATH),
                "analysis_function": "build_meta_results; dersimonian_laird",
                "analysis_code_line": analysis_line,
            }
        )
    result = pd.DataFrame(rows)
    result["adjusted_p_value"] = bh_adjust(result["p_value"].to_numpy(float))
    result["direction_matches_discovery"] = np.where(
        result["meta_status"].eq("META_ANALYZED"),
        result["pooled_direction"].eq(result["discovery_direction"]),
        pd.NA,
    )
    return result


def compare_stage12(stage12: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    old = stage12.copy()
    old["stage12_result_source_row"] = old.index.to_numpy() + 2
    old = old.rename(
        columns={
            "random_effect": "stage12_raw_scale_effect",
            "SE": "stage12_se",
            "P.Value": "stage12_p_value",
            "FDR": "stage12_fdr",
            "I2": "stage12_i2_percent",
            "tau2": "stage12_tau2",
        }
    )
    new = meta.rename(
        columns={
            "pooled_hedges_g": "m8_hedges_g",
            "pooled_se": "m8_se",
            "p_value": "m8_p_value",
            "adjusted_p_value": "m8_fdr",
            "i2_percent": "m8_i2_percent",
            "tau2": "m8_tau2",
        }
    )
    columns_old = [
        "gene_symbol",
        "stage12_raw_scale_effect",
        "stage12_se",
        "stage12_p_value",
        "stage12_fdr",
        "stage12_i2_percent",
        "stage12_tau2",
        "stage12_result_source_row",
    ]
    columns_new = [
        "gene_symbol",
        "m8_hedges_g",
        "m8_se",
        "m8_p_value",
        "m8_fdr",
        "m8_i2_percent",
        "m8_tau2",
        "ci_95_low",
        "ci_95_high",
        "candidate_source_row",
    ]
    comparison = old[columns_old].merge(new[columns_new], on="gene_symbol", how="left")
    comparison["stage12_direction"] = np.where(
        comparison["stage12_raw_scale_effect"] > 0, "UP", "DOWN"
    )
    comparison["m8_direction"] = np.where(comparison["m8_hedges_g"] > 0, "UP", "DOWN")
    comparison["direction_reversal"] = comparison["stage12_direction"].ne(
        comparison["m8_direction"]
    )
    comparison["stage12_fdr_lt_0_05"] = comparison["stage12_fdr"].lt(0.05)
    comparison["m8_fdr_lt_0_05"] = comparison["m8_fdr"].lt(0.05)
    comparison["conclusion_change"] = comparison.apply(
        lambda row: classify_conclusion_change(
            row["stage12_raw_scale_effect"],
            row["stage12_fdr"],
            row["m8_hedges_g"],
            row["m8_fdr"],
        ),
        axis=1,
    )
    comparison["scale_comparability_note"] = (
        "Stage12 raw-scale pooled effect and M8 Hedges g differ in units; compare direction and inference status, not magnitude"
    )
    comparison["analysis_code_path"] = str(SCRIPT_PATH)
    comparison["analysis_function"] = "compare_stage12; classify_conclusion_change"
    comparison["analysis_code_line"] = source_line(compare_stage12)
    return comparison


def build_validation_direction_and_correlation(
    candidates: pd.DataFrame,
    validations: dict[str, pd.DataFrame],
    per_cohort: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    direction_rows: list[dict[str, Any]] = []
    analysis_line = source_line(build_validation_direction_and_correlation)
    correlation_rows: list[dict[str, Any]] = []
    combined_rows: list[dict[str, Any]] = []

    for dataset, result in validations.items():
        mapped = result.loc[result["mapping_status"].eq("MAPPED_TESTED")].copy()
        success = int(mapped["direction_matches_discovery"].eq(True).sum())
        test = direction_binom_test(success, len(mapped))
        direction_rows.append(
            {
                "validation_scope": dataset,
                "combination_rule": "Single validation cohort",
                "mapped_genes": len(mapped),
                "direction_concordant_genes": success,
                "direction_discordant_genes": len(mapped) - success,
                "concordance_fraction": success / len(mapped) if len(mapped) else np.nan,
                "null_probability": 0.5,
                "alternative": "greater",
                "p_value": test["p_value"],
                "test_method": "one-sided exact binomial test",
                "dependence_caveat": "Candidate genes are biologically correlated; exact binomial independence is an approximation",
                "analysis_code_path": str(SCRIPT_PATH),
                "analysis_function": "build_validation_direction_and_correlation; direction_binom_test",
                "analysis_code_line": analysis_line,
            }
        )
        correlation_rows.extend(
            correlation_records(
                candidates,
                dataset,
                result.set_index("gene_symbol")["log2FC"],
                "validation log2-scale mean difference",
            )
        )

    validation_effects = per_cohort.loc[
        per_cohort["cohort"].isin(["GSE30528", "GSE30529"])
    ]
    for candidate in candidates.itertuples(index=False):
        sub = validation_effects.loc[
            validation_effects["gene_symbol"].eq(candidate.gene_symbol)
        ]
        fixed = fixed_effect_pool(sub["hedges_g"], sub["variance"])
        effect = fixed["effect"]
        direction = "UP" if effect > 0 else "DOWN" if effect < 0 else "NOT_ESTIMATED"
        combined_rows.append(
            {
                "gene_symbol": candidate.gene_symbol,
                "discovery_effect": float(candidate.discovery_effect),
                "discovery_direction": candidate.discovery_direction,
                "combined_validation_hedges_g": effect,
                "combined_validation_se": fixed["se"],
                "n_validation_cohorts": int(fixed["k"]),
                "combined_direction": direction,
                "direction_matches_discovery": (
                    direction == candidate.discovery_direction
                    if direction != "NOT_ESTIMATED"
                    else pd.NA
                ),
            }
        )
    combined = pd.DataFrame(combined_rows)
    combined_mapped = combined.loc[combined["combined_direction"].ne("NOT_ESTIMATED")]
    combined_success = int(
        combined_mapped["direction_matches_discovery"].eq(True).sum()
    )
    combined_test = direction_binom_test(combined_success, len(combined_mapped))
    direction_rows.append(
        {
            "validation_scope": "GSE30528_GSE30529_combined",
            "combination_rule": "Per-gene inverse-variance fixed-effect pool of validation Hedges g",
            "mapped_genes": len(combined_mapped),
            "direction_concordant_genes": combined_success,
            "direction_discordant_genes": len(combined_mapped) - combined_success,
            "concordance_fraction": (
                combined_success / len(combined_mapped) if len(combined_mapped) else np.nan
            ),
            "null_probability": 0.5,
            "alternative": "greater",
            "p_value": combined_test["p_value"],
            "test_method": "one-sided exact binomial test",
            "dependence_caveat": "Candidate genes are biologically correlated; exact binomial independence is an approximation",
            "analysis_code_path": str(SCRIPT_PATH),
            "analysis_function": "build_validation_direction_and_correlation; fixed_effect_pool; direction_binom_test",
            "analysis_code_line": analysis_line,
        }
    )
    correlation_rows.extend(
        correlation_records(
            candidates,
            "GSE30528_GSE30529_combined",
            combined.set_index("gene_symbol")["combined_validation_hedges_g"],
            "inverse-variance fixed-effect validation Hedges g",
        )
    )
    return pd.DataFrame(direction_rows), pd.DataFrame(correlation_rows), combined


def correlation_records(
    candidates: pd.DataFrame,
    dataset: str,
    validation_effect: pd.Series,
    validation_scale: str,
) -> list[dict[str, Any]]:
    points = candidates[["gene_symbol", "discovery_effect"]].copy()
    points["validation_effect"] = points["gene_symbol"].map(validation_effect)
    valid = points.dropna(subset=["discovery_effect", "validation_effect"])
    if len(valid) >= 3:
        correlation = spearmanr(
            valid["discovery_effect"], valid["validation_effect"], nan_policy="omit"
        )
        rho = float(correlation.statistic)
        p_value = float(correlation.pvalue)
    else:
        rho = np.nan
        p_value = np.nan
    common = {
        "validation_dataset": dataset,
        "spearman_rho": rho,
        "correlation_p_value": p_value,
        "n_mapped": len(valid),
        "correlation_method": "Spearman rank correlation (two-sided)",
        "selection_rationale": "Cross-platform effects have different raw scales; rank correlation was prespecified",
        "discovery_scale": "GSE142025 advanced_logFC_limma from frozen candidate table",
        "validation_scale": validation_scale,
        "analysis_code_path": str(SCRIPT_PATH),
        "analysis_function": "correlation_records; scipy.stats.spearmanr",
        "analysis_code_line": source_line(correlation_records),
    }
    records: list[dict[str, Any]] = [
        {
            "record_type": "SUMMARY",
            "gene_symbol": "",
            "discovery_effect": np.nan,
            "validation_effect": np.nan,
            "mapped": len(valid),
            **common,
        }
    ]
    for row in points.itertuples(index=False):
        records.append(
            {
                "record_type": "SCATTER_POINT",
                "gene_symbol": row.gene_symbol,
                "discovery_effect": row.discovery_effect,
                "validation_effect": row.validation_effect,
                "mapped": bool(np.isfinite(row.validation_effect)),
                **common,
            }
        )
    return records


def setup_logger() -> logging.Logger:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("stage21_m8")
    logger.handlers.clear()
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(LOG_DIR / "m8_run.log", mode="w", encoding="utf-8")
    handler.setFormatter(UTCFormatter("%(asctime)sZ | %(levelname)s | %(message)s"))
    logger.addHandler(handler)
    return logger


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    DOC_DIR.mkdir(parents=True, exist_ok=True)
    logger = setup_logger()
    np.random.seed(RANDOM_SEED)
    started = datetime.now(timezone.utc)
    logger.info("M8_START utc=%s seed=%d python=%s", started.isoformat(), RANDOM_SEED, sys.version)
    requirements_path = SCRIPT_PATH.parent / "requirements.txt"
    if requirements_path.exists():
        requirements_lines = len(requirements_path.read_text(encoding="utf-8").splitlines())
        logger.info(
            "REQUIREMENTS path=%s source=python_-m_pip_freeze lines=%d",
            requirements_path,
            requirements_lines,
        )
    else:
        logger.warning("REQUIREMENTS_MISSING path=%s", requirements_path)

    paths = {
        "candidates": ROOT / "tables" / "phase3_5_limma_audit" / "core_candidate_gene_table_limma_audited.csv",
        "gse142025_expression": ROOT / "data_processed" / "GSE142025" / "expression_matrix_clean.csv",
        "gse142025_metadata": ROOT / "data_processed" / "GSE142025" / "sample_annotation_clean.csv",
        "gse96804_expression": ROOT / "data_processed" / "GSE96804" / "expression_matrix_clean.csv",
        "gse96804_metadata": ROOT / "data_processed" / "GSE96804" / "sample_annotation_clean.csv",
        "gse30528_soft": ROOT / "data_raw" / "geo_soft" / "GSE30528" / "GSE30528_family.soft.gz",
        "gse30528_metadata": ROOT / "data_raw" / "geo_soft" / "GSE30528" / "samples.csv",
        "gse30529_soft": ROOT / "data_raw" / "geo_soft" / "GSE30529" / "GSE30529_family.soft.gz",
        "gse30529_metadata": ROOT / "data_raw" / "geo_soft" / "GSE30529" / "samples.csv",
        "gse111154_soft": ROOT / "data_raw" / "geo_soft" / "GSE111154" / "GSE111154_family.soft.gz",
        "gse111154_metadata": ROOT / "data_raw" / "geo_soft" / "GSE111154" / "samples.csv",
        "stage10_results": ROOT / "tables" / "stage10_scireports_evidence_expansion" / "external_validation_gene_level_results.csv",
        "stage10_script": ROOT / "scripts" / "stage10_scireports_evidence_expansion" / "run_stage10_external_validation.py",
        "stage12_results": ROOT / "tables" / "stage12_high_confidence_upgrade" / "gene_level_random_effects_meta_analysis.csv",
        "stage12_script": ROOT / "scripts" / "stage12_high_confidence_upgrade" / "run_stage12_upgrade.py",
        "stage14_results": ROOT / "tables" / "stage14_5to7if_route_rescue" / "GSE111154_pathology_stage_signature_results.csv",
        "stage14_script": ROOT / "scripts" / "stage14_5to7if_route_rescue" / "run_stage14_rescue_audit.py",
    }
    missing_inputs = [str(path) for path in paths.values() if not path.exists()]
    if missing_inputs:
        raise FileNotFoundError(f"Missing frozen inputs: {missing_inputs}")
    manifest_rows = []
    for label, path in paths.items():
        manifest_rows.append(
            {
                "input_id": label,
                "absolute_path": str(path.resolve()),
                "size_bytes": path.stat().st_size,
                "modified_time": datetime.fromtimestamp(
                    path.stat().st_mtime, timezone.utc
                ).isoformat(),
                "sha256": sha256_file(path),
            }
        )
    input_manifest = pd.DataFrame(manifest_rows)
    write_table(input_manifest, "input_file_manifest.csv", logger)

    candidates = load_candidates(paths["candidates"])
    logger.info(
        "CANDIDATE_DENOMINATOR rows=%d unique=%d direction_disagreements=%d",
        len(candidates),
        candidates["gene_symbol"].nunique(),
        int((~candidates["discovery_direction_matches_gse96804"]).sum()),
    )
    gse142025_expression = load_clean_expression(paths["gse142025_expression"])
    gse142025_metadata = pd.read_csv(paths["gse142025_metadata"])
    if set(gse142025_metadata["sample_id"]) != set(gse142025_expression.columns):
        raise AssertionError("GSE142025 metadata and expression sample IDs differ")
    contrasts = build_gse142025_contrasts(gse142025_metadata)
    c1_members = contrasts["C1_advanced_vs_control"]
    c2_members = contrasts["C2_advanced_vs_early_plus_control"]
    if (len(c1_members["case"]), len(c1_members["control"])) != (21, 9):
        raise AssertionError("C1 membership is not 21 Advanced versus 9 Control")
    if (len(c2_members["case"]), len(c2_members["control"])) != (21, 15):
        raise AssertionError("C2 membership is not 21 Advanced versus 15 Early+Control")
    c1 = run_welch_contrast(
        gse142025_expression,
        c1_members["case"],
        c1_members["control"],
        candidates,
        dataset="GSE142025",
        contrast_name="C1_advanced_vs_control",
    )
    c2 = run_welch_contrast(
        gse142025_expression,
        c2_members["case"],
        c2_members["control"],
        candidates,
        dataset="GSE142025",
        contrast_name="C2_advanced_vs_early_plus_control",
    )
    for frame, source in ((c1, paths["gse142025_expression"]), (c2, paths["gse142025_expression"])):
        frame["expression_input_path"] = str(source.resolve())
        frame["metadata_input_path"] = str(paths["gse142025_metadata"].resolve())
    stage10 = pd.read_csv(paths["stage10_results"])
    c2_check = build_reproduction_check(c2, stage10)
    c2_check["stage10_result_path"] = str(paths["stage10_results"].resolve())
    c2_check["comparison_code_path"] = str(SCRIPT_PATH)
    c2_check["comparison_function"] = "build_reproduction_check"
    c2_check["comparison_code_line"] = source_line(build_reproduction_check)
    write_table(c1, "contrast_C1_results.csv", logger)
    write_table(c2, "contrast_C2_results.csv", logger)
    write_table(c2_check, "C2_vs_stage10_reproduction_check.csv", logger)

    validation_results: dict[str, pd.DataFrame] = {}
    selected_matrices: dict[str, pd.DataFrame] = {}
    validation_members: dict[str, dict[str, list[str]]] = {}
    validation_platform_headers: dict[str, list[str]] = {}
    for accession in ("GSE30528", "GSE30529"):
        probe_matrix, probe_to_gene, platform_headers = parse_geo_soft(
            paths[f"{accession.lower()}_soft"]
        )
        selected = select_highest_mean_probes(probe_matrix, probe_to_gene)
        expression = selected.drop(columns=["selected_probe_id", "selected_probe_mean"])
        metadata = pd.read_csv(paths[f"{accession.lower()}_metadata"])
        is_dkd = metadata["characteristics"].str.contains(
            "diabetic kidney disease", case=False, na=False
        )
        case_ids = metadata.loc[is_dkd, "gsm_accession"].astype(str).tolist()
        control_ids = metadata.loc[~is_dkd, "gsm_accession"].astype(str).tolist()
        expected_counts = (9, 13) if accession == "GSE30528" else (10, 12)
        if (len(case_ids), len(control_ids)) != expected_counts:
            raise AssertionError(f"Unexpected {accession} groups")
        result = run_welch_contrast(
            expression,
            case_ids,
            control_ids,
            candidates,
            dataset=accession,
            contrast_name="DKD_vs_Control",
            selected_probe_metadata=selected[["selected_probe_id", "selected_probe_mean"]],
        )
        result = add_stage10_comparison(result, stage10, accession)
        finite_test_count = int(result["p_value"].notna().sum())
        result["multiple_testing"] = (
            f"Benjamini-Hochberg across {finite_test_count} finite P values; "
            "frozen 61-row denominator retained"
        )
        result["expression_input_path"] = str(paths[f"{accession.lower()}_soft"].resolve())
        result["metadata_input_path"] = str(paths[f"{accession.lower()}_metadata"].resolve())
        result["analysis_function"] = "parse_geo_soft; select_highest_mean_probes; run_welch_contrast"
        result["analysis_code_line"] = (
            f"parse_geo_soft:{source_line(parse_geo_soft)};"
            f"select_highest_mean_probes:{source_line(select_highest_mean_probes)};"
            f"run_welch_contrast:{source_line(run_welch_contrast)}"
        )
        validation_results[accession] = result
        selected_matrices[accession] = selected
        validation_members[accession] = {"case": case_ids, "control": control_ids}
        validation_platform_headers[accession] = platform_headers

    gse96804_expression = load_clean_expression(paths["gse96804_expression"])
    gse96804_metadata = pd.read_csv(paths["gse96804_metadata"])
    if set(gse96804_metadata["sample_id"]) != set(gse96804_expression.columns):
        raise AssertionError("GSE96804 metadata and expression sample IDs differ")
    gse96804_case = gse96804_metadata.loc[
        gse96804_metadata["disease_group_clean"].eq("DKD"), "sample_id"
    ].astype(str).tolist()
    gse96804_control = gse96804_metadata.loc[
        gse96804_metadata["disease_group_clean"].eq("Control"), "sample_id"
    ].astype(str).tolist()
    if (len(gse96804_case), len(gse96804_control)) != (41, 20):
        raise AssertionError("Unexpected GSE96804 group counts")

    cohort_definitions = [
        {
            "cohort": "GSE96804",
            "expression": gse96804_expression,
            "case_ids": gse96804_case,
            "control_ids": gse96804_control,
            "platform": "GPL17586; Affymetrix Human Transcriptome Array 2.0",
            "tissue": "glomeruli",
            "processing": "project frozen clean matrix; supplied normalized expression",
            "expression_scale": "log2-scale normalized microarray expression",
            "contrast": "DKD_vs_Control",
            "input_path": paths["gse96804_expression"].resolve(),
            "probe_metadata": None,
        },
        {
            "cohort": "GSE30528",
            "expression": selected_matrices["GSE30528"].drop(
                columns=["selected_probe_id", "selected_probe_mean"]
            ),
            "case_ids": validation_members["GSE30528"]["case"],
            "control_ids": validation_members["GSE30528"]["control"],
            "platform": "GPL571; Affymetrix Human Genome U133A 2.0",
            "tissue": "glomeruli",
            "processing": "cached GEO SOFT normalized values; highest all-sample-mean probe",
            "expression_scale": "supplied log2-like normalized microarray intensity",
            "contrast": "DKD_vs_Control",
            "input_path": paths["gse30528_soft"].resolve(),
            "probe_metadata": selected_matrices["GSE30528"][[
                "selected_probe_id",
                "selected_probe_mean",
            ]],
        },
        {
            "cohort": "GSE30529",
            "expression": selected_matrices["GSE30529"].drop(
                columns=["selected_probe_id", "selected_probe_mean"]
            ),
            "case_ids": validation_members["GSE30529"]["case"],
            "control_ids": validation_members["GSE30529"]["control"],
            "platform": "GPL571; Affymetrix Human Genome U133A 2.0",
            "tissue": "tubuli",
            "processing": "cached GEO SOFT normalized values; highest all-sample-mean probe",
            "expression_scale": "supplied log2-like normalized microarray intensity",
            "contrast": "DKD_vs_Control",
            "input_path": paths["gse30529_soft"].resolve(),
            "probe_metadata": selected_matrices["GSE30529"][[
                "selected_probe_id",
                "selected_probe_mean",
            ]],
        },
        {
            "cohort": "GSE142025_C1",
            "expression": gse142025_expression,
            "case_ids": c1_members["case"],
            "control_ids": c1_members["control"],
            "platform": "GPL20301; Illumina HiSeq 4000 RNA-seq",
            "tissue": "whole kidney biopsy",
            "processing": "project frozen clean matrix; source log2 transform and quantile normalization",
            "expression_scale": "log2-transformed quantile-normalized expression",
            "contrast": "Advanced_DKD_vs_Control_only",
            "input_path": paths["gse142025_expression"].resolve(),
            "probe_metadata": None,
        },
    ]
    per_cohort = build_per_cohort_effect_sizes(candidates, cohort_definitions)
    expected_cohorts = [cohort["cohort"] for cohort in cohort_definitions]
    meta = build_meta_results(candidates, per_cohort, expected_cohorts)
    stage12 = pd.read_csv(paths["stage12_results"])
    stage12_comparison = compare_stage12(stage12, meta)
    stage12_comparison["stage12_result_path"] = str(paths["stage12_results"].resolve())
    stage12_comparison["m8_meta_result_path"] = str(
        (TABLE_DIR / "meta_analysis_unified_scale_61genes.csv").resolve()
    )
    direction_test, correlation_table, combined_validation = (
        build_validation_direction_and_correlation(
            candidates, validation_results, per_cohort
        )
    )
    for accession, result in validation_results.items():
        write_table(result, f"{accession}_61gene_validation.csv", logger)
    write_table(direction_test, "direction_concordance_test.csv", logger)
    write_table(correlation_table, "effect_size_correlation.csv", logger)
    write_table(per_cohort, "per_cohort_effect_sizes.csv", logger)
    write_table(meta, "meta_analysis_unified_scale_61genes.csv", logger)
    write_table(stage12_comparison, "stage12_vs_m8_comparison.csv", logger)

    _, gse111154_mapping, gse111154_headers = parse_geo_soft(paths["gse111154_soft"])
    gse111154_metadata = pd.read_csv(paths["gse111154_metadata"])
    unique_signature_genes = sorted({gene for genes in SIGNATURES.values() for gene in genes})
    gse111154_stage14 = pd.read_csv(paths["stage14_results"])
    if len(gse111154_mapping) != 0 or "Gene Symbol" in gse111154_headers:
        raise AssertionError("GSE111154 audit expected the Stage 14 Gene Symbol mapping failure")

    reports = build_reports(
        paths=paths,
        candidates=candidates,
        metadata142025=gse142025_metadata,
        c1=c1,
        c2=c2,
        c2_check=c2_check,
        validation_results=validation_results,
        validation_members=validation_members,
        direction_test=direction_test,
        correlation_table=correlation_table,
        per_cohort=per_cohort,
        meta=meta,
        stage12_comparison=stage12_comparison,
        gse111154_metadata=gse111154_metadata,
        gse111154_headers=gse111154_headers,
        gse111154_stage14=gse111154_stage14,
        unique_signature_genes=unique_signature_genes,
    )
    for filename, content in reports.items():
        write_report(content, filename, logger)

    settings = {
        "stage": "Stage 21-M8",
        "execution_utc": started.isoformat(),
        "random_seed": RANDOM_SEED,
        "candidate_denominator": 61,
        "within_cohort_test": "Welch two-sample t-test; two-sided; equal_var=False",
        "multiple_testing": "Benjamini-Hochberg; explicit 61-gene family",
        "meta_effect": "Hedges g from statsmodels.effectsize_smd",
        "meta_model": "DerSimonian-Laird random effects",
        "direction_test": "one-sided exact binomial, p0=0.5",
        "correlation": "Spearman, two-sided",
        "gse142025_meta_contrast": "C1 Advanced_DKD vs Control only",
        "not_run": ["WGCNA", "single-nucleus", "CellPhoneDB", "clinical association"],
    }
    settings_path = TABLE_DIR / "analysis_settings.json"
    settings_path.write_text(json.dumps(settings, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    logger.info("OUTPUT_SETTINGS path=%s", settings_path)
    finished = datetime.now(timezone.utc)
    logger.info("M8_COMPLETE utc=%s elapsed_seconds=%.3f", finished.isoformat(), (finished-started).total_seconds())
    reversals = int(stage12_comparison["direction_reversal"].sum())
    inconsistency_count = 6 + int(
        stage12_comparison["conclusion_change"].ne("UNCHANGED").sum()
    )
    print(f"C1 | adjusted P < 0.05 | {int(c1['adjusted_p_value'].lt(0.05).sum())}")
    print(f"C2 | adjusted P < 0.05 | {int(c2['adjusted_p_value'].lt(0.05).sum())}")
    print(f"META | direction reversals | {reversals}")
    for row in direction_test.itertuples(index=False):
        print(
            f"{row.validation_scope} | direction concordance | "
            f"{row.direction_concordant_genes}/{row.mapped_genes} | P={row.p_value:.15g}"
        )
    print(f"MANUSCRIPT_INCONSISTENCIES | {inconsistency_count}")
    print(
        "OUTPUT_PATHS | "
        f"tables={TABLE_DIR.resolve()} | docs={DOC_DIR.resolve()} | "
        f"scripts={SCRIPT_PATH.parent.resolve()} | logs={LOG_DIR.resolve()}"
    )


def build_reports(**data: Any) -> dict[str, str]:
    paths = data["paths"]
    candidates = data["candidates"]
    metadata = data["metadata142025"]
    c1 = data["c1"]
    c2 = data["c2"]
    c2_check = data["c2_check"]
    validations = data["validation_results"]
    members = data["validation_members"]
    direction_test = data["direction_test"]
    correlation_table = data["correlation_table"]
    per_cohort = data["per_cohort"]
    meta = data["meta"]
    stage12_comparison = data["stage12_comparison"]
    metadata111154 = data["gse111154_metadata"]
    headers111154 = data["gse111154_headers"]
    stage14_result = data["gse111154_stage14"]
    unique_signature_genes = data["unique_signature_genes"]

    group_counts = (
        metadata.groupby("disease_group_clean", dropna=False)
        .size()
        .rename("n")
        .reset_index()
    )
    group_ids = (
        metadata.groupby("disease_group_clean")["sample_id"]
        .apply(lambda values: ";".join(values.astype(str)))
        .rename("sample_ids")
        .reset_index()
    )
    group_table = group_counts.merge(group_ids, on="disease_group_clean")
    c1_sig = int(c1["adjusted_p_value"].lt(0.05).sum())
    c2_sig = int(c2["adjusted_p_value"].lt(0.05).sum())
    c1_c2 = c1[["gene_symbol", "log2FC", "adjusted_p_value", "direction"]].merge(
        c2[["gene_symbol", "log2FC", "adjusted_p_value", "direction"]],
        on="gene_symbol",
        suffixes=("_C1", "_C2"),
    )
    c1_c2["direction_changed"] = c1_c2["direction_C1"].ne(c1_c2["direction_C2"])
    c1_c2["fdr_status_changed"] = c1_c2["adjusted_p_value_C1"].lt(0.05).ne(
        c1_c2["adjusted_p_value_C2"].lt(0.05)
    )
    history_available = c2_check["historical_comparator_available"].eq(True)
    raw_reproduced = int(
        (
            history_available
            & c2_check["effect_matches_1e_12"].eq(True)
            & c2_check["raw_p_matches_1e_12"].eq(True)
        ).sum()
    )
    stage10_group_line = find_line(paths["stage10_script"], "gg=pd.Series")
    stage10_group_code = paths["stage10_script"].read_text(encoding="utf-8").splitlines()[
        stage10_group_line - 1
    ].strip()
    contrast_report = f"""# Contrast definition correction

## Sample composition

Source phenotype field: `disease_group_clean` in `{paths['gse142025_metadata'].resolve()}`; this field was derived from GEO `characteristics_ch1`. Total samples: **{len(metadata)}**.

{markdown_table(group_table, ['disease_group_clean', 'n', 'sample_ids'])}

## Historical Stage 10 grouping

`{paths['stage10_script'].resolve()}:{stage10_group_line}` assigns `Advanced_DKD` to the DKD side and every other sample to the control side. The historical file label `GSE142025_advanced_vs_control` therefore represents 21 Advanced_DKD versus 6 Early_DKD plus 9 Control samples, not Advanced_DKD versus healthy Control only.

Exact grouping line:

```python
{stage10_group_code}
```

## Frozen M8 contrasts

- C1: Advanced_DKD versus Control only; n=21 versus n=9. The healthy-control side is not below 5, so the requested `<5` warning is not triggered.
- C2: Advanced_DKD versus Early_DKD plus Control; n=21 versus n=15. This reproduces Stage 10's actual grouping.
- Test: two-sided Welch independent-samples t-test, `equal_var=False`, `nan_policy='omit'`.
- Multiplicity: explicit Benjamini-Hochberg correction across the frozen 61 candidates per contrast.
- Effect: case minus control mean on the supplied log2-scale matrix, reported as `log2FC`.

## Results and reproduction

- C1 adjusted P<0.05: **{c1_sig}/61**.
- C2 adjusted P<0.05: **{c2_sig}/61**.
- C1/C2 direction changes: **{int(c1_c2['direction_changed'].sum())}**.
- C1/C2 FDR<0.05 status changes: **{int(c1_c2['fdr_status_changed'].sum())}**.
- Stage 10 historical comparators available: **{int(history_available.sum())}/61**; the remaining {61-int(history_available.sum())} candidates were never calculated by Stage 10.
- Historical effect and raw P reproduced within 1e-12 for **{raw_reproduced}/{int(history_available.sum())}** available genes.
- Adjusted P values are not expected to reproduce because Stage 10 corrected 8 tests and M8 corrects the prespecified 61-test family. The reproduction table reports this separately instead of treating it as raw-code disagreement.

Complete sample IDs and statistics are in `contrast_C1_results.csv`, `contrast_C2_results.csv`, and `C2_vs_stage10_reproduction_check.csv`.
"""

    scale_rows = (
        per_cohort[
            ["cohort", "platform", "tissue", "processing", "expression_scale", "contrast"]
        ]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    stage12_weight_line = find_line(paths["stage12_script"], "ext['se']")
    stage12_meta_line = find_line(paths["stage12_script"], "def meta")
    stage12_weight_code = paths["stage12_script"].read_text(encoding="utf-8").splitlines()[
        stage12_weight_line - 1
    ].strip()
    stage12_meta_code = paths["stage12_script"].read_text(encoding="utf-8").splitlines()[
        stage12_meta_line - 1
    ].strip()
    reversals = int(stage12_comparison["direction_reversal"].sum())
    losses = int(stage12_comparison["conclusion_change"].str.contains("SIGNIFICANCE_LOSS").sum())
    gains = int(stage12_comparison["conclusion_change"].str.contains("SIGNIFICANCE_GAIN").sum())
    effect_report = f"""# Effect-scale unification

## Existing Stage 12 diagnosis

{markdown_table(scale_rows, ['cohort', 'platform', 'tissue', 'processing', 'expression_scale', 'contrast'])}

The four raw mean-difference scales are not directly interchangeable: cohorts span microarray and RNA-seq-derived matrices, separate normalization procedures, and different numerical dispersions. Stage 12 read the raw mean differences, reconstructed external SE values from effect/P at `{paths['stage12_script'].resolve()}:{stage12_weight_line}`, and applied inverse-variance DerSimonian-Laird calculations beginning at line {stage12_meta_line}. Inverse-variance weighting does not make unlike effect units commensurable.

Exact Stage 12 code excerpts:

```python
{stage12_meta_code}
{stage12_weight_code}
```

## Frozen M8 method

M8 selected Hedges' g before inspecting results because it standardizes each cohort's case-control mean difference by its within-cohort pooled dispersion and applies a small-sample correction. Sampling variance comes from `statsmodels.stats.meta_analysis.effectsize_smd`. Per-gene synthesis uses DerSimonian-Laird tau-squared, random inverse-variance weights, normal 95% CI and pooled two-sided P, Cochran Q and chi-square Q P, I2, tau2, and BH correction across 61 genes. GSE142025 contributes C1 only; C2 is excluded because it mixes progression and disease-presence contrasts.

## New versus old eight-gene conclusions

- Direction reversals: **{reversals}**.
- FDR<0.05 losses: **{losses}**.
- FDR<0.05 gains: **{gains}**.

{markdown_table(stage12_comparison, ['gene_symbol', 'stage12_raw_scale_effect', 'stage12_fdr', 'm8_hedges_g', 'm8_fdr', 'direction_reversal', 'conclusion_change'])}

The old and new pooled magnitudes use different units and must not be compared numerically. Direction, interval, and inferential-status changes are reported without suppression in `stage12_vs_m8_comparison.csv`.
"""

    mapping_rows = []
    for accession, result in validations.items():
        mapped = int(result["mapping_status"].eq("MAPPED_TESTED").sum())
        unmapped = result.loc[result["mapping_status"].eq("NOT_MAPPED"), "gene_symbol"].tolist()
        mapping_rows.append(
            {
                "dataset": accession,
                "n_case": len(members[accession]["case"]),
                "n_control": len(members[accession]["control"]),
                "mapped": mapped,
                "unmapped": 61 - mapped,
                "unmapped_genes": ";".join(unmapped) if unmapped else "NONE",
            }
        )
    mapping_frame = pd.DataFrame(mapping_rows)
    correlation_summary = correlation_table.loc[
        correlation_table["record_type"].eq("SUMMARY")
    ][["validation_dataset", "spearman_rho", "correlation_p_value", "n_mapped"]]
    stage10_matches = []
    main_gene_frames: list[pd.DataFrame] = []
    for accession, result in validations.items():
        available = result["stage10_comparator_available"].eq(True)
        exact = (
            available
            & result["stage10_effect_matches_1e_12"].eq(True)
            & result["stage10_raw_p_matches_1e_12"].eq(True)
        )
        stage10_matches.append(f"{accession}: {int(exact.sum())}/{int(available.sum())}")
        main_gene_frame = result.loc[result["gene_symbol"].isin(MAIN_GENES), [
            "dataset",
            "gene_symbol",
            "log2FC",
            "p_value",
            "adjusted_p_value",
            "direction",
            "direction_matches_discovery",
            "stage10_effect_matches_1e_12",
            "stage10_raw_p_matches_1e_12",
        ]].copy()
        main_gene_frame["gene_order"] = main_gene_frame["gene_symbol"].map(
            {gene: index for index, gene in enumerate(MAIN_GENES)}
        )
        main_gene_frames.append(main_gene_frame.sort_values("gene_order"))
    main_gene_results = pd.concat(main_gene_frames, ignore_index=True)
    validation_report = f"""# Independent 61-gene validation

## Cohorts and mapping

- GSE30528: GPL571 glomeruli; DKD 9, Control 13.
- GSE30529: GPL571 tubuli; DKD 10, Control 12.
- Mapping follows the frozen Stage 10 rule: platform `Gene Symbol`, first ` /// ` symbol, uppercase, then the probe with the highest across-all-samples mean. Probe selection does not use group labels.

{markdown_table(mapping_frame, ['dataset', 'n_case', 'n_control', 'mapped', 'unmapped', 'unmapped_genes'])}

Every candidate remains in each output even when unmapped. Welch P values are available for 60 mapped genes; BH correction is applied to those 60 finite P values while the frozen 61-row denominator retains C1QC explicitly as `NOT_MAPPED` with missing statistics.

## Direction concordance

Discovery direction is the frozen `advanced_logFC_limma` sign; its sign agrees with GSE96804 for 61/61 candidates before validation. Each validation cohort and the per-gene combined validation Hedges' g are tested with a one-sided exact binomial test against p0=0.5. Genes are biologically correlated, so the test's independence assumption is an explicit limitation.

{markdown_table(direction_test, ['validation_scope', 'mapped_genes', 'direction_concordant_genes', 'direction_discordant_genes', 'concordance_fraction', 'p_value'])}

## Effect correlation

Spearman correlation was fixed in advance because raw effects have different cross-platform scales.

{markdown_table(correlation_summary, ['validation_dataset', 'spearman_rho', 'correlation_p_value', 'n_mapped'])}

Point-level scatter data are in `effect_size_correlation.csv`; no figure was drawn.

## Stage 10 eight-gene reproduction

Effect and raw P matches within 1e-12: **{'; '.join(stage10_matches)}**. Adjusted P is separately labeled because Stage 10 used an 8-gene family and M8 uses the 60 mapped tests while retaining all 61 candidates as rows.

{markdown_table(main_gene_results, ['dataset', 'gene_symbol', 'log2FC', 'p_value', 'adjusted_p_value', 'direction', 'direction_matches_discovery', 'stage10_effect_matches_1e_12', 'stage10_raw_p_matches_1e_12'])}
"""

    stage14_parser_line = find_line(paths["stage14_script"], "z.get('Gene Symbol'")
    stage14_score_line = find_line(paths["stage14_script"], "q=scores")
    edn_count = int(
        metadata111154["characteristics"].str.contains(
            "early diabetic", case=False, na=False
        ).sum()
    )
    control_count = len(metadata111154) - edn_count
    all_na = int(stage14_result[["early_minus_control", "P.Value"]].isna().all(axis=1).sum())
    negative_report = f"""# GSE111154 non-estimability report

## What Stage 14 actually did

Stage 14 read the cached GPL17586 family SOFT, assigned 4 early diabetic nephropathy and 4 Control samples from the `characteristics` field, attempted to collapse probes to genes, z-scored each mapped gene across samples, averaged predefined signature members, and attempted Welch early-DKD-versus-Control tests. It did not have individual pathology-stage labels.

## Mapping failure

At `{paths['stage14_script'].resolve()}:{stage14_parser_line}`, Stage 14 requested a platform column named `Gene Symbol`. The actual GPL17586 header is:

`{' | '.join(headers111154)}`

It contains `gene_assignment` but not `Gene Symbol`; therefore the coded mapping dictionary contained zero probes/genes. At line {stage14_score_line}, all 7 unique predefined signature genes (`{'; '.join(unique_signature_genes)}`) were reindexed against an empty gene matrix. Mapping rate under the Stage 14 implementation was **0/7 (0%)**, and all {all_na}/4 signature result rows had missing effect and P values.

This is an implementation-limited non-estimability result, not evidence that the biological signature is absent and not a platform-level proof that the genes cannot be mapped. Correctly parsing `gene_assignment` would be a new analysis and was prohibited in M8 Task 4, so it was not attempted.

## Verifiable boundaries and possible contributors

- Primary demonstrated cause: platform-column mismatch in the Stage 14 parser.
- Sample size: 4 early-DKD and 4 Control samples would yield limited power even after a corrected mapping.
- Platform/context: GPL17586 and kidney cortical tissue differ from the discovery inputs; the consequence was not estimated here.
- Signature granularity: the four signatures contain only 2-4 members and 7 unique genes, making coverage sensitive to mapping loss.
- Pathology: GEO sample metadata provides early-DKD/Control groups, not individual RPS stage, so no individual pathology-stage association can be made.

No trend, near-significance, or positive implication is assigned to this empty result.
"""

    changed = stage12_comparison.loc[
        stage12_comparison["conclusion_change"].ne("UNCHANGED")
    ]
    inconsistencies = [
        "Stage 10 label GSE142025_advanced_vs_control represented Advanced_DKD versus Early_DKD plus Control, not healthy Control only.",
        "Stage 12 pooled raw mean differences from unlike cohort scales; M8 replaces this with Hedges' g and the Stage 12 pooled magnitude must not be treated as a common-unit effect.",
        "GSE142025's contribution to the disease-versus-control meta-analysis is now C1 (Advanced_DKD versus Control only); the historical C2 mixed progression and disease-presence contrasts.",
        "Stage 10 external_validation_effect_size_forestplot.png is a grouped bar chart without confidence intervals, not a forest plot.",
        "Stage 10 gene-set values are averages of member-gene effects with no gene-set P value, not sample-level standardized expression scores.",
        "GSE111154 Stage 14 non-estimability was caused by a platform-column parser mismatch; it is not an interpretable biological negative validation.",
    ]
    for row in changed.itertuples(index=False):
        inconsistencies.append(
            f"{row.gene_symbol}: Stage 12 versus M8 conclusion status is {row.conclusion_change}."
        )
    inconsistency_lines = "\n".join(
        f"{index}. {item}" for index, item in enumerate(inconsistencies, 1)
    )
    direction_lines = "\n".join(
        f"- {row.validation_scope}: {row.direction_concordant_genes}/{row.mapped_genes}, exact-binomial P={row.p_value:.6g}"
        for row in direction_test.itertuples(index=False)
    )
    correlation_lines = "\n".join(
        f"- {row.validation_dataset}: Spearman rho={row.spearman_rho:.6g}, P={row.correlation_p_value:.6g}, n={row.n_mapped}"
        for row in correlation_summary.itertuples(index=False)
    )
    summary = f"""# Stage 21-M8 summary

## 1. Contrast correction

C1 (21 Advanced_DKD versus 9 Control) produced **{c1_sig}** genes with adjusted P<0.05. C2 (21 Advanced_DKD versus 15 Early_DKD+Control) produced **{c2_sig}**. They differed in direction for **{int(c1_c2['direction_changed'].sum())}** genes and in FDR<0.05 status for **{int(c1_c2['fdr_status_changed'].sum())}** genes. The historical C2 raw effect/P values reproduced for {raw_reproduced}/8 available Stage 10 genes; BH values differ where the correction family changed from 8 to 61.

## 2. Unified-scale meta-analysis

The frozen replacement uses Hedges' g and DerSimonian-Laird random effects across GSE96804, GSE30528, GSE30529, and GSE142025 C1. Among the original eight Stage 12 genes: direction reversals **{reversals}**, FDR<0.05 losses **{losses}**, gains **{gains}**. Every change is itemized in `stage12_vs_m8_comparison.csv` and below.

## 3. Complete 61-gene validation

{direction_lines}

{correlation_lines}

The 61-gene denominator was frozen before analysis; no candidate was selected or removed based on validation results.

## 4. GSE111154

Stage 14's `Gene Symbol` parser did not match GPL17586's `gene_assignment` header. Zero of 7 unique signature genes entered scoring; all four signature result rows were non-estimable. This is an implementation failure with 4 early-DKD and 4 Control samples, not a biological negative trend and not an individual pathology-stage analysis.

## 5. Inconsistencies with existing manuscript-era conclusions or labels

Count: **{len(inconsistencies)}**.

{inconsistency_lines}

## 6. Analyses not performed

- WGCNA: frozen prior results retained; Rscript unavailable and M8 prohibited rerun.
- Single-nucleus localization: frozen prior results retained; Scanpy/AnnData unavailable and M8 prohibited rerun.
- CellPhoneDB: frozen prior results retained; M8 prohibited rerun.
- GSE142025 clinical association: not attempted because eGFR, proteinuria, albuminuria, age, and sex are 0/36 complete.
- GSE111154 corrected `gene_assignment` analysis: not attempted because Task 4 allowed organization of the existing result only.
"""
    return {
        "contrast_definition_correction.md": contrast_report,
        "effect_scale_unification.md": effect_report,
        "independent_validation_report.md": validation_report,
        "GSE111154_negative_result_report.md": negative_report,
        "M8_SUMMARY.md": summary,
    }


if __name__ == "__main__":
    main()
