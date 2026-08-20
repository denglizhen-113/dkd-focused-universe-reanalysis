from __future__ import annotations

import hashlib
import importlib.util
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import t, ttest_ind


ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = Path(__file__).resolve()
TABLE_DIR = ROOT / "tables" / "stage21_m18_remediation"
DOC_DIR = ROOT / "docs" / "stage21_m18_remediation"
LOG_DIR = ROOT / "logs" / "stage21_m18_remediation"
M8_PATH = ROOT / "scripts" / "stage21_m8_validation" / "run_m8_analysis.py"
RANDOM_SEED = 20260815
N_PERMUTATIONS = 10_000
PRIMARY_COHORTS = ["GSE142025_C1", "GSE96804", "GSE30528"]
DISCOVERY_COHORTS = ["GSE142025_C1", "GSE96804"]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bh_family(p_values: pd.Series, family_size: int) -> pd.Series:
    """BH adjustment while retaining missing/non-estimable members as P=1."""

    if len(p_values) != family_size:
        raise ValueError(f"Expected family size {family_size}, found {len(p_values)}")
    p = pd.to_numeric(p_values, errors="coerce").fillna(1.0).clip(0.0, 1.0).to_numpy()
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * family_size / np.arange(1, family_size + 1))[::-1])[::-1]
    result = np.empty(family_size, dtype=float)
    result[order] = np.minimum(adjusted, 1.0)
    return pd.Series(result, index=p_values.index)


def restricted_loglik(tau2: float, y: np.ndarray, v: np.ndarray) -> float:
    total_variance = v + tau2
    weights = 1.0 / total_variance
    mean = float(np.sum(weights * y) / np.sum(weights))
    residual = float(np.sum(weights * np.square(y - mean)))
    # Additive constants independent of tau2 are omitted.
    return -0.5 * (
        float(np.sum(np.log(total_variance)))
        + math.log(float(np.sum(weights)))
        + residual
    )


def reml_hartung_knapp(y: np.ndarray, v: np.ndarray) -> dict[str, float]:
    """REML random-effects mean with unmodified and modified HKSJ uncertainty."""

    y = np.asarray(y, dtype=float)
    v = np.asarray(v, dtype=float)
    valid = np.isfinite(y) & np.isfinite(v) & (v > 0)
    y = y[valid]
    v = v[valid]
    k = len(y)
    empty = {
        "k": float(k),
        "pooled_effect": math.nan,
        "tau2_reml": math.nan,
        "hk_scale": math.nan,
        "pooled_se_hk": math.nan,
        "pooled_se_modified_hk": math.nan,
        "ci_95_low_hk": math.nan,
        "ci_95_high_hk": math.nan,
        "ci_95_low_modified_hk": math.nan,
        "ci_95_high_modified_hk": math.nan,
        "p_value_hk": math.nan,
        "p_value_modified_hk": math.nan,
        "prediction_interval_95_low": math.nan,
        "prediction_interval_95_high": math.nan,
        "prediction_interval_95_low_modified_hk": math.nan,
        "prediction_interval_95_high_modified_hk": math.nan,
        "cochran_q": math.nan,
        "i2_percent": math.nan,
    }
    if k < 2:
        return empty

    upper = max(10.0, float(np.var(y, ddof=1) * 100.0 + np.max(v)))
    optimized = minimize_scalar(
        lambda tau2: -restricted_loglik(float(tau2), y, v),
        bounds=(0.0, upper),
        method="bounded",
        options={"xatol": 1e-12},
    )
    tau2 = float(max(0.0, optimized.x))
    # A boundary solution should be represented as exactly zero.
    if tau2 < 1e-10 or restricted_loglik(0.0, y, v) >= restricted_loglik(tau2, y, v) - 1e-10:
        tau2 = 0.0

    weights = 1.0 / (v + tau2)
    weight_sum = float(np.sum(weights))
    pooled = float(np.sum(weights * y) / weight_sum)
    df = k - 1
    hk_scale = float(np.sum(weights * np.square(y - pooled)) / df)
    se_hk = math.sqrt(max(hk_scale, 0.0) / weight_sum)
    se_modified = math.sqrt(max(hk_scale, 1.0) / weight_sum)
    critical = float(t.ppf(0.975, df))
    p_value = (
        float(2.0 * t.sf(abs(pooled / se_hk), df))
        if se_hk > 0
        else (0.0 if pooled != 0 else 1.0)
    )
    p_value_modified = (
        float(2.0 * t.sf(abs(pooled / se_modified), df))
        if se_modified > 0
        else (0.0 if pooled != 0 else 1.0)
    )

    fixed_weights = 1.0 / v
    fixed_mean = float(np.sum(fixed_weights * y) / np.sum(fixed_weights))
    q = float(np.sum(fixed_weights * np.square(y - fixed_mean)))
    i2 = max(0.0, (q - df) / q * 100.0) if q > 0 else 0.0

    pi_low = math.nan
    pi_high = math.nan
    pi_low_modified = math.nan
    pi_high_modified = math.nan
    if k >= 3:
        pi_critical = float(t.ppf(0.975, k - 2))
        pi_half_width = pi_critical * math.sqrt(tau2 + se_hk**2)
        pi_half_width_modified = pi_critical * math.sqrt(tau2 + se_modified**2)
        pi_low = pooled - pi_half_width
        pi_high = pooled + pi_half_width
        pi_low_modified = pooled - pi_half_width_modified
        pi_high_modified = pooled + pi_half_width_modified

    return {
        "k": float(k),
        "pooled_effect": pooled,
        "tau2_reml": tau2,
        "hk_scale": hk_scale,
        "pooled_se_hk": se_hk,
        "pooled_se_modified_hk": se_modified,
        "ci_95_low_hk": pooled - critical * se_hk,
        "ci_95_high_hk": pooled + critical * se_hk,
        "ci_95_low_modified_hk": pooled - critical * se_modified,
        "ci_95_high_modified_hk": pooled + critical * se_modified,
        "p_value_hk": p_value,
        "p_value_modified_hk": p_value_modified,
        "prediction_interval_95_low": pi_low,
        "prediction_interval_95_high": pi_high,
        "prediction_interval_95_low_modified_hk": pi_low_modified,
        "prediction_interval_95_high_modified_hk": pi_high_modified,
        "cochran_q": q,
        "i2_percent": i2,
    }


def read_focused_universe() -> tuple[pd.DataFrame, dict[str, list[str]]]:
    phase3_path = ROOT / "tables" / "phase3" / "core_candidate_gene_table.csv"
    phase3 = pd.read_csv(phase3_path)
    phase3["gene_symbol"] = phase3["gene_symbol"].astype(str).str.upper()
    if len(phase3) != 145 or phase3["gene_symbol"].nunique() != 145:
        raise AssertionError("The focused universe must contain 145 unique genes")
    candidates = pd.DataFrame(
        {
            "gene_symbol": phase3["gene_symbol"],
            "theme": phase3["focused_pathways"],
            "candidate_source_row": np.arange(2, len(phase3) + 2),
            "discovery_effect": pd.to_numeric(
                phase3["GSE142025_advanced_vs_control_effect_size_log_expression"],
                errors="coerce",
            ),
        }
    )
    candidates["discovery_direction"] = np.where(
        candidates["discovery_effect"] > 0,
        "UP",
        np.where(candidates["discovery_effect"] < 0, "DOWN", "ZERO"),
    )

    pathways: dict[str, list[str]] = {}
    gmt_path = ROOT / "config" / "phase3_focus_gene_sets.gmt"
    for line in gmt_path.read_text(encoding="utf-8").splitlines():
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 3:
            pathways[fields[0]] = [gene.upper() for gene in fields[2:]]
    if len(pathways) != 9 or len({gene for genes in pathways.values() for gene in genes}) != 145:
        raise AssertionError("GMT must contain 9 pathways and 145 unique genes")
    return candidates, pathways


def prepare_cohorts(m8, candidates: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    expression_142 = m8.load_clean_expression(
        ROOT / "data_processed" / "GSE142025" / "expression_matrix_clean.csv"
    )
    metadata_142 = pd.read_csv(
        ROOT / "data_processed" / "GSE142025" / "sample_annotation_clean.csv"
    )
    c1 = m8.build_gse142025_contrasts(metadata_142)["C1_advanced_vs_control"]

    expression_96804 = m8.load_clean_expression(
        ROOT / "data_processed" / "GSE96804" / "expression_matrix_clean.csv"
    )
    metadata_96804 = pd.read_csv(
        ROOT / "data_processed" / "GSE96804" / "sample_annotation_clean.csv"
    )
    case_96804 = metadata_96804.loc[
        metadata_96804["disease_group_clean"].eq("DKD"), "sample_id"
    ].astype(str).tolist()
    control_96804 = metadata_96804.loc[
        metadata_96804["disease_group_clean"].eq("Control"), "sample_id"
    ].astype(str).tolist()

    selected: dict[str, pd.DataFrame] = {}
    members: dict[str, dict[str, list[str]]] = {}
    metadata_geo: dict[str, pd.DataFrame] = {}
    for accession in ("GSE30528", "GSE30529"):
        soft_path = ROOT / "data_raw" / "geo_soft" / accession / f"{accession}_family.soft.gz"
        probe_matrix, probe_to_gene, _ = m8.parse_geo_soft(soft_path)
        selected[accession] = m8.select_highest_mean_probes(probe_matrix, probe_to_gene)
        metadata = pd.read_csv(ROOT / "data_raw" / "geo_soft" / accession / "samples.csv")
        metadata_geo[accession] = metadata
        is_dkd = metadata["characteristics"].str.contains(
            "diabetic kidney disease", case=False, na=False
        )
        members[accession] = {
            "case": metadata.loc[is_dkd, "gsm_accession"].astype(str).tolist(),
            "control": metadata.loc[~is_dkd, "gsm_accession"].astype(str).tolist(),
        }

    definitions = [
        {
            "cohort": "GSE142025_C1",
            "expression": expression_142,
            "case_ids": c1["case"],
            "control_ids": c1["control"],
            "platform": "GPL20301; Illumina HiSeq 4000 RNA-seq",
            "tissue": "whole kidney biopsy",
            "processing": "project frozen clean matrix",
            "expression_scale": "log2-transformed quantile-normalized expression",
            "contrast": "Advanced_DKD_vs_Control_only",
            "input_path": ROOT / "data_processed" / "GSE142025" / "expression_matrix_clean.csv",
            "probe_metadata": None,
        },
        {
            "cohort": "GSE96804",
            "expression": expression_96804,
            "case_ids": case_96804,
            "control_ids": control_96804,
            "platform": "GPL17586; Affymetrix Human Transcriptome Array 2.0",
            "tissue": "glomeruli",
            "processing": "project frozen clean matrix",
            "expression_scale": "log2-scale normalized microarray expression",
            "contrast": "DKD_vs_Control",
            "input_path": ROOT / "data_processed" / "GSE96804" / "expression_matrix_clean.csv",
            "probe_metadata": None,
        },
    ]
    for accession, tissue in (("GSE30528", "glomeruli"), ("GSE30529", "tubulointerstitium")):
        definitions.append(
            {
                "cohort": accession,
                "expression": selected[accession].drop(
                    columns=["selected_probe_id", "selected_probe_mean"]
                ),
                "case_ids": members[accession]["case"],
                "control_ids": members[accession]["control"],
                "platform": "GPL571; Affymetrix Human Genome U133A 2.0",
                "tissue": tissue,
                "processing": "cached GEO SOFT normalized values; highest all-sample-mean probe",
                "expression_scale": "supplied log2-like normalized microarray intensity",
                "contrast": "DKD_vs_Control",
                "input_path": ROOT / "data_raw" / "geo_soft" / accession / f"{accession}_family.soft.gz",
                "probe_metadata": selected[accession][["selected_probe_id", "selected_probe_mean"]],
            }
        )
    per_cohort = m8.build_per_cohort_effect_sizes(candidates, definitions)
    context = {
        "selected": selected,
        "members": members,
        "metadata_geo": metadata_geo,
        "definitions": definitions,
    }
    return per_cohort, context


def discovery_direction_table(per_cohort: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for gene, sub in per_cohort.loc[per_cohort["cohort"].isin(DISCOVERY_COHORTS)].groupby(
        "gene_symbol", sort=False
    ):
        valid = sub.loc[
            np.isfinite(sub["hedges_g"])
            & np.isfinite(sub["variance"])
            & sub["variance"].gt(0)
        ]
        if len(valid) == 2:
            weights = 1.0 / valid["variance"].to_numpy(float)
            effect = float(np.sum(weights * valid["hedges_g"].to_numpy(float)) / np.sum(weights))
        else:
            effect = math.nan
        rows.append(
            {
                "gene_symbol": gene,
                "discovery_effect_fixed": effect,
                "discovery_direction": (
                    "UP" if effect > 0 else "DOWN" if effect < 0 else "NOT_ESTIMABLE"
                ),
                "n_discovery_studies": len(valid),
            }
        )
    return pd.DataFrame(rows)


def run_primary_meta(
    candidates: pd.DataFrame, per_cohort: pd.DataFrame, directions: pd.DataFrame
) -> pd.DataFrame:
    direction_index = directions.set_index("gene_symbol")
    rows: list[dict[str, Any]] = []
    for candidate in candidates.itertuples(index=False):
        sub = per_cohort.loc[
            per_cohort["gene_symbol"].eq(candidate.gene_symbol)
            & per_cohort["cohort"].isin(PRIMARY_COHORTS)
        ]
        valid = sub.loc[
            np.isfinite(sub["hedges_g"])
            & np.isfinite(sub["variance"])
            & sub["variance"].gt(0)
        ]
        meta = reml_hartung_knapp(
            valid["hedges_g"].to_numpy(float), valid["variance"].to_numpy(float)
        )
        complete = len(valid) == 3
        discovery_direction = direction_index.loc[candidate.gene_symbol, "discovery_direction"]
        pooled_direction = (
            "UP"
            if meta["pooled_effect"] > 0
            else "DOWN"
            if meta["pooled_effect"] < 0
            else "NOT_ESTIMABLE"
        )
        rows.append(
            {
                "gene_symbol": candidate.gene_symbol,
                "focused_pathways": candidate.theme,
                "analysis_family": "analysis_defined_145_gene_focused_universe",
                "primary_estimand": "three-study synthesis with GSE30122 represented by glomerular GSE30528",
                "included_studies": ";".join(valid["cohort"].tolist()),
                "complete_three_study_mapping": complete,
                "discovery_direction": discovery_direction,
                "pooled_direction": pooled_direction,
                "direction_matches_discovery": pooled_direction == discovery_direction,
                **meta,
            }
        )
    frame = pd.DataFrame(rows)
    frame["p_value_for_145_family"] = np.where(
        frame["complete_three_study_mapping"], frame["p_value_modified_hk"], 1.0
    )
    frame["fdr_bh_145_family"] = bh_family(frame["p_value_for_145_family"], 145)
    frame["fdr_lt_0_05"] = frame["fdr_bh_145_family"].lt(0.05)
    frame["model"] = "REML random effects with modified Hartung-Knapp primary t inference"
    frame["inference_df"] = np.where(frame["complete_three_study_mapping"], 2, np.nan)
    frame["prediction_interval_df"] = np.where(frame["complete_three_study_mapping"], 1, np.nan)
    return frame


def run_external_gene_validation(
    candidates: pd.DataFrame,
    per_cohort: pd.DataFrame,
    directions: pd.DataFrame,
) -> pd.DataFrame:
    validation = per_cohort.loc[per_cohort["cohort"].eq("GSE30528")].copy()
    validation = validation.merge(directions, on="gene_symbol", how="left", validate="one_to_one")
    validation["validation_direction"] = np.where(
        validation["hedges_g"] > 0,
        "UP",
        np.where(validation["hedges_g"] < 0, "DOWN", "NOT_ESTIMABLE"),
    )
    # Welch P values are calculated directly over all 145 genes because the
    # frozen M8 helper intentionally asserts an exactly 61-row denominator.
    m8 = load_module("stage21_m8_for_external", M8_PATH)
    selected = prepare_external_matrix(m8, "GSE30528")
    rows: list[dict[str, Any]] = []
    for gene in candidates["gene_symbol"]:
        if gene not in selected["expression"].index:
            rows.append(
                {
                    "gene_symbol": gene,
                    "external_welch_p_value": math.nan,
                    "external_log2_mean_difference": math.nan,
                    "external_mapping_status": "NOT_MAPPED",
                }
            )
            continue
        case = selected["expression"].loc[gene, selected["case_ids"]].to_numpy(float)
        control = selected["expression"].loc[gene, selected["control_ids"]].to_numpy(float)
        test = ttest_ind(case, control, equal_var=False, nan_policy="omit")
        rows.append(
            {
                "gene_symbol": gene,
                "external_welch_p_value": float(test.pvalue),
                "external_log2_mean_difference": float(np.nanmean(case) - np.nanmean(control)),
                "external_mapping_status": "MAPPED_TESTED",
            }
        )
    direct = pd.DataFrame(rows)
    validation = validation.merge(direct, on="gene_symbol", how="left", validate="one_to_one")
    # A pair of NOT_ESTIMABLE labels must not be counted as directionally
    # concordant.  Concordance is defined only among genes that were mapped
    # and actually received an external Welch test.
    validation["direction_concordant"] = validation["external_welch_p_value"].notna() & (
        validation["validation_direction"].eq(validation["discovery_direction"])
    )
    validation["external_p_value_for_145_family"] = validation["external_welch_p_value"].fillna(1.0)
    validation["external_fdr_bh_145_family"] = bh_family(
        validation["external_p_value_for_145_family"], 145
    )
    validation["external_fdr_lt_0_05"] = validation["external_fdr_bh_145_family"].lt(0.05)
    validation["validation_scope"] = "one external source study, prespecified glomerular compartment"
    return validation


def prepare_external_matrix(m8, accession: str) -> dict[str, Any]:
    soft = ROOT / "data_raw" / "geo_soft" / accession / f"{accession}_family.soft.gz"
    probe_matrix, probe_to_gene, _ = m8.parse_geo_soft(soft)
    selected = m8.select_highest_mean_probes(probe_matrix, probe_to_gene)
    metadata = pd.read_csv(ROOT / "data_raw" / "geo_soft" / accession / "samples.csv")
    is_dkd = metadata["characteristics"].str.contains(
        "diabetic kidney disease", case=False, na=False
    )
    return {
        "expression": selected.drop(columns=["selected_probe_id", "selected_probe_mean"]),
        "probe_metadata": selected[["selected_probe_id", "selected_probe_mean"]],
        "case_ids": metadata.loc[is_dkd, "gsm_accession"].astype(str).tolist(),
        "control_ids": metadata.loc[~is_dkd, "gsm_accession"].astype(str).tolist(),
    }


def hedges_g_matrix(values: np.ndarray, case_mask: np.ndarray) -> np.ndarray:
    case = values[:, case_mask]
    control = values[:, ~case_mask]
    n_case = case.shape[1]
    n_control = control.shape[1]
    case_mean = np.mean(case, axis=1)
    control_mean = np.mean(control, axis=1)
    case_var = np.var(case, axis=1, ddof=1)
    control_var = np.var(control, axis=1, ddof=1)
    pooled_var = ((n_case - 1) * case_var + (n_control - 1) * control_var) / (
        n_case + n_control - 2
    )
    correction = 1.0 - 3.0 / (4.0 * (n_case + n_control) - 9.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        effects = correction * (case_mean - control_mean) / np.sqrt(pooled_var)
    effects[~np.isfinite(effects)] = np.nan
    return effects


def run_pathway_permutation(
    pathways: dict[str, list[str]], directions: pd.DataFrame
) -> pd.DataFrame:
    m8 = load_module("stage21_m8_for_permutation", M8_PATH)
    prepared = prepare_external_matrix(m8, "GSE30528")
    expression = prepared["expression"]
    ordered_samples = [*prepared["case_ids"], *prepared["control_ids"]]
    expression = expression.loc[:, ordered_samples]
    direction_index = directions.set_index("gene_symbol")["discovery_direction"]
    usable_genes = [
        gene
        for gene in expression.index.astype(str)
        if gene in direction_index.index and direction_index.loc[gene] in {"UP", "DOWN"}
    ]
    values = expression.loc[usable_genes].to_numpy(float)
    signs = np.array([1.0 if direction_index.loc[gene] == "UP" else -1.0 for gene in usable_genes])
    observed_mask = np.array(
        [sample in set(prepared["case_ids"]) for sample in ordered_samples], dtype=bool
    )
    observed_effects = hedges_g_matrix(values, observed_mask) * signs
    gene_to_index = {gene: idx for idx, gene in enumerate(usable_genes)}
    pathway_indexes = {
        pathway: np.array([gene_to_index[g] for g in genes if g in gene_to_index], dtype=int)
        for pathway, genes in pathways.items()
    }
    observed = {
        pathway: float(np.nanmean(observed_effects[indexes]))
        for pathway, indexes in pathway_indexes.items()
    }
    exceed = {pathway: 0 for pathway in pathways}
    rng = np.random.default_rng(RANDOM_SEED)
    for _ in range(N_PERMUTATIONS):
        permuted_case = np.zeros(len(ordered_samples), dtype=bool)
        permuted_case[rng.choice(len(ordered_samples), observed_mask.sum(), replace=False)] = True
        permuted_effects = hedges_g_matrix(values, permuted_case) * signs
        for pathway, indexes in pathway_indexes.items():
            statistic = float(np.nanmean(permuted_effects[indexes]))
            if statistic >= observed[pathway] - 1e-15:
                exceed[pathway] += 1
    rows = []
    for pathway, genes in pathways.items():
        mapped = [gene for gene in genes if gene in gene_to_index]
        rows.append(
            {
                "pathway": pathway,
                "genes_in_gmt": len(genes),
                "genes_mapped_with_discovery_direction": len(mapped),
                "observed_mean_direction_aligned_hedges_g": observed[pathway],
                "permutation_exceedances": exceed[pathway],
                "permutations": N_PERMUTATIONS,
                "one_sided_permutation_p": (exceed[pathway] + 1) / (N_PERMUTATIONS + 1),
                "permutation_unit": "GSE30528 sample disease labels; fixed 9-versus-13 allocation",
                "correlation_handling": "all pathway genes recomputed jointly for each label permutation",
            }
        )
    frame = pd.DataFrame(rows)
    frame["fdr_bh_9_pathways"] = bh_family(frame["one_sided_permutation_p"], 9)
    return frame.sort_values("fdr_bh_9_pathways", kind="stable").reset_index(drop=True)


def input_paths() -> list[Path]:
    return [
        SCRIPT_PATH,
        M8_PATH,
        ROOT / "config" / "phase3_focus_gene_sets.gmt",
        ROOT / "tables" / "phase3" / "core_candidate_gene_table.csv",
        ROOT / "data_processed" / "GSE142025" / "expression_matrix_clean.csv",
        ROOT / "data_processed" / "GSE142025" / "sample_annotation_clean.csv",
        ROOT / "data_processed" / "GSE96804" / "expression_matrix_clean.csv",
        ROOT / "data_processed" / "GSE96804" / "sample_annotation_clean.csv",
        ROOT / "data_raw" / "geo_soft" / "GSE30528" / "GSE30528_family.soft.gz",
        ROOT / "data_raw" / "geo_soft" / "GSE30528" / "samples.csv",
        ROOT / "data_raw" / "geo_soft" / "GSE30529" / "GSE30529_family.soft.gz",
        ROOT / "data_raw" / "geo_soft" / "GSE30529" / "samples.csv",
    ]


def write_reports(
    candidates: pd.DataFrame,
    primary: pd.DataFrame,
    external: pd.DataFrame,
    pathways: pd.DataFrame,
    hashes: pd.DataFrame,
) -> None:
    complete = int(primary["complete_three_study_mapping"].sum())
    significant = int(primary["fdr_lt_0_05"].sum())
    external_mapped = int(external["external_welch_p_value"].notna().sum())
    external_fdr = int(external["external_fdr_lt_0_05"].sum())
    concordant = int(external.loc[external["external_welch_p_value"].notna(), "direction_concordant"].sum())
    pathway_fdr = int(pathways["fdr_bh_9_pathways"].lt(0.05).sum())
    selected_61 = int(
        pd.read_csv(ROOT / "tables" / "phase3" / "core_candidate_gene_table.csv")[
            "prioritization_status"
        ].eq("prioritized_primary_trend_and_directional_glomerular_consistency").sum()
    )
    decision = f"""# Stage 21-M18 statistical remediation decision

## Decisions fixed before the M18 result tables were inspected

1. The analysis family is the complete 145-gene union encoded by `config/phase3_focus_gene_sets.gmt`. Repository lineage shows that this union preceded the outcome-dependent 145-to-61 filter. Because the project has no Git history, external registration, or trusted third-party timestamp, M18 calls it an **analysis-defined focused universe**, not a preregistered universe.
2. The primary synthesis uses GSE142025 C1, GSE96804, and GSE30528. GSE30528 is the prespecified glomerular contribution from source study GSE30122. GSE30529 is a secondary tubulointerstitial sensitivity dataset from the same source study and is not entered as a fourth independent study.
3. Random-effects estimation uses REML. Because only three studies contribute and their precisions differ, primary mean uncertainty uses modified Hartung-Knapp with `max(1, q)` and t degrees of freedom `k-1`; the unmodified HKSJ interval is also reported as a sensitivity. Prediction intervals use t degrees of freedom `k-2` when three studies are available.
4. BH correction retains all 145 genes in the family. Missing or non-estimable primary P values are assigned P=1 rather than silently reducing the multiplicity denominator.
5. External gene-level validation uses only GSE30528 and corrects across all 145 focused genes. Pathway validation permutes the 9-versus-13 disease labels at sample level and recomputes all genes jointly, preserving observed gene-gene correlation.
6. The 61-gene subset is retained only as a documented outcome-dependent subset (`n={selected_61}`); its conventional pooled P values/FDR are not the primary evidence.
7. Single-nucleus findings are not part of the primary statistical claim. They may be presented only as supplementary, provisional context until reference mapping, doublet detection, expanded QC, and donor-level pseudobulk are completed.

## Result closure

- Complete three-study mapping: {complete}/145 genes.
- Primary REML-HK FDR<0.05 within the 145-gene family: {significant}/145 genes.
- External GSE30528 mapping: {external_mapped}/145 genes.
- External GSE30528 FDR<0.05 within the 145-gene family: {external_fdr}/145 genes.
- Direction concordance among mapped external genes: {concordant}/{external_mapped}.
- Pathways with sample-label permutation FDR<0.05: {pathway_fdr}/9.

These counts describe an analysis-defined focused universe and one external source study; they are not evidence of replication across multiple independent validation studies.
"""
    (DOC_DIR / "M18_STATISTICAL_DECISION.md").write_text(decision, encoding="utf-8")

    methods = """# M18 statistical methods insert

The focused analysis family comprised the 145 unique genes in nine biologically defined complement/coagulation, vascular, inflammatory, extracellular-matrix, hypoxia, and transforming-growth-factor-beta gene sets. Repository lineage placed this union before the outcome-dependent filter that selected 61 genes; in the absence of an external registration timestamp, the union is described as analysis-defined rather than preregistered. All 145 genes were carried into inference, and missing or non-estimable tests were retained as P=1 for Benjamini-Hochberg adjustment.

For the primary synthesis, each source study contributed one Hedges' g standardized mean difference: advanced DKD versus control in GSE142025, DKD versus control in GSE96804, and the glomerular GSE30528 compartment of source study GSE30122. GSE30529 was not entered as an additional independent study because it is the tubulointerstitial companion series from GSE30122 and includes repeated donors. Between-study variance was estimated by restricted maximum likelihood. The pooled mean used inverse total-variance weights. Because only three studies contributed and their precisions differed, primary uncertainty used the modified Hartung-Knapp scale max(1, q) with a t distribution and k-1 degrees of freedom; the unmodified HKSJ interval was retained as a sensitivity. A 95% prediction interval based on the modified uncertainty was calculated with k-2 degrees of freedom when all three studies were estimable.

External validation was defined before inspecting M18 results as the glomerular GSE30528 contrast. Gene-level Welch tests were adjusted over the complete 145-gene family. For each of the nine pathways, discovery direction was fixed from an inverse-variance synthesis of GSE142025 and GSE96804. The pathway statistic was the mean direction-aligned Hedges' g in GSE30528. Disease labels were permuted 10,000 times while retaining the 9-versus-13 allocation; all pathway genes were recomputed jointly within each permutation so that sample structure and gene-gene correlation were preserved. One-sided permutation P values used the plus-one correction and were BH-adjusted across nine pathways.
"""
    (DOC_DIR / "M18_METHODS_INSERT.md").write_text(methods, encoding="utf-8")

    audit = pd.DataFrame(
        [
            {
                "artifact": path.relative_to(ROOT).as_posix(),
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
            for path in [
                TABLE_DIR / "focused_universe_145.csv",
                TABLE_DIR / "per_cohort_effect_sizes_145.csv",
                TABLE_DIR / "primary_reml_hk_145.csv",
                TABLE_DIR / "external_gse30528_gene_validation_145.csv",
                TABLE_DIR / "external_gse30528_pathway_permutation.csv",
            ]
        ]
    )
    audit.to_csv(TABLE_DIR / "output_manifest.csv", index=False)


def main() -> None:
    for directory in (TABLE_DIR, DOC_DIR, LOG_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc)
    inputs = input_paths()
    before = {path: sha256_file(path) for path in inputs}

    m8 = load_module("stage21_m8_frozen_for_m18", M8_PATH)
    candidates, pathway_map = read_focused_universe()
    per_cohort, _ = prepare_cohorts(m8, candidates)
    directions = discovery_direction_table(per_cohort)
    primary = run_primary_meta(candidates, per_cohort, directions)
    external = run_external_gene_validation(candidates, per_cohort, directions)
    pathway_results = run_pathway_permutation(pathway_map, directions)

    candidates.to_csv(TABLE_DIR / "focused_universe_145.csv", index=False)
    per_cohort.to_csv(TABLE_DIR / "per_cohort_effect_sizes_145.csv", index=False)
    directions.to_csv(TABLE_DIR / "discovery_direction_145.csv", index=False)
    primary.to_csv(TABLE_DIR / "primary_reml_hk_145.csv", index=False)
    external.to_csv(TABLE_DIR / "external_gse30528_gene_validation_145.csv", index=False)
    pathway_results.to_csv(TABLE_DIR / "external_gse30528_pathway_permutation.csv", index=False)

    after = {path: sha256_file(path) for path in inputs}
    hashes = pd.DataFrame(
        [
            {
                "path": path.relative_to(ROOT).as_posix(),
                "sha256_before": before[path],
                "sha256_after": after[path],
                "unchanged": before[path] == after[path],
            }
            for path in inputs
        ]
    )
    hashes.to_csv(TABLE_DIR / "input_hash_verification.csv", index=False)
    if not hashes["unchanged"].all():
        raise AssertionError("A frozen input changed during M18")

    settings = {
        "stage": "Stage 21-M18 statistical remediation",
        "execution_utc": started.isoformat(),
        "random_seed": RANDOM_SEED,
        "permutations": N_PERMUTATIONS,
        "analysis_family": 145,
        "primary_studies": PRIMARY_COHORTS,
        "primary_gse30122_estimand": "GSE30528 glomerular compartment",
        "secondary_same_study_compartment": "GSE30529 tubulointerstitium",
        "meta_model": "REML random effects; modified Hartung-Knapp primary t inference; unmodified HKSJ sensitivity",
        "prediction_interval": "t(k-2) when k>=3",
        "multiplicity": "BH over exactly 145 genes; missing/non-estimable P=1",
        "pathway_test": "GSE30528 sample-label permutation preserving 9-versus-13 allocation and gene correlation",
        "single_nucleus_role": "supplementary provisional context only",
    }
    (TABLE_DIR / "analysis_settings.json").write_text(
        json.dumps(settings, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    write_reports(candidates, primary, external, pathway_results, hashes)

    summary = {
        "complete_three_study_mapping": int(primary["complete_three_study_mapping"].sum()),
        "primary_fdr_lt_0_05": int(primary["fdr_lt_0_05"].sum()),
        "external_mapped": int(external["external_welch_p_value"].notna().sum()),
        "external_fdr_lt_0_05": int(external["external_fdr_lt_0_05"].sum()),
        "external_direction_concordant": int(
            external.loc[external["external_welch_p_value"].notna(), "direction_concordant"].sum()
        ),
        "pathway_fdr_lt_0_05": int(pathway_results["fdr_bh_9_pathways"].lt(0.05).sum()),
    }
    (LOG_DIR / "m18_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
