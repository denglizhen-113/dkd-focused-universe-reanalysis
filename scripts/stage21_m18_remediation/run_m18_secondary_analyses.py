from __future__ import annotations

import importlib.util
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, ttest_ind


ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = ROOT / "tables" / "stage21_m18_remediation"
M8_PATH = ROOT / "scripts" / "stage21_m8_validation" / "run_m8_analysis.py"
M18_PATH = ROOT / "scripts" / "stage21_m18_remediation" / "run_m18_remediation.py"
SHARED_DONORS = {"62", "67", "164", "168", "178", "76", "77", "81", "82"}
RANDOM_SEED = 20260815
N_PERMUTATIONS = 10_000


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def donor_from_title(title: str) -> str:
    match = re.search(r"Kidney\s+(\d+)", str(title), flags=re.IGNORECASE)
    return match.group(1) if match else ""


def prepare_donor_disjoint_gse30529(m8):
    prepared = load_module("m18_secondary_helpers", M18_PATH).prepare_external_matrix(m8, "GSE30529")
    metadata = pd.read_csv(ROOT / "data_raw" / "geo_soft" / "GSE30529" / "samples.csv")
    metadata["donor_id_inferred"] = metadata["title"].map(donor_from_title)
    metadata["shared_with_gse30528"] = metadata["donor_id_inferred"].isin(SHARED_DONORS)
    retained = metadata.loc[~metadata["shared_with_gse30528"]].copy()
    is_dkd = retained["characteristics"].str.contains(
        "diabetic kidney disease", case=False, na=False
    )
    prepared["case_ids"] = retained.loc[is_dkd, "gsm_accession"].astype(str).tolist()
    prepared["control_ids"] = retained.loc[~is_dkd, "gsm_accession"].astype(str).tolist()
    prepared["metadata"] = metadata
    prepared["retained_metadata"] = retained
    return prepared


def donor_disjoint_gene_results(m8, m18, candidates: pd.DataFrame, directions: pd.DataFrame):
    prepared = prepare_donor_disjoint_gse30529(m8)
    expression = prepared["expression"]
    direction_index = directions.set_index("gene_symbol")
    rows = []
    for gene in candidates["gene_symbol"].astype(str):
        base = {
            "gene_symbol": gene,
            "n_case": len(prepared["case_ids"]),
            "n_control": len(prepared["control_ids"]),
            "removed_shared_donors": ";".join(sorted(SHARED_DONORS, key=int)),
            "discovery_direction": direction_index.loc[gene, "discovery_direction"],
            "discovery_effect_fixed": direction_index.loc[gene, "discovery_effect_fixed"],
        }
        if gene not in expression.index:
            rows.append({**base, "mapping_status": "NOT_MAPPED"})
            continue
        case = expression.loc[gene, prepared["case_ids"]].to_numpy(float)
        control = expression.loc[gene, prepared["control_ids"]].to_numpy(float)
        effect = m8.hedges_g_from_samples(case, control)
        test = ttest_ind(case, control, equal_var=False, nan_policy="omit")
        direction = "UP" if effect["hedges_g"] > 0 else "DOWN" if effect["hedges_g"] < 0 else "NOT_ESTIMABLE"
        rows.append(
            {
                **base,
                "mapping_status": "MAPPED_TESTED",
                "case_mean": float(np.nanmean(case)),
                "control_mean": float(np.nanmean(control)),
                "raw_mean_difference": float(np.nanmean(case) - np.nanmean(control)),
                "hedges_g": effect["hedges_g"],
                "variance": effect["variance"],
                "se": effect["se"],
                "welch_p_value": float(test.pvalue),
                "validation_direction": direction,
                "direction_concordant": direction == base["discovery_direction"],
            }
        )
    frame = pd.DataFrame(rows)
    frame["p_for_145_family"] = frame["welch_p_value"].fillna(1.0)
    frame["fdr_bh_145_family"] = m18.bh_family(frame["p_for_145_family"], 145)
    frame["fdr_lt_0_05"] = frame["fdr_bh_145_family"].lt(0.05)
    frame["analysis_scope"] = (
        "secondary donor-disjoint tubulointerstitial sensitivity within GSE30122; "
        "not an additional independent validation study"
    )
    return frame, prepared


def donor_disjoint_pathways(m18, pathways, directions, prepared):
    expression = prepared["expression"]
    ordered_samples = [*prepared["case_ids"], *prepared["control_ids"]]
    expression = expression.loc[:, ordered_samples]
    direction_index = directions.set_index("gene_symbol")["discovery_direction"]
    usable_genes = [
        gene for gene in expression.index.astype(str)
        if gene in direction_index.index and direction_index.loc[gene] in {"UP", "DOWN"}
    ]
    values = expression.loc[usable_genes].to_numpy(float)
    signs = np.array([1.0 if direction_index.loc[g] == "UP" else -1.0 for g in usable_genes])
    observed_mask = np.array([sample in set(prepared["case_ids"]) for sample in ordered_samples])
    observed_effects = m18.hedges_g_matrix(values, observed_mask) * signs
    gene_to_index = {gene: index for index, gene in enumerate(usable_genes)}
    pathway_indexes = {
        pathway: np.array([gene_to_index[g] for g in genes if g in gene_to_index], dtype=int)
        for pathway, genes in pathways.items()
    }
    observed = {p: float(np.nanmean(observed_effects[idx])) for p, idx in pathway_indexes.items()}
    exceed = {pathway: 0 for pathway in pathways}
    rng = np.random.default_rng(RANDOM_SEED + 1)
    for _ in range(N_PERMUTATIONS):
        permuted_case = np.zeros(len(ordered_samples), dtype=bool)
        permuted_case[rng.choice(len(ordered_samples), observed_mask.sum(), replace=False)] = True
        permuted_effects = m18.hedges_g_matrix(values, permuted_case) * signs
        for pathway, indexes in pathway_indexes.items():
            if float(np.nanmean(permuted_effects[indexes])) >= observed[pathway] - 1e-15:
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
                "permutation_unit": "GSE30529 donor-disjoint sample labels; fixed 5-versus-8 allocation",
                "analysis_scope": "secondary same-source compartment sensitivity",
            }
        )
    frame = pd.DataFrame(rows)
    frame["fdr_bh_9_pathways"] = m18.bh_family(frame["one_sided_permutation_p"], 9)
    return frame.sort_values("fdr_bh_9_pathways", kind="stable").reset_index(drop=True)


def paired_compartment_correlations(m8, candidates: pd.DataFrame):
    glom = load_module("m18_paired_helper_g", M18_PATH).prepare_external_matrix(m8, "GSE30528")
    tub = load_module("m18_paired_helper_t", M18_PATH).prepare_external_matrix(m8, "GSE30529")
    glom_meta = pd.read_csv(ROOT / "data_raw" / "geo_soft" / "GSE30528" / "samples.csv")
    tub_meta = pd.read_csv(ROOT / "data_raw" / "geo_soft" / "GSE30529" / "samples.csv")
    glom_meta["donor"] = glom_meta["title"].map(donor_from_title)
    tub_meta["donor"] = tub_meta["title"].map(donor_from_title)
    glom_map = glom_meta.set_index("donor")["gsm_accession"].to_dict()
    tub_map = tub_meta.set_index("donor")["gsm_accession"].to_dict()
    donors = sorted(SHARED_DONORS, key=int)
    rows = []
    for gene in candidates["gene_symbol"].astype(str):
        if gene not in glom["expression"].index or gene not in tub["expression"].index:
            rows.append({"gene_symbol": gene, "n_paired_donors": 0, "pearson_r": math.nan, "p_value": math.nan})
            continue
        x = np.array([glom["expression"].loc[gene, glom_map[d]] for d in donors], dtype=float)
        y = np.array([tub["expression"].loc[gene, tub_map[d]] for d in donors], dtype=float)
        valid = np.isfinite(x) & np.isfinite(y)
        if valid.sum() >= 3 and np.std(x[valid], ddof=1) > 0 and np.std(y[valid], ddof=1) > 0:
            result = pearsonr(x[valid], y[valid])
            r_value, p_value = float(result.statistic), float(result.pvalue)
        else:
            r_value, p_value = math.nan, math.nan
        rows.append(
            {
                "gene_symbol": gene,
                "n_paired_donors": int(valid.sum()),
                "paired_donor_ids": ";".join(donors),
                "pearson_r": r_value,
                "p_value": p_value,
                "interpretation_limit": (
                    "expression correlation across inferred cross-compartment donor pairs; "
                    "not an estimate of covariance between cohort-level disease effects"
                ),
            }
        )
    return pd.DataFrame(rows)


def wgcna_integration(candidates: pd.DataFrame, external: pd.DataFrame):
    modules = pd.read_csv(ROOT / "tables" / "phase5_wgcna" / "GSE96804_gene_module_assignment.csv")
    modules["gene_symbol"] = modules["gene_symbol"].astype(str).str.upper()
    columns = [
        "gene_symbol", "module", "gene_module_membership", "gene_significance_for_DKD",
        "is_module_hub_candidate", "module_trait_correlation_with_DKD", "module_trait_pvalue",
    ]
    base = candidates[["gene_symbol", "theme"]].copy()
    merged = base.merge(modules[columns], on="gene_symbol", how="left", validate="one_to_one")
    ext = external[["gene_symbol", "direction_concordant", "external_fdr_lt_0_05", "hedges_g"]].rename(
        columns={"hedges_g": "gse30528_hedges_g"}
    )
    merged = merged.merge(ext, on="gene_symbol", how="left", validate="one_to_one")
    merged["mapped_to_wgcna_top_5000"] = merged["module"].notna()
    merged["in_black_module"] = merged["module"].eq("black")
    merged["descriptive_class"] = np.select(
        [
            merged["in_black_module"] & merged["external_fdr_lt_0_05"].fillna(False),
            merged["in_black_module"] & ~merged["external_fdr_lt_0_05"].fillna(False),
            ~merged["in_black_module"] & merged["external_fdr_lt_0_05"].fillna(False),
        ],
        [
            "BLACK_AND_GSE30528_FDR",
            "BLACK_ONLY",
            "GSE30528_FDR_OUTSIDE_BLACK",
        ],
        default="NEITHER_OR_NOT_IN_WGCNA_UNIVERSE",
    )
    merged["interpretation_limit"] = (
        "WGCNA was derived inside GSE96804 and is descriptive internal corroboration, not validation"
    )
    summary = (
        merged.groupby("descriptive_class", dropna=False)
        .agg(n_genes=("gene_symbol", "size"), genes=("gene_symbol", lambda x: ";".join(sorted(x))))
        .reset_index()
    )
    return merged, summary


def compartment_comparison(glom: pd.DataFrame, tub: pd.DataFrame):
    left = glom[[
        "gene_symbol", "external_mapping_status", "direction_concordant",
        "external_fdr_bh_145_family", "external_fdr_lt_0_05", "hedges_g",
    ]].rename(columns=lambda c: f"gse30528_{c}" if c != "gene_symbol" else c)
    right = tub[[
        "gene_symbol", "mapping_status", "direction_concordant", "fdr_bh_145_family",
        "fdr_lt_0_05", "hedges_g",
    ]].rename(columns=lambda c: f"gse30529_disjoint_{c}" if c != "gene_symbol" else c)
    merged = left.merge(right, on="gene_symbol", how="outer", validate="one_to_one")
    g = merged["gse30528_external_fdr_lt_0_05"].fillna(False)
    t = merged["gse30529_disjoint_fdr_lt_0_05"].fillna(False)
    merged["significance_switch_class"] = np.select(
        [g & t, g & ~t, ~g & t],
        ["BOTH", "GSE30528_ONLY", "GSE30529_DISJOINT_ONLY"],
        default="NEITHER",
    )
    merged["display_color_code"] = merged["significance_switch_class"].map(
        {"BOTH": "PURPLE", "GSE30528_ONLY": "RED", "GSE30529_DISJOINT_ONLY": "BLUE", "NEITHER": "GRAY"}
    )
    return merged


def primary_outlier_exclusion_sensitivity(m8, m18, candidates, directions):
    diagnostics_path = TABLE_DIR / "sample_correlation_outlier_diagnostics.csv"
    if not diagnostics_path.exists():
        return pd.DataFrame()
    diagnostics = pd.read_csv(diagnostics_path)
    flagged = diagnostics.loc[diagnostics["predefined_outlier_flag"].eq(True)]
    flagged_by_dataset = {
        "GSE142025_C1": set(flagged.loc[flagged["dataset"].eq("GSE142025 C1"), "sample_id"]),
        "GSE96804": set(flagged.loc[flagged["dataset"].eq("GSE96804"), "sample_id"]),
        "GSE30528": set(flagged.loc[flagged["dataset"].eq("GSE30528"), "sample_id"]),
    }
    _, context = m18.prepare_cohorts(m8, candidates)
    filtered_definitions = []
    for definition in context["definitions"]:
        if definition["cohort"] not in m18.PRIMARY_COHORTS:
            continue
        remove = flagged_by_dataset.get(definition["cohort"], set())
        updated = dict(definition)
        updated["case_ids"] = [sample for sample in definition["case_ids"] if sample not in remove]
        updated["control_ids"] = [sample for sample in definition["control_ids"] if sample not in remove]
        updated["contrast"] = definition["contrast"] + "_correlation_outliers_excluded"
        filtered_definitions.append(updated)
    effects = m8.build_per_cohort_effect_sizes(candidates, filtered_definitions)
    result = m18.run_primary_meta(candidates, effects, directions)
    result["excluded_sample_ids"] = ";".join(sorted(set(flagged["sample_id"].astype(str))))
    result["excluded_n"] = len(flagged)
    result["sensitivity_scope"] = (
        "excludes samples with predefined robust-z median-correlation flag; "
        "potential biological and technical causes are not distinguishable"
    )
    return result


def main():
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    m8 = load_module("stage21_m8_secondary", M8_PATH)
    m18 = load_module("stage21_m18_secondary", M18_PATH)
    candidates, pathways = m18.read_focused_universe()
    directions = pd.read_csv(TABLE_DIR / "discovery_direction_145.csv")
    external = pd.read_csv(TABLE_DIR / "external_gse30528_gene_validation_145.csv")

    tub, prepared = donor_disjoint_gene_results(m8, m18, candidates, directions)
    tub_pathways = donor_disjoint_pathways(m18, pathways, directions, prepared)
    paired = paired_compartment_correlations(m8, candidates)
    wgcna, wgcna_summary = wgcna_integration(candidates, external)
    comparison = compartment_comparison(external, tub)
    outlier_sensitivity = primary_outlier_exclusion_sensitivity(m8, m18, candidates, directions)

    tub.to_csv(TABLE_DIR / "gse30529_donor_disjoint_gene_validation_145.csv", index=False)
    tub_pathways.to_csv(TABLE_DIR / "gse30529_donor_disjoint_pathway_permutation.csv", index=False)
    paired.to_csv(TABLE_DIR / "paired_compartment_expression_correlations_145.csv", index=False)
    wgcna.to_csv(TABLE_DIR / "wgcna_descriptive_integration_145.csv", index=False)
    wgcna_summary.to_csv(TABLE_DIR / "wgcna_descriptive_summary.csv", index=False)
    comparison.to_csv(TABLE_DIR / "compartment_significance_switch_summary_145.csv", index=False)
    if not outlier_sensitivity.empty:
        outlier_sensitivity.to_csv(
            TABLE_DIR / "primary_outlier_exclusion_sensitivity_145.csv", index=False
        )

    mapped = tub["welch_p_value"].notna()
    print(
        "gse30529_disjoint "
        f"mapped={int(mapped.sum())}/145 "
        f"direction_concordant={int(tub.loc[mapped, 'direction_concordant'].sum())}/{int(mapped.sum())} "
        f"fdr={int(tub['fdr_lt_0_05'].sum())}/145 "
        f"pathways_fdr={int(tub_pathways['fdr_bh_9_pathways'].lt(0.05).sum())}/9"
    )
    valid_r = paired["pearson_r"].dropna()
    print(f"paired_compartment_r_median={valid_r.median():.4f} n_genes={len(valid_r)}")
    if not outlier_sensitivity.empty:
        print(
            "outlier_exclusion_primary "
            f"complete={int(outlier_sensitivity['complete_three_study_mapping'].sum())}/145 "
            f"fdr={int(outlier_sensitivity['fdr_lt_0_05'].sum())}/145"
        )


if __name__ == "__main__":
    main()
