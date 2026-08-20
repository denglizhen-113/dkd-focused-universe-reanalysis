#!/usr/bin/env python3
"""M19 compartment-stratified DKD analysis with canonical Reactome pathways.

This analysis avoids a cross-compartment pooled estimand. The primary gene-level
meta-analysis is restricted to three glomerular source studies. Whole-kidney,
tubulointerstitial, and interstitial-only estimates are reported separately.
Canonical Reactome pathways replace the manually curated project groupings for
the primary pathway family. Sample-label tests are two-sided and maxT-adjusted;
they remain explicitly conditional on exchangeability of the archived samples.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import itertools
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import t as t_dist
from scipy.stats import ttest_ind


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "tables" / "stage21_m19_scientific_reports_revision"
FIG = ROOT / "figures" / "stage21_m19_scientific_reports_revision"
DOC = ROOT / "docs" / "stage21_m19_scientific_reports_revision"
M8_PATH = ROOT / "scripts" / "stage21_m8_validation" / "run_m8_analysis.py"
M18_PATH = ROOT / "scripts" / "stage21_m18_remediation" / "run_m18_remediation.py"
M9_GSE111154_PATH = ROOT / "scripts" / "stage21_m9_gse111154" / "run_gse111154_rerun.py"
REACTOME_GMT = (
    ROOT
    / "data_raw"
    / "reference_gene_sets_m19"
    / "reactome_gmt"
    / "ReactomePathways.gmt"
)
GOA_GAF = ROOT / "data_raw" / "reference_gene_sets_m19" / "goa_human.gaf.gz"
SEED = 20260821
MONTE_CARLO_PERMUTATIONS = 10_000
PRIMARY_GLOMERULAR_COHORTS = ["GSE96804", "GSE30528", "GSE104948_H7"]

CANONICAL_PATHWAYS = {
    "Complement cascade": "R-HSA-166658",
    "Coagulation pathway": "R-HSA-9769740",
    "Cell surface interactions at the vascular wall": "R-HSA-202733",
    "Chemokine receptors bind chemokines": "R-HSA-380108",
    "Extracellular matrix organization": "R-HSA-1474244",
    "Cellular response to hypoxia": "R-HSA-1234174",
    "Signaling by TGF-beta Receptor Complex": "R-HSA-170834",
}

COHORT_SEEDS = {
    "GSE96804": SEED + 1,
    "GSE30528": SEED + 2,
    "GSE104948_H7": SEED + 3,
    "GSE142025_advanced": SEED + 4,
    "GSE166239": SEED + 5,
    "GSE30529": SEED + 6,
    "GSE104954_H7": SEED + 7,
    "GSE163603_interstitium": SEED + 8,
    "GSE1009_donor_averaged": SEED + 9,
    "GSE111154": SEED + 10,
    "GSE199838": SEED + 11,
}


@dataclass
class Cohort:
    name: str
    source_study: str
    compartment: str
    platform: str
    expression: pd.DataFrame
    case_ids: list[str]
    control_ids: list[str]
    processing: str
    control_definition: str
    covariates_available: str


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_human_symbols() -> set[str]:
    symbols: set[str] = set()
    with gzip.open(GOA_GAF, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) > 2:
                symbol = fields[2].upper()
                if re.fullmatch(r"[A-Z0-9][A-Z0-9.-]*", symbol):
                    symbols.add(symbol)
    return symbols


def load_canonical_pathways() -> dict[str, list[str]]:
    human_symbols = load_human_symbols()
    by_id: dict[str, tuple[str, list[str]]] = {}
    for line in REACTOME_GMT.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) >= 3:
            by_id[fields[1]] = (fields[0], fields[2:])

    pathways: dict[str, list[str]] = {}
    provenance_rows = []
    for expected_name, stable_id in CANONICAL_PATHWAYS.items():
        actual_name, raw_genes = by_id[stable_id]
        if actual_name != expected_name:
            raise AssertionError(f"Reactome name drift for {stable_id}: {actual_name}")
        genes = sorted({gene.upper() for gene in raw_genes if gene.upper() in human_symbols})
        if len(genes) < 5:
            raise AssertionError(f"Too few human genes for {stable_id}")
        pathways[stable_id] = genes
        provenance_rows.append(
            {
                "reactome_id": stable_id,
                "reactome_name": actual_name,
                "human_gene_count_after_GOA_filter": len(genes),
                "genes": ";".join(genes),
                "source_url": "https://reactome.org/download/current/ReactomePathways.gmt.zip",
                "source_sha256": sha256(
                    ROOT / "data_raw" / "reference_gene_sets_m19" / "ReactomePathways.gmt.zip"
                ),
                "human_symbol_filter": "GOA human GAF current release",
                "human_symbol_source_sha256": sha256(GOA_GAF),
                "retrieval_date": "2026-08-21",
                "inferential_role": "fixed seven-pathway primary family",
            }
        )
    pd.DataFrame(provenance_rows).to_csv(OUT / "canonical_pathway_provenance.csv", index=False)
    with (OUT / "canonical_pathways.gmt").open("w", encoding="utf-8", newline="\n") as handle:
        for stable_id, genes in pathways.items():
            name = next(name for name, sid in CANONICAL_PATHWAYS.items() if sid == stable_id)
            handle.write("\t".join([stable_id, name, *genes]) + "\n")
    return pathways


def load_gene_matrix(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    key = "gene_symbol" if "gene_symbol" in frame.columns else frame.columns[0]
    frame[key] = frame[key].astype(str).str.upper()
    frame = frame.drop_duplicates(key, keep="first").set_index(key)
    return frame.apply(pd.to_numeric, errors="coerce")


def parse_platform_symbol_map(soft_path: Path) -> pd.Series:
    in_table = False
    header: list[str] | None = None
    mapping: dict[str, str] = {}
    with gzip.open(soft_path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith("!platform_table_begin"):
                in_table = True
                header = None
                continue
            if line.startswith("!platform_table_end"):
                break
            if not in_table:
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                continue
            row = dict(zip(header, fields))
            probe = row.get("ID", "")
            symbol = row.get("Gene Symbol", row.get("Symbol", ""))
            symbol = symbol.split(" /// ")[0].split(";")[0].strip().upper()
            if probe and re.fullmatch(r"[A-Z0-9][A-Z0-9.-]*", symbol):
                mapping[probe] = symbol
    return pd.Series(mapping, name="gene_symbol", dtype=str)


def aggregate_probe_matrix(
    probe_matrix: pd.DataFrame, mapping: pd.DataFrame, aggregation: str
) -> pd.DataFrame:
    mapping = mapping.loc[mapping["probe_id"].isin(probe_matrix.index)].drop_duplicates()
    if aggregation == "highest_mean":
        means = probe_matrix.mean(axis=1, skipna=True)
        mapping = mapping.assign(_mean=mapping["probe_id"].map(means))
        selected = (
            mapping.sort_values(["gene_symbol", "_mean", "probe_id"], ascending=[True, False, True])
            .drop_duplicates("gene_symbol", keep="first")
        )
        result = probe_matrix.loc[selected["probe_id"]].copy()
        result.index = selected["gene_symbol"].to_numpy()
    elif aggregation == "median":
        expanded = probe_matrix.loc[mapping["probe_id"]].copy()
        expanded.index = mapping["gene_symbol"].to_numpy()
        result = expanded.groupby(level=0).median(numeric_only=True)
    else:
        raise ValueError(aggregation)
    result.index.name = "gene_symbol"
    return result.sort_index()


def load_soft_gene_matrix(
    m8, accession: str, aggregation: str = "highest_mean", m9_111154=None
) -> pd.DataFrame:
    soft = ROOT / "data_raw" / "geo_soft" / accession / f"{accession}_family.soft.gz"
    probe_matrix, probe_map, _ = m8.parse_geo_soft(soft)
    if probe_map.empty:
        probe_map = parse_platform_symbol_map(soft)
    if probe_map.empty and accession == "GSE111154" and m9_111154 is not None:
        probe_matrix, mapping, _ = m9_111154.parse_gpl17586_soft(soft)
        return aggregate_probe_matrix(probe_matrix, mapping, aggregation)
    common = probe_matrix.index.intersection(probe_map.index)
    annotated = probe_matrix.loc[common].copy()
    symbols = probe_map.loc[common].astype(str).str.upper()
    if aggregation == "highest_mean":
        means = annotated.mean(axis=1, skipna=True)
        order = means.sort_values(ascending=False, kind="stable").index
        ordered_symbols = symbols.loc[order]
        selected = order[~ordered_symbols.duplicated(keep="first")]
        result = annotated.loc[selected].copy()
        result.index = symbols.loc[selected].to_numpy()
    elif aggregation == "median":
        result = annotated.assign(_gene=symbols.to_numpy()).groupby("_gene").median(numeric_only=True)
    else:
        raise ValueError(aggregation)
    result.index.name = "gene_symbol"
    return result.sort_index()


def prepare_cohorts(m8, m9_111154) -> dict[str, Cohort]:
    cohorts: dict[str, Cohort] = {}

    expr = load_gene_matrix(ROOT / "data_processed" / "GSE96804" / "expression_matrix_clean.csv")
    meta = pd.read_csv(ROOT / "data_processed" / "GSE96804" / "sample_annotation_clean.csv")
    cohorts["GSE96804"] = Cohort(
        "GSE96804", "GSE96804", "glomerular", "GPL17586 microarray", expr,
        meta.loc[meta["disease_group_clean"].eq("DKD"), "sample_id"].astype(str).tolist(),
        meta.loc[meta["disease_group_clean"].eq("Control"), "sample_id"].astype(str).tolist(),
        "project-frozen gene-level normalized matrix", "unaffected tumor-nephrectomy glomeruli",
        "sex complete; batch and age unavailable",
    )

    for accession, compartment in (("GSE30528", "glomerular"), ("GSE30529", "tubulointerstitial")):
        expr = load_soft_gene_matrix(m8, accession)
        meta = pd.read_csv(ROOT / "data_raw" / "geo_soft" / accession / "samples.csv")
        case = meta["characteristics"].str.contains("diabetic kidney disease", case=False, na=False)
        cohorts[accession] = Cohort(
            accession, "GSE30122", compartment, "GPL571 microarray", expr,
            meta.loc[case, "gsm_accession"].astype(str).tolist(),
            meta.loc[~case, "gsm_accession"].astype(str).tolist(),
            "GEO normalized SOFT; highest all-sample-mean probe per gene",
            "living-donor kidney compartment", "no common demographic or batch covariates",
        )

    for accession, compartment, token in (
        ("GSE104948", "glomerular", "Glom"),
        ("GSE104954", "tubulointerstitial", "Tub"),
    ):
        expr = load_soft_gene_matrix(m8, accession)
        meta = pd.read_csv(ROOT / "data_raw" / "geo_soft" / accession / "samples.csv")
        h7 = meta["title"].str.startswith(f"H7-{token}-", na=False)
        case = h7 & meta["title"].str.contains("-DN", case=False, na=False)
        control = h7 & meta["title"].str.contains("-LD", case=False, na=False)
        name = f"{accession}_H7"
        cohorts[name] = Cohort(
            name, "ERCB_H7", compartment, "custom-CDF microarray", expr,
            meta.loc[case, "gsm_accession"].astype(str).tolist(),
            meta.loc[control, "gsm_accession"].astype(str).tolist(),
            "GEO normalized SOFT; H7 batch only; highest all-sample-mean probe per gene",
            "H7 living-donor compartment", "batch controlled by H7 restriction; demographics unavailable",
        )

    expr = load_gene_matrix(ROOT / "data_processed" / "GSE142025" / "expression_matrix_clean.csv")
    meta = pd.read_csv(ROOT / "data_processed" / "GSE142025" / "sample_annotation_clean.csv")
    cohorts["GSE142025_advanced"] = Cohort(
        "GSE142025_advanced", "GSE142025", "whole/cortical kidney", "GPL20301 RNA-seq", expr,
        meta.loc[meta["disease_group_clean"].eq("Advanced_DKD"), "sample_id"].astype(str).tolist(),
        meta.loc[meta["disease_group_clean"].eq("Control"), "sample_id"].astype(str).tolist(),
        "project-frozen log2 quantile-normalized gene matrix", "nephrectomy normal kidney",
        "no common demographic or batch covariates",
    )

    expr = load_gene_matrix(ROOT / "data_processed" / "GSE166239" / "expression_matrix_clean.csv")
    meta = pd.read_csv(ROOT / "data_processed" / "GSE166239" / "sample_annotation_clean.csv")
    cohorts["GSE166239"] = Cohort(
        "GSE166239", "GSE166239", "whole/cortical kidney", "GPL18573 RNA-seq", expr,
        meta.loc[meta["disease_group_clean"].eq("T2DN"), "sample_id"].astype(str).tolist(),
        meta.loc[meta["disease_group_clean"].eq("Control"), "sample_id"].astype(str).tolist(),
        "project-frozen logCPM-normalized gene matrix", "non-diseased renal tissue",
        "age and sex complete; small n=6 per group",
    )

    counts = pd.read_csv(
        ROOT / "data_raw" / "geo_supp_m19" / "GSE163603_processed_data.csv.gz",
        compression="gzip",
    ).set_index("Geneid")
    counts.index = counts.index.astype(str).str.upper()
    case_cols = [c for c in counts.columns if c.startswith("DM") and c.endswith("Interstitium")]
    control_cols = [c for c in counts.columns if not c.startswith("DM") and c.endswith("Interstitium")]
    use = case_cols + control_cols
    cpm = counts[use].divide(counts[use].sum(axis=0), axis=1) * 1_000_000.0
    expr = np.log2(cpm + 0.5)
    cohorts["GSE163603_interstitium"] = Cohort(
        "GSE163603_interstitium", "GSE163603", "interstitium only", "GPL16791 RNA-seq", expr,
        case_cols, control_cols, "author processed counts; log2(CPM+0.5)",
        "laser-microdissected nephrectomy interstitium",
        "sex available; all DKD samples male; disease and sex partly confounded",
    )

    # Small-sample studies are retained as sensitivities instead of being excluded
    # by an arbitrary minimum-n threshold.
    expr = load_soft_gene_matrix(m8, "GSE1009")
    meta = pd.read_csv(ROOT / "data_raw" / "geo_soft" / "GSE1009" / "samples.csv")
    title_to_gsm = dict(zip(meta["title"], meta["gsm_accession"]))
    donor_averaged = pd.DataFrame(
        {
            "Control_1": expr[[title_to_gsm["Control 1a"], title_to_gsm["Control 1b"]]].mean(axis=1),
            "Control_2": expr[title_to_gsm["Control 2"]],
            "Diabetes_1": expr[[title_to_gsm["Diabetes 1a"], title_to_gsm["Diabetes 1b"]]].mean(axis=1),
            "Diabetes_2": expr[title_to_gsm["Diabetes 2"]],
        }
    )
    cohorts["GSE1009_donor_averaged"] = Cohort(
        "GSE1009_donor_averaged", "GSE1009", "glomerular", "legacy microarray", donor_averaged,
        ["Diabetes_1", "Diabetes_2"], ["Control_1", "Control_2"],
        "GEO normalized SOFT; apparent a/b technical replicates averaged at donor label level",
        "non-diabetic glomeruli", "two apparent donors per group; descriptive sensitivity only",
    )

    expr = load_soft_gene_matrix(m8, "GSE111154", m9_111154=m9_111154)
    meta = pd.read_csv(ROOT / "data_raw" / "geo_soft" / "GSE111154" / "samples.csv")
    case = meta["characteristics"].str.contains("early diabetic nephropathy", case=False, na=False)
    cohorts["GSE111154"] = Cohort(
        "GSE111154", "GSE111154", "whole/cortical kidney", "GPL17586 microarray", expr,
        meta.loc[case, "gsm_accession"].astype(str).tolist(),
        meta.loc[~case, "gsm_accession"].astype(str).tolist(),
        "GEO normalized SOFT; corrected GPL17586 gene-assignment parser",
        "postmortem non-diabetic kidney", "four samples per group; descriptive sensitivity only",
    )

    raw = pd.read_csv(
        ROOT / "data_raw" / "geo_supp_m19" / "GSE199838_all_gene.xls.gz",
        sep="\t", compression="gzip",
    )
    sample_cols = ["N-4-Y", "N-1-H", "N-2-G", "FN-1-Y", "FN-2-L", "FN-3-C"]
    raw["gene"] = raw["gene"].astype(str).str.upper()
    raw = raw.drop_duplicates("gene", keep="first").set_index("gene")
    counts = raw[sample_cols].apply(pd.to_numeric, errors="coerce")
    cpm = counts.divide(counts.sum(axis=0), axis=1) * 1_000_000.0
    expr = np.log2(cpm + 0.5)
    cohorts["GSE199838"] = Cohort(
        "GSE199838", "GSE199838", "whole/cortical kidney", "RNA-seq", expr,
        ["FN-1-Y", "FN-2-L", "FN-3-C"], ["N-4-Y", "N-1-H", "N-2-G"],
        "author gene-count table; log2(CPM+0.5)",
        "non-diabetic renal-cancer cortex", "three samples per group; descriptive sensitivity only",
    )
    return cohorts


def hedges_effects(expression: pd.DataFrame, case_ids: list[str], control_ids: list[str]) -> pd.DataFrame:
    case = expression[case_ids].to_numpy(float)
    control = expression[control_ids].to_numpy(float)
    n1, n0 = case.shape[1], control.shape[1]
    mean1 = np.nanmean(case, axis=1)
    mean0 = np.nanmean(control, axis=1)
    sd1 = np.nanstd(case, axis=1, ddof=1)
    sd0 = np.nanstd(control, axis=1, ddof=1)
    pooled_var = ((n1 - 1) * sd1**2 + (n0 - 1) * sd0**2) / (n1 + n0 - 2)
    pooled_sd = np.sqrt(pooled_var)
    correction = 1.0 - 3.0 / (4.0 * (n1 + n0) - 9.0)
    g = correction * (mean1 - mean0) / pooled_sd
    variance = (n1 + n0) / (n1 * n0) + g**2 / (2.0 * (n1 + n0 - 2.0))
    g[~np.isfinite(g)] = np.nan
    variance[~np.isfinite(g)] = np.nan
    p = ttest_ind(case, control, axis=1, equal_var=False, nan_policy="omit").pvalue
    return pd.DataFrame(
        {
            "gene_symbol": expression.index,
            "n_case": n1,
            "n_control": n0,
            "case_mean": mean1,
            "control_mean": mean0,
            "raw_mean_difference": mean1 - mean0,
            "hedges_g": g,
            "variance": variance,
            "welch_p_value": p,
        }
    )


def bh(values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(values), dtype=float)
    p = np.where(np.isfinite(p), np.clip(p, 0.0, 1.0), 1.0)
    order = np.argsort(p, kind="stable")
    ranked = p[order]
    q = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty(len(p), dtype=float)
    out[order] = np.minimum(q, 1.0)
    return out


def effect_vector(values: np.ndarray, case_mask: np.ndarray) -> np.ndarray:
    case = values[:, case_mask]
    control = values[:, ~case_mask]
    n1, n0 = case.shape[1], control.shape[1]
    mean1, mean0 = np.mean(case, axis=1), np.mean(control, axis=1)
    sd1, sd0 = np.std(case, axis=1, ddof=1), np.std(control, axis=1, ddof=1)
    pooled = np.sqrt(((n1 - 1) * sd1**2 + (n0 - 1) * sd0**2) / (n1 + n0 - 2))
    correction = 1.0 - 3.0 / (4.0 * (n1 + n0) - 9.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        g = correction * (mean1 - mean0) / pooled
    return np.where(np.isfinite(g), g, np.nan)


def allocation_masks(n: int, n_case: int, rng: np.random.Generator):
    total = math.comb(n, n_case)
    if total <= MONTE_CARLO_PERMUTATIONS:
        for indexes in itertools.combinations(range(n), n_case):
            mask = np.zeros(n, dtype=bool)
            mask[list(indexes)] = True
            yield mask, "exact", total
    else:
        for _ in range(MONTE_CARLO_PERMUTATIONS):
            mask = np.zeros(n, dtype=bool)
            mask[rng.choice(n, size=n_case, replace=False)] = True
            yield mask, "monte_carlo", MONTE_CARLO_PERMUTATIONS


def pathway_permutation(cohort: Cohort, pathways: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_ids = cohort.case_ids + cohort.control_ids
    union = sorted({g for genes in pathways.values() for g in genes if g in cohort.expression.index})
    matrix = cohort.expression.loc[union, sample_ids]
    values = matrix.to_numpy(float)
    observed_mask = np.array([s in set(cohort.case_ids) for s in sample_ids], dtype=bool)
    observed_g = effect_vector(values, observed_mask)
    gene_index = {gene: i for i, gene in enumerate(union)}
    indexes = {
        sid: np.array([gene_index[g] for g in genes if g in gene_index], dtype=int)
        for sid, genes in pathways.items()
    }
    observed = {sid: float(np.nanmean(observed_g[idx])) for sid, idx in indexes.items()}
    exceed = {sid: 0 for sid in pathways}
    max_exceed = {sid: 0 for sid in pathways}
    rng = np.random.default_rng(COHORT_SEEDS[cohort.name])
    method = ""
    denominator = 0
    for mask, method, denominator in allocation_masks(len(sample_ids), len(cohort.case_ids), rng):
        perm_g = effect_vector(values, mask)
        stats = {sid: float(np.nanmean(perm_g[idx])) for sid, idx in indexes.items()}
        max_abs = max(abs(value) for value in stats.values())
        for sid in pathways:
            if abs(stats[sid]) >= abs(observed[sid]) - 1e-15:
                exceed[sid] += 1
            if max_abs >= abs(observed[sid]) - 1e-15:
                max_exceed[sid] += 1

    rows = []
    loo_rows = []
    for sid, genes in pathways.items():
        idx = indexes[sid]
        if method == "exact":
            p = exceed[sid] / denominator
            max_p = max_exceed[sid] / denominator
        else:
            p = (exceed[sid] + 1) / (denominator + 1)
            max_p = (max_exceed[sid] + 1) / (denominator + 1)
        variances = []
        effects = []
        for gene in genes:
            if gene in gene_index:
                i = gene_index[gene]
                if np.isfinite(observed_g[i]):
                    effects.append(observed_g[i])
        for gene in genes:
            if gene not in gene_index or len(effects) <= 1:
                continue
            keep = [g for g in genes if g in gene_index and g != gene]
            loo = float(np.nanmean(observed_g[[gene_index[g] for g in keep]]))
            loo_rows.append(
                {
                    "cohort": cohort.name,
                    "reactome_id": sid,
                    "removed_gene": gene,
                    "full_mean_hedges_g": observed[sid],
                    "leave_one_gene_out_mean_hedges_g": loo,
                    "absolute_change": abs(loo - observed[sid]),
                }
            )
        rows.append(
            {
                "cohort": cohort.name,
                "source_study": cohort.source_study,
                "compartment": cohort.compartment,
                "reactome_id": sid,
                "reactome_name": next(name for name, value in CANONICAL_PATHWAYS.items() if value == sid),
                "genes_in_reference_set": len(genes),
                "genes_mapped": len(idx),
                "observed_mean_hedges_g": observed[sid],
                "direction": "UP" if observed[sid] > 0 else "DOWN" if observed[sid] < 0 else "ZERO",
                "permutation_method": method,
                "permutation_allocations": denominator,
                "two_sided_permutation_p": p,
                "maxT_fwer_p": max_p,
                "exchangeability_limit": (
                    "Conditional on exchangeability of archived case/control labels; maxT controls the "
                    "fixed seven-pathway family but cannot remove unrecorded confounding."
                ),
            }
        )
    result = pd.DataFrame(rows)
    result["bh_fdr_7_pathways"] = bh(result["two_sided_permutation_p"])
    return result, pd.DataFrame(loo_rows)


def meta_table(m18, per_cohort: pd.DataFrame, cohort_names: list[str], universe: list[str], label: str) -> pd.DataFrame:
    rows = []
    for gene in universe:
        sub = per_cohort.loc[
            per_cohort["cohort"].isin(cohort_names) & per_cohort["gene_symbol"].eq(gene)
        ]
        valid = sub.loc[np.isfinite(sub["hedges_g"]) & np.isfinite(sub["variance"])]
        meta = m18.reml_hartung_knapp(valid["hedges_g"].to_numpy(float), valid["variance"].to_numpy(float))
        rows.append(
            {
                "gene_symbol": gene,
                "compartment_model": label,
                "eligible_cohorts": ";".join(cohort_names),
                "contributing_cohorts": ";".join(valid["cohort"].tolist()),
                **meta,
            }
        )
    frame = pd.DataFrame(rows)
    frame["p_value_for_family"] = frame["p_value_modified_hk"].fillna(1.0)
    frame["fdr_bh_canonical_union"] = bh(frame["p_value_for_family"])
    frame["fdr_lt_0_05"] = frame["fdr_bh_canonical_union"] < 0.05
    return frame


def sample_qc(cohorts: dict[str, Cohort]) -> pd.DataFrame:
    rows = []
    for cohort in cohorts.values():
        samples = cohort.case_ids + cohort.control_ids
        expr = cohort.expression[samples]
        variable = expr.var(axis=1).sort_values(ascending=False).head(min(1000, len(expr))).index
        corr = expr.loc[variable].corr()
        off = corr.mask(np.eye(len(corr), dtype=bool))
        med = off.median(axis=1)
        center = float(med.median())
        mad = float(np.median(np.abs(med - center)))
        z = (med - center) / (1.4826 * mad) if mad > 0 else pd.Series(0.0, index=med.index)
        for sample in samples:
            rows.append(
                {
                    "cohort": cohort.name,
                    "sample_id": sample,
                    "group": "DKD" if sample in set(cohort.case_ids) else "CONTROL",
                    "median_sample_correlation": med[sample],
                    "robust_z_median_correlation": z[sample],
                    "correlation_outlier_flag": bool(z[sample] < -3.0),
                    "rule": "robust z of median Pearson correlation < -3 using top 1000 variable genes",
                    "interpretation": (
                        "Screening flag only. A case-specific global expression pattern can lower correlations; "
                        "samples are not excluded without independent technical evidence."
                    ),
                }
            )
    return pd.DataFrame(rows)


def ols_disease_coefficient(expression: pd.DataFrame, design: pd.DataFrame) -> pd.DataFrame:
    """Return standardized disease coefficients from a fixed OLS design.

    The outcome is standardized gene by gene across the included samples. This is
    a sensitivity analysis for measured covariates, not a substitute for a
    prespecified covariate-adjusted primary model.
    """
    sample_ids = design.index.tolist()
    y = expression[sample_ids].to_numpy(float).T
    means = np.nanmean(y, axis=0)
    scales = np.nanstd(y, axis=0, ddof=1)
    valid_gene = np.isfinite(scales) & (scales > 0)
    y = (y - means) / scales
    x = design.to_numpy(float)
    # The project runtime has an unstable external BLAS/LAPACK installation on
    # this Windows host. The design has at most four columns, so a transparent
    # pivoted Gauss-Jordan inverse plus einsum avoids an opaque native crash.
    xtx = np.einsum("ni,nj->ij", x, x)
    augmented = np.concatenate([xtx.copy(), np.eye(xtx.shape[0])], axis=1)
    for column in range(xtx.shape[0]):
        pivot = column + int(np.argmax(np.abs(augmented[column:, column])))
        if abs(augmented[pivot, column]) < 1e-12:
            raise ValueError("Covariate design is rank deficient")
        if pivot != column:
            augmented[[column, pivot]] = augmented[[pivot, column]]
        augmented[column] /= augmented[column, column]
        for row in range(xtx.shape[0]):
            if row != column:
                augmented[row] -= augmented[row, column] * augmented[column]
    xtx_inv = augmented[:, xtx.shape[0] :]
    rank = xtx.shape[0]
    df = x.shape[0] - rank
    if df <= 0 or "disease" not in design.columns:
        raise ValueError("Covariate design has no residual degrees of freedom or disease column")
    xty = np.einsum("ni,ng->ig", x, y)
    beta = np.einsum("ij,jg->ig", xtx_inv, xty)
    residual = y - np.einsum("ni,ig->ng", x, beta)
    sigma2 = np.nansum(residual**2, axis=0) / df
    disease_index = design.columns.get_loc("disease")
    se = np.sqrt(np.maximum(sigma2 * xtx_inv[disease_index, disease_index], 0.0))
    coefficient = beta[disease_index]
    with np.errstate(divide="ignore", invalid="ignore"):
        statistic = coefficient / se
    p_value = 2.0 * t_dist.sf(np.abs(statistic), df)
    coefficient[~valid_gene] = np.nan
    se[~valid_gene] = np.nan
    p_value[~valid_gene] = np.nan
    return pd.DataFrame(
        {
            "gene_symbol": expression.index,
            "standardized_disease_beta": coefficient,
            "standard_error": se,
            "t_statistic": statistic,
            "degrees_freedom": df,
            "p_value": p_value,
        }
    )


def covariate_sensitivity(
    cohorts: dict[str, Cohort], pathways: dict[str, list[str]], universe: list[str]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    specifications: list[tuple[str, pd.DataFrame, list[str]]] = []

    for cohort_name, metadata_path, covariates in (
        ("GSE96804", ROOT / "data_processed" / "GSE96804" / "sample_annotation_clean.csv", ["sex"]),
        (
            "GSE166239",
            ROOT / "data_processed" / "GSE166239" / "sample_annotation_clean.csv",
            ["age", "sex"],
        ),
    ):
        cohort = cohorts[cohort_name]
        metadata = pd.read_csv(metadata_path).set_index("sample_id")
        samples = cohort.case_ids + cohort.control_ids
        cov = metadata.loc[samples, covariates].copy()
        cov["disease"] = [1.0 if sample in set(cohort.case_ids) else 0.0 for sample in samples]
        cov["sex"] = cov["sex"].astype(str).str.lower().map({"male": 1.0, "female": 0.0})
        if "age" in cov:
            cov["age"] = pd.to_numeric(cov["age"], errors="coerce")
            cov["age"] = (cov["age"] - cov["age"].mean()) / cov["age"].std(ddof=1)
        cov.insert(0, "intercept", 1.0)
        specifications.append((cohort_name, cov[["intercept", "disease", *covariates]], covariates))

    # GSE163603 count-table columns can be linked deterministically to the GEO
    # sample titles, which carry sex. Disease is partly confounded with sex.
    cohort = cohorts["GSE163603_interstitium"]
    metadata = pd.read_csv(ROOT / "data_raw" / "geo_soft" / "GSE163603" / "samples.csv")
    interstitium = metadata.loc[metadata["title"].str.contains("Interstitium", na=False)].copy()
    interstitium["sex"] = interstitium["characteristics"].str.extract(
        r"gender:\s*(Male|Female)", expand=False
    ).str.lower().map({"male": 1.0, "female": 0.0})
    sex_by_title = interstitium.set_index("title")["sex"]
    samples = cohort.case_ids + cohort.control_ids
    cov = pd.DataFrame(index=samples)
    cov["intercept"] = 1.0
    cov["disease"] = [1.0 if sample in set(cohort.case_ids) else 0.0 for sample in samples]
    cov["sex"] = cov.index.to_series().map(sex_by_title)
    specifications.append(("GSE163603_interstitium", cov, ["sex"]))

    gene_frames = []
    pathway_rows = []
    for cohort_name, adjusted_design, covariates in specifications:
        print(f"covariate sensitivity: starting {cohort_name}", flush=True)
        cohort = cohorts[cohort_name]
        if adjusted_design.isna().any().any():
            raise AssertionError(f"Missing covariate values in {cohort_name}")
        genes = sorted(set(universe) & set(cohort.expression.index))
        expression = cohort.expression.loc[genes]
        unadjusted_design = adjusted_design[["intercept", "disease"]]
        unadjusted = ols_disease_coefficient(expression, unadjusted_design)
        print(f"covariate sensitivity: unadjusted {cohort_name}", flush=True)
        adjusted = ols_disease_coefficient(expression, adjusted_design)
        print(f"covariate sensitivity: adjusted {cohort_name}", flush=True)
        merged = unadjusted.merge(adjusted, on="gene_symbol", suffixes=("_unadjusted", "_adjusted"))
        merged.insert(0, "cohort", cohort_name)
        merged.insert(1, "compartment", cohort.compartment)
        merged.insert(2, "adjustment", "+".join(covariates))
        merged["direction_changed"] = (
            np.sign(merged["standardized_disease_beta_unadjusted"])
            != np.sign(merged["standardized_disease_beta_adjusted"])
        )
        merged["adjusted_p_for_family"] = merged["p_value_adjusted"].fillna(1.0)
        merged["adjusted_fdr_bh_canonical_union"] = bh(merged["adjusted_p_for_family"])
        gene_frames.append(merged)
        by_gene = merged.set_index("gene_symbol")
        for sid, members in pathways.items():
            mapped = [gene for gene in members if gene in by_gene.index]
            before = float(by_gene.loc[mapped, "standardized_disease_beta_unadjusted"].mean())
            after = float(by_gene.loc[mapped, "standardized_disease_beta_adjusted"].mean())
            pathway_rows.append(
                {
                    "cohort": cohort_name,
                    "compartment": cohort.compartment,
                    "adjustment": "+".join(covariates),
                    "reactome_id": sid,
                    "reactome_name": next(name for name, value in CANONICAL_PATHWAYS.items() if value == sid),
                    "genes_mapped": len(mapped),
                    "mean_standardized_beta_unadjusted": before,
                    "mean_standardized_beta_adjusted": after,
                    "direction_changed": np.sign(before) != np.sign(after),
                    "limitation": (
                        "Measured-covariate OLS sensitivity only; residual confounding and disease-sex "
                        "partial confounding in GSE163603 remain possible."
                    ),
                }
            )
    return pd.concat(gene_frames, ignore_index=True), pd.DataFrame(pathway_rows)


def leave_one_sample_out(
    cohorts: dict[str, Cohort], pathways: dict[str, list[str]]
) -> pd.DataFrame:
    rows = []
    for cohort in cohorts.values():
        full_samples = cohort.case_ids + cohort.control_ids
        if min(len(cohort.case_ids), len(cohort.control_ids)) < 3:
            for sid in pathways:
                rows.append(
                    {
                        "cohort": cohort.name,
                        "reactome_id": sid,
                        "reactome_name": next(name for name, value in CANONICAL_PATHWAYS.items() if value == sid),
                        "removed_sample": "NOT_EVALUABLE",
                        "removed_group": "NOT_EVALUABLE",
                        "full_mean_hedges_g": np.nan,
                        "leave_one_sample_out_mean_hedges_g": np.nan,
                        "direction_changed": np.nan,
                        "reason": "At least three samples per group are required before removal.",
                    }
                )
            continue
        genes = sorted({gene for members in pathways.values() for gene in members if gene in cohort.expression.index})
        full_effect = hedges_effects(cohort.expression.loc[genes], cohort.case_ids, cohort.control_ids).set_index(
            "gene_symbol"
        )["hedges_g"]
        for removed in full_samples:
            case = [sample for sample in cohort.case_ids if sample != removed]
            control = [sample for sample in cohort.control_ids if sample != removed]
            loo_effect = hedges_effects(cohort.expression.loc[genes], case, control).set_index("gene_symbol")[
                "hedges_g"
            ]
            for sid, members in pathways.items():
                mapped = [gene for gene in members if gene in full_effect.index]
                full_value = float(full_effect.loc[mapped].mean())
                loo_value = float(loo_effect.loc[mapped].mean())
                rows.append(
                    {
                        "cohort": cohort.name,
                        "reactome_id": sid,
                        "reactome_name": next(name for name, value in CANONICAL_PATHWAYS.items() if value == sid),
                        "removed_sample": removed,
                        "removed_group": "DKD" if removed in set(cohort.case_ids) else "CONTROL",
                        "full_mean_hedges_g": full_value,
                        "leave_one_sample_out_mean_hedges_g": loo_value,
                        "absolute_change": abs(loo_value - full_value),
                        "direction_changed": np.sign(full_value) != np.sign(loo_value),
                        "reason": "",
                    }
                )
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    DOC.mkdir(parents=True, exist_ok=True)
    m8 = load_module("m19_m8_helpers", M8_PATH)
    m18 = load_module("m19_m18_meta", M18_PATH)
    m9_111154 = load_module("m19_m9_gse111154", M9_GSE111154_PATH)
    pathways = load_canonical_pathways()
    universe = sorted({gene for genes in pathways.values() for gene in genes})
    cohorts = prepare_cohorts(m8, m9_111154)

    characteristics = []
    effect_frames = []
    pathway_frames = []
    loo_frames = []
    for cohort in cohorts.values():
        missing = set(cohort.case_ids + cohort.control_ids) - set(cohort.expression.columns)
        if missing:
            raise AssertionError(f"{cohort.name} missing samples: {sorted(missing)}")
        characteristics.append(
            {
                "cohort": cohort.name,
                "source_study": cohort.source_study,
                "compartment": cohort.compartment,
                "platform": cohort.platform,
                "n_case": len(cohort.case_ids),
                "n_control": len(cohort.control_ids),
                "control_definition": cohort.control_definition,
                "processing": cohort.processing,
                "covariates_available": cohort.covariates_available,
                "analysis_tier": (
                    "primary glomerular meta-analysis"
                    if cohort.name in set(PRIMARY_GLOMERULAR_COHORTS)
                    else "small-sample glomerular sensitivity"
                    if cohort.name == "GSE1009_donor_averaged"
                    else "compartment-specific contextual analysis"
                ),
            }
        )
        available = [gene for gene in universe if gene in cohort.expression.index]
        effects = hedges_effects(cohort.expression.loc[available], cohort.case_ids, cohort.control_ids)
        complete = pd.DataFrame({"gene_symbol": universe}).merge(effects, on="gene_symbol", how="left")
        complete.insert(0, "cohort", cohort.name)
        complete.insert(1, "source_study", cohort.source_study)
        complete.insert(2, "compartment", cohort.compartment)
        complete["welch_p_for_canonical_union"] = complete["welch_p_value"].fillna(1.0)
        complete["welch_fdr_bh_canonical_union"] = bh(complete["welch_p_for_canonical_union"])
        effect_frames.append(complete)
        pathway, loo = pathway_permutation(cohort, pathways)
        pathway_frames.append(pathway)
        loo_frames.append(loo)
        print(f"completed cohort {cohort.name}", flush=True)

    characteristics_frame = pd.DataFrame(characteristics)
    effects_all = pd.concat(effect_frames, ignore_index=True)
    pathways_all = pd.concat(pathway_frames, ignore_index=True)
    loo_all = pd.concat(loo_frames, ignore_index=True)
    tier_by_cohort = characteristics_frame.set_index("cohort")["analysis_tier"]
    pathways_all["analysis_tier"] = pathways_all["cohort"].map(tier_by_cohort)
    effects_all["analysis_tier"] = effects_all["cohort"].map(tier_by_cohort)

    glomerular = meta_table(
        m18, effects_all, PRIMARY_GLOMERULAR_COHORTS, universe,
        "primary glomerular three-source-study synthesis",
    )
    glomerular_small_n = meta_table(
        m18,
        effects_all,
        ["GSE96804", "GSE30528", "GSE104948_H7", "GSE1009_donor_averaged"],
        universe,
        "small-sample glomerular sensitivity including donor-averaged GSE1009",
    )
    tubulo = meta_table(
        m18, effects_all, ["GSE30529", "GSE104954_H7"], universe,
        "descriptive tubulointerstitial two-source-study synthesis",
    )
    whole = meta_table(
        m18,
        effects_all,
        ["GSE142025_advanced", "GSE166239", "GSE111154", "GSE199838"],
        universe,
        "descriptive whole/cortical-kidney sensitivity including small studies",
    )
    print("completed compartment-specific meta-analyses", flush=True)

    replication = (
        pathways_all.groupby(["compartment", "reactome_id", "reactome_name"], dropna=False)
        .agg(
            n_source_studies=("source_study", "nunique"),
            n_cohorts=("cohort", "nunique"),
            n_positive=("observed_mean_hedges_g", lambda x: int((x > 0).sum())),
            n_negative=("observed_mean_hedges_g", lambda x: int((x < 0).sum())),
            n_maxT_fwer_lt_0_05=("maxT_fwer_p", lambda x: int((x < 0.05).sum())),
            min_maxT_fwer_p=("maxT_fwer_p", "min"),
        )
        .reset_index()
    )
    replication["interpretation"] = (
        "All-tier study-wise association summary. Primary replication calls are restricted to the "
        "prespecified three-source glomerular analysis and are reported separately."
    )

    primary_pathways = pathways_all.loc[pathways_all["cohort"].isin(PRIMARY_GLOMERULAR_COHORTS)].copy()
    primary_replication_rows = []
    for (sid, name), group in primary_pathways.groupby(["reactome_id", "reactome_name"], sort=False):
        significant = group.loc[group["maxT_fwer_p"] < 0.05]
        significant_concordant = len(significant) >= 2 and significant["direction"].nunique() == 1
        primary_replication_rows.append(
            {
                "reactome_id": sid,
                "reactome_name": name,
                "n_primary_source_studies": group["source_study"].nunique(),
                "n_positive": int((group["observed_mean_hedges_g"] > 0).sum()),
                "n_negative": int((group["observed_mean_hedges_g"] < 0).sum()),
                "n_maxT_fwer_lt_0_05": len(significant),
                "significant_directions": ";".join(sorted(significant["direction"].unique())),
                "primary_replication_call": bool(significant_concordant),
                "replication_rule": (
                    "At least two of the three independent glomerular source studies have maxT FWER<0.05 "
                    "with the same net-effect direction."
                ),
                "study_effects": ";".join(
                    f"{row.cohort}:g={row.observed_mean_hedges_g:.6g},maxT={row.maxT_fwer_p:.6g}"
                    for row in group.itertuples()
                ),
            }
        )
    primary_replication = pd.DataFrame(primary_replication_rows)
    print("completed pathway replication summaries", flush=True)

    qc = sample_qc(cohorts)
    print("completed sample-correlation QC", flush=True)
    covariate_genes, covariate_pathways = covariate_sensitivity(cohorts, pathways, universe)
    print("completed measured-covariate sensitivities", flush=True)
    sample_loo = leave_one_sample_out(cohorts, pathways)
    print("completed leave-one-sample-out sensitivities", flush=True)
    characteristics_frame.to_csv(OUT / "cohort_characteristics.csv", index=False)
    effects_all.to_csv(OUT / "per_cohort_canonical_gene_effects.csv", index=False)
    glomerular.to_csv(OUT / "primary_glomerular_gene_meta.csv", index=False)
    glomerular_small_n.to_csv(OUT / "small_sample_glomerular_gene_meta_sensitivity.csv", index=False)
    tubulo.to_csv(OUT / "descriptive_tubulointerstitial_gene_meta.csv", index=False)
    whole.to_csv(OUT / "descriptive_whole_kidney_gene_meta.csv", index=False)
    pathways_all.to_csv(OUT / "canonical_pathway_permutation_results.csv", index=False)
    replication.to_csv(OUT / "canonical_pathway_replication_summary.csv", index=False)
    primary_replication.to_csv(OUT / "primary_glomerular_pathway_replication_summary.csv", index=False)
    loo_all.to_csv(OUT / "canonical_pathway_leave_one_gene_out.csv", index=False)
    sample_loo.to_csv(OUT / "canonical_pathway_leave_one_sample_out.csv", index=False)
    qc.to_csv(OUT / "sample_correlation_qc.csv", index=False)
    covariate_genes.to_csv(OUT / "covariate_adjustment_gene_sensitivity.csv", index=False)
    covariate_pathways.to_csv(OUT / "covariate_adjustment_pathway_sensitivity.csv", index=False)

    # Probe aggregation sensitivity for the two GSE30122 compartments.
    sensitivity_rows = []
    pathway_sensitivity_rows = []
    for accession in ("GSE30528", "GSE30529"):
        cohort = cohorts[accession]
        median_expr = load_soft_gene_matrix(m8, accession, aggregation="median")
        genes = sorted(set(universe) & set(cohort.expression.index) & set(median_expr.index))
        primary = hedges_effects(cohort.expression.loc[genes], cohort.case_ids, cohort.control_ids)
        median = hedges_effects(median_expr.loc[genes], cohort.case_ids, cohort.control_ids)
        merged = primary[["gene_symbol", "hedges_g"]].merge(
            median[["gene_symbol", "hedges_g"]], on="gene_symbol", suffixes=("_highest_mean", "_median_probe")
        )
        merged.insert(0, "cohort", accession)
        merged["direction_changed"] = np.sign(merged["hedges_g_highest_mean"]) != np.sign(merged["hedges_g_median_probe"])
        sensitivity_rows.append(merged)
        by_gene = merged.set_index("gene_symbol")
        pearson = float(merged[["hedges_g_highest_mean", "hedges_g_median_probe"]].corr(method="pearson").iloc[0, 1])
        spearman = float(merged[["hedges_g_highest_mean", "hedges_g_median_probe"]].corr(method="spearman").iloc[0, 1])
        for sid, members in pathways.items():
            mapped = [gene for gene in members if gene in by_gene.index]
            high = float(by_gene.loc[mapped, "hedges_g_highest_mean"].mean())
            median_value = float(by_gene.loc[mapped, "hedges_g_median_probe"].mean())
            pathway_sensitivity_rows.append(
                {
                    "cohort": accession,
                    "reactome_id": sid,
                    "reactome_name": next(name for name, value in CANONICAL_PATHWAYS.items() if value == sid),
                    "genes_mapped": len(mapped),
                    "mean_hedges_g_highest_mean_probe": high,
                    "mean_hedges_g_median_probe": median_value,
                    "direction_changed": np.sign(high) != np.sign(median_value),
                    "gene_effect_pearson_across_canonical_union": pearson,
                    "gene_effect_spearman_across_canonical_union": spearman,
                }
            )
    pd.concat(sensitivity_rows, ignore_index=True).to_csv(
        OUT / "probe_aggregation_sensitivity.csv", index=False
    )
    pd.DataFrame(pathway_sensitivity_rows).to_csv(
        OUT / "probe_aggregation_pathway_sensitivity.csv", index=False
    )

    summary = {
        "canonical_pathways": len(pathways),
        "canonical_union_genes": len(universe),
        "cohorts": len(cohorts),
        "independent_source_studies": characteristics_frame["source_study"].nunique(),
        "primary_glomerular_complete_genes": int(glomerular["k"].eq(3).sum()),
        "primary_glomerular_fdr_lt_0_05": int(glomerular["fdr_lt_0_05"].sum()),
        "pathway_rows": len(pathways_all),
        "study_pathways_maxT_fwer_lt_0_05": int(pathways_all["maxT_fwer_p"].lt(0.05).sum()),
        "primary_glomerular_replicated_pathways": int(primary_replication["primary_replication_call"].sum()),
        "primary_glomerular_replicated_pathway_names": primary_replication.loc[
            primary_replication["primary_replication_call"], "reactome_name"
        ].tolist(),
        "sample_correlation_outliers": int(qc["correlation_outlier_flag"].sum()),
        "covariate_sensitivity_cohorts": sorted(covariate_genes["cohort"].unique().tolist()),
        "leave_one_sample_direction_changes": int(sample_loo["direction_changed"].fillna(False).astype(bool).sum()),
        "interpretive_boundary": (
            "No cross-compartment pooled effect and no independent pathway replication claim without "
            "two same-compartment source studies passing maxT FWER with concordant direction."
        ),
    }
    (OUT / "m19_analysis_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
