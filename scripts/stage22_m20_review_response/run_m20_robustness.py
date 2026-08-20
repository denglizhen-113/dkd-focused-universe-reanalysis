#!/usr/bin/env python3
"""Design-aware M20 robustness analysis for the Scientific Reports revision.

The primary pathway estimand is the mean Hedges' g over genes measurable in
all three primary glomerular cohorts.  Family-wise inference uses a
studentized maxT statistic.  GSE96804 permutations are restricted within sex;
the two smaller cohorts have no usable person-level design factors and are
therefore enumerated without restriction, with that limitation retained in
the output.  The script also reports Monte Carlo precision, paired
cohort-specific mapped-set sensitivities, pathway bootstrap intervals, and
leave-one-source-out gene syntheses.

Use ``--prepare-inputs`` once in the full repository to create the frozen,
minimal primary input bundle.  Thereafter the analysis is self-contained in
that bundle and does not import any earlier M8/M9/M18 project helper.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import t as t_dist


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data_processed" / "m20_primary_reproduction"
OUT = ROOT / "tables" / "stage22_m20_review_response"
SEED = 20260821
MONTE_CARLO_ALLOCATIONS = 100_000
BOOTSTRAP_REPLICATES = 5_000
BATCH_SIZE = 500
COHORTS = ("GSE96804", "GSE30528", "GSE104948_H7")


@dataclass(frozen=True)
class CohortInput:
    name: str
    expression: pd.DataFrame
    metadata: pd.DataFrame
    restriction: str
    exchangeability_scope: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def prepare_inputs() -> None:
    """Freeze the three primary matrices and design metadata for clean replay."""
    m19_path = ROOT / "scripts" / "stage21_m19_scientific_reports_revision" / "run_m19_compartment_analysis.py"
    m8_path = ROOT / "scripts" / "stage21_m8_validation" / "run_m8_analysis.py"
    m9_path = ROOT / "scripts" / "stage21_m9_gse111154" / "run_gse111154_rerun.py"
    m19 = load_module("m20_m19_export", m19_path)
    m8 = load_module("m20_m8_export", m8_path)
    m9 = load_module("m20_m9_export", m9_path)
    pathways = m19.load_canonical_pathways()
    universe = sorted({gene for genes in pathways.values() for gene in genes})
    cohorts = m19.prepare_cohorts(m8, m9)
    INPUT.mkdir(parents=True, exist_ok=True)

    frozen_gmt = INPUT / "canonical_pathways.gmt"
    with frozen_gmt.open("w", encoding="utf-8", newline="\n") as handle:
        for sid, genes in pathways.items():
            name = next(name for name, value in m19.CANONICAL_PATHWAYS.items() if value == sid)
            handle.write("\t".join([sid, name, *genes]) + "\n")

    source_rows = []
    sex_meta = pd.read_csv(ROOT / "data_processed" / "GSE96804" / "sample_annotation_clean.csv")
    sex_by_sample = sex_meta.set_index("sample_id")["sex"].astype(str).str.lower()
    for name in COHORTS:
        cohort = cohorts[name]
        samples = cohort.case_ids + cohort.control_ids
        expression = cohort.expression.loc[cohort.expression.index.intersection(universe), samples].copy()
        expression.index.name = "gene_symbol"
        expression_path = INPUT / f"{name}__canonical_expression.csv.gz"
        expression.to_csv(expression_path, compression="gzip")
        metadata = pd.DataFrame(
            {
                "sample_id": samples,
                "group": ["DKD" if sample in set(cohort.case_ids) else "CONTROL" for sample in samples],
                "stratum": [
                    sex_by_sample.get(sample, "unknown") if name == "GSE96804" else "single_stratum"
                    for sample in samples
                ],
            }
        )
        metadata_path = INPUT / f"{name}__design.csv"
        metadata.to_csv(metadata_path, index=False)
        source_rows.extend(
            [
                {"file": expression_path.name, "bytes": expression_path.stat().st_size, "sha256": sha256(expression_path)},
                {"file": metadata_path.name, "bytes": metadata_path.stat().st_size, "sha256": sha256(metadata_path)},
            ]
        )
    source_rows.append({"file": frozen_gmt.name, "bytes": frozen_gmt.stat().st_size, "sha256": sha256(frozen_gmt)})
    pd.DataFrame(source_rows).to_csv(INPUT / "input_manifest.csv", index=False)


def load_pathways() -> tuple[dict[str, list[str]], dict[str, str]]:
    pathways: dict[str, list[str]] = {}
    names: dict[str, str] = {}
    path = INPUT / "canonical_pathways.gmt"
    for line in path.read_text(encoding="utf-8").splitlines():
        sid, name, *genes = line.split("\t")
        pathways[sid] = genes
        names[sid] = name
    if len(pathways) != 7:
        raise AssertionError(f"Expected seven pathways, found {len(pathways)}")
    return pathways, names


def load_cohorts() -> dict[str, CohortInput]:
    result = {}
    for name in COHORTS:
        expr = pd.read_csv(INPUT / f"{name}__canonical_expression.csv.gz", index_col=0)
        expr.index = expr.index.astype(str).str.upper()
        design = pd.read_csv(INPUT / f"{name}__design.csv")
        if name == "GSE96804":
            restriction = "case count preserved within archived sex strata"
            scope = "Sex-restricted exchangeability; age, batch, renal function, medication, ancestry, and control-origin confounding remain untestable."
        elif name == "GSE104948_H7":
            restriction = "H7 batch fixed before analysis; unrestricted within H7"
            scope = "Conditional on within-H7 label exchangeability; individual demographics and clinical covariates were unavailable."
        else:
            restriction = "unrestricted within cohort; no person-level design factor available"
            scope = "Conditional on cohort-wide label exchangeability; clinical covariates and control-selection confounding were unavailable."
        result[name] = CohortInput(name, expr, design, restriction, scope)
    return result


def hedges_g(values: np.ndarray, case_mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    case = values[:, case_mask]
    control = values[:, ~case_mask]
    n1, n0 = case.shape[1], control.shape[1]
    mean1, mean0 = case.mean(axis=1), control.mean(axis=1)
    sd1, sd0 = case.std(axis=1, ddof=1), control.std(axis=1, ddof=1)
    pooled = np.sqrt(((n1 - 1) * sd1**2 + (n0 - 1) * sd0**2) / (n1 + n0 - 2))
    correction = 1.0 - 3.0 / (4.0 * (n1 + n0) - 9.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        g = correction * (mean1 - mean0) / pooled
    variance = (n1 + n0) / (n1 * n0) + g**2 / (2.0 * (n1 + n0 - 2.0))
    g[~np.isfinite(g)] = np.nan
    variance[~np.isfinite(g)] = np.nan
    return g, variance


def batched_hedges_g(values: np.ndarray, masks: np.ndarray) -> np.ndarray:
    """Return gene x allocation Hedges' g using matrix sufficient statistics."""
    weights = masks.astype(float)
    n1 = weights.sum(axis=1)
    n0 = values.shape[1] - n1
    sums1 = values @ weights.T
    sums_sq1 = np.square(values) @ weights.T
    total = values.sum(axis=1)[:, None]
    total_sq = np.square(values).sum(axis=1)[:, None]
    sums0 = total - sums1
    sums_sq0 = total_sq - sums_sq1
    means1 = sums1 / n1
    means0 = sums0 / n0
    vars1 = np.maximum((sums_sq1 - sums1**2 / n1) / (n1 - 1), 0.0)
    vars0 = np.maximum((sums_sq0 - sums0**2 / n0) / (n0 - 1), 0.0)
    pooled = np.sqrt(((n1 - 1) * vars1 + (n0 - 1) * vars0) / (n1 + n0 - 2))
    correction = 1.0 - 3.0 / (4.0 * (n1 + n0) - 9.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        g = correction * (means1 - means0) / pooled
    return np.where(np.isfinite(g), g, 0.0)


def exact_mask_batches(n: int, n_case: int) -> Iterator[np.ndarray]:
    batch: list[tuple[int, ...]] = []
    for combination in itertools.combinations(range(n), n_case):
        batch.append(combination)
        if len(batch) == BATCH_SIZE:
            masks = np.zeros((len(batch), n), dtype=bool)
            for row, indexes in enumerate(batch):
                masks[row, list(indexes)] = True
            yield masks
            batch = []
    if batch:
        masks = np.zeros((len(batch), n), dtype=bool)
        for row, indexes in enumerate(batch):
            masks[row, list(indexes)] = True
        yield masks


def monte_carlo_mask_batches(metadata: pd.DataFrame, rng: np.random.Generator) -> Iterator[np.ndarray]:
    n = len(metadata)
    remaining = MONTE_CARLO_ALLOCATIONS
    strata = []
    observed = metadata["group"].eq("DKD").to_numpy()
    for _, index in metadata.groupby("stratum", sort=True).groups.items():
        positions = metadata.index.get_indexer(index)
        strata.append((positions, int(observed[positions].sum())))
    while remaining:
        size = min(BATCH_SIZE, remaining)
        masks = np.zeros((size, n), dtype=bool)
        for positions, n_case in strata:
            scores = rng.random((size, len(positions)))
            chosen = np.argpartition(scores, n_case - 1, axis=1)[:, :n_case]
            rows = np.repeat(np.arange(size), n_case)
            cols = positions[chosen.reshape(-1)]
            masks[rows, cols] = True
        yield masks
        remaining -= size


def wilson_interval(exceed: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    p = exceed / total
    denominator = 1.0 + z * z / total
    center = (p + z * z / (2.0 * total)) / denominator
    half = z * math.sqrt(p * (1.0 - p) / total + z * z / (4.0 * total * total)) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def pathway_weights(genes: list[str], members: dict[str, list[str]]) -> np.ndarray:
    position = {gene: i for i, gene in enumerate(genes)}
    weights = np.zeros((len(members), len(genes)), dtype=float)
    for row, pathway_genes in enumerate(members.values()):
        indexes = [position[g] for g in pathway_genes if g in position]
        weights[row, indexes] = 1.0 / len(indexes)
    return weights


def bootstrap_ci(
    values: np.ndarray,
    case_mask: np.ndarray,
    weights: np.ndarray,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    case_values = values[:, case_mask]
    control_values = values[:, ~case_mask]
    n1, n0 = case_values.shape[1], control_values.shape[1]
    collected = []
    remaining = BOOTSTRAP_REPLICATES
    correction = 1.0 - 3.0 / (4.0 * (n1 + n0) - 9.0)
    while remaining:
        size = min(BATCH_SIZE, remaining)
        case_draw = rng.integers(0, n1, size=(size, n1))
        control_draw = rng.integers(0, n0, size=(size, n0))
        case_counts = np.zeros((size, n1), dtype=float)
        control_counts = np.zeros((size, n0), dtype=float)
        np.add.at(case_counts, (np.repeat(np.arange(size), n1), case_draw.ravel()), 1.0)
        np.add.at(control_counts, (np.repeat(np.arange(size), n0), control_draw.ravel()), 1.0)
        sum1 = case_values @ case_counts.T
        sum0 = control_values @ control_counts.T
        ss1 = np.square(case_values) @ case_counts.T
        ss0 = np.square(control_values) @ control_counts.T
        mean1, mean0 = sum1 / n1, sum0 / n0
        var1 = np.maximum((ss1 - sum1**2 / n1) / (n1 - 1), 0.0)
        var0 = np.maximum((ss0 - sum0**2 / n0) / (n0 - 1), 0.0)
        pooled = np.sqrt(((n1 - 1) * var1 + (n0 - 1) * var0) / (n1 + n0 - 2))
        with np.errstate(divide="ignore", invalid="ignore"):
            g = correction * (mean1 - mean0) / pooled
        g = np.where(np.isfinite(g), g, 0.0)
        collected.append((weights @ g).T)
        remaining -= size
    stats = np.concatenate(collected, axis=0)
    return np.quantile(stats, 0.025, axis=0), np.quantile(stats, 0.975, axis=0)


def permutation_analysis(
    cohort: CohortInput,
    pathways: dict[str, list[str]],
    names: dict[str, str],
    common: dict[str, list[str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = cohort.metadata.reset_index(drop=True)
    samples = metadata["sample_id"].tolist()
    genes = cohort.expression.index.tolist()
    values = cohort.expression.loc[genes, samples].to_numpy(float)
    if not np.isfinite(values).all():
        raise AssertionError(f"Non-finite values in frozen primary matrix: {cohort.name}")
    observed_mask = metadata["group"].eq("DKD").to_numpy()
    # Use the same sufficient-statistic implementation for observed and
    # permuted allocations so exact-enumeration equality is not lost to
    # harmless floating-point differences.
    observed_g = batched_hedges_g(values, observed_mask[None, :])[:, 0]
    mapped_members = {sid: [g for g in members if g in set(genes)] for sid, members in pathways.items()}
    common_members = {sid: [g for g in members if g in set(genes)] for sid, members in common.items()}
    weights_by_estimand = {
        "common_three_cohort_intersection": pathway_weights(genes, common_members),
        "cohort_specific_mapped_set_sensitivity": pathway_weights(genes, mapped_members),
    }
    observed = {key: weights @ np.nan_to_num(observed_g, nan=0.0) for key, weights in weights_by_estimand.items()}

    n = len(samples)
    n_case = int(observed_mask.sum())
    total_unrestricted = math.comb(n, n_case)
    exact = cohort.name != "GSE96804" and total_unrestricted <= 600_000
    if exact:
        mask_batches = exact_mask_batches(n, n_case)
        method = "exact enumeration"
        denominator = total_unrestricted
    else:
        mask_batches = monte_carlo_mask_batches(metadata, np.random.default_rng(SEED + COHORTS.index(cohort.name) + 1))
        method = "restricted Monte Carlo" if cohort.name == "GSE96804" else "Monte Carlo"
        denominator = MONTE_CARLO_ALLOCATIONS

    permutation_stats = {key: [] for key in weights_by_estimand}
    for masks in mask_batches:
        gene_g = batched_hedges_g(values, masks)
        for key, weights in weights_by_estimand.items():
            permutation_stats[key].append((weights @ gene_g).T)
    permutation_stats = {key: np.concatenate(parts, axis=0) for key, parts in permutation_stats.items()}
    if any(len(stats) != denominator for stats in permutation_stats.values()):
        raise AssertionError("Permutation allocation count mismatch")

    ci_low, ci_high = bootstrap_ci(
        values,
        observed_mask,
        weights_by_estimand["common_three_cohort_intersection"],
        np.random.default_rng(SEED + 100 + COHORTS.index(cohort.name)),
    )

    rows = []
    precision_rows = []
    for estimand, stats in permutation_stats.items():
        null_mean = stats.mean(axis=0)
        null_sd = stats.std(axis=0, ddof=1)
        obs = observed[estimand]
        obs_studentized = (obs - null_mean) / null_sd
        perm_studentized = (stats - null_mean) / null_sd
        max_studentized = np.max(np.abs(perm_studentized), axis=1)
        max_raw = np.max(np.abs(stats), axis=1)
        members_for_estimand = common_members if estimand.startswith("common") else mapped_members
        for index, sid in enumerate(pathways):
            tolerance = 1e-12
            exceed_single = int(np.count_nonzero(np.abs(perm_studentized[:, index]) >= abs(obs_studentized[index]) - tolerance))
            exceed_max = int(np.count_nonzero(max_studentized >= abs(obs_studentized[index]) - tolerance))
            exceed_raw = int(np.count_nonzero(max_raw >= abs(obs[index]) - tolerance))
            if exact:
                # The observed allocation is one member of the exact reference
                # set, hence a valid exact tail count cannot be zero.
                exceed_single = max(exceed_single, 1)
                exceed_max = max(exceed_max, 1)
                exceed_raw = max(exceed_raw, 1)
                p_single = exceed_single / denominator
                p_max = exceed_max / denominator
                p_raw = exceed_raw / denominator
                mc_se = mc_low = mc_high = np.nan
            else:
                p_single = (exceed_single + 1) / (denominator + 1)
                p_max = (exceed_max + 1) / (denominator + 1)
                p_raw = (exceed_raw + 1) / (denominator + 1)
                mc_se = math.sqrt(p_max * (1.0 - p_max) / (denominator + 1))
                mc_low, mc_high = wilson_interval(exceed_max, denominator)
            rows.append(
                {
                    "cohort": cohort.name,
                    "reactome_id": sid,
                    "reactome_name": names[sid],
                    "estimand": estimand,
                    "genes_in_reference_set": len(pathways[sid]),
                    "genes_used": len(members_for_estimand[sid]),
                    "observed_mean_hedges_g": obs[index],
                    "bootstrap_95_ci_low": ci_low[index] if estimand.startswith("common") else np.nan,
                    "bootstrap_95_ci_high": ci_high[index] if estimand.startswith("common") else np.nan,
                    "permutation_null_mean": null_mean[index],
                    "permutation_null_sd": null_sd[index],
                    "observed_studentized_statistic": obs_studentized[index],
                    "two_sided_studentized_permutation_p": p_single,
                    "studentized_maxT_fwer_p": p_max,
                    "raw_mean_maxT_fwer_p_sensitivity": p_raw,
                    "permutation_method": method,
                    "permutation_allocations": denominator,
                    "design_restriction": cohort.restriction,
                    "exchangeability_scope": cohort.exchangeability_scope,
                    "maxT_alpha_two_sided": 0.05,
                }
            )
            precision_rows.append(
                {
                    "cohort": cohort.name,
                    "reactome_id": sid,
                    "estimand": estimand,
                    "maxT_exceedances": exceed_max,
                    "allocations": denominator,
                    "reported_maxT_p": p_max,
                    "monte_carlo_se": mc_se,
                    "monte_carlo_wilson_95_low": mc_low,
                    "monte_carlo_wilson_95_high": mc_high,
                    "precision_note": "exact enumeration; no Monte Carlo error" if exact else "Wilson interval describes simulation uncertainty, not biological uncertainty",
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(precision_rows)


def restricted_loglik(tau2: float, y: np.ndarray, v: np.ndarray) -> float:
    total_variance = v + tau2
    weights = 1.0 / total_variance
    mean = float(np.sum(weights * y) / np.sum(weights))
    residual = float(np.sum(weights * np.square(y - mean)))
    return -0.5 * (float(np.sum(np.log(total_variance))) + math.log(float(np.sum(weights))) + residual)


def reml_hk(y: np.ndarray, v: np.ndarray) -> dict[str, float]:
    valid = np.isfinite(y) & np.isfinite(v) & (v > 0)
    y, v = y[valid], v[valid]
    k = len(y)
    if k < 2:
        return {key: np.nan for key in ("pooled_effect", "tau2_reml", "ci_low_modified_hk", "ci_high_modified_hk", "p_modified_hk", "i2_percent") } | {"k": k}
    upper = max(10.0, float(np.var(y, ddof=1) * 100.0 + np.max(v)))
    optimized = minimize_scalar(lambda x: -restricted_loglik(float(x), y, v), bounds=(0.0, upper), method="bounded", options={"xatol": 1e-12})
    tau2 = float(max(0.0, optimized.x))
    if tau2 < 1e-10 or restricted_loglik(0.0, y, v) >= restricted_loglik(tau2, y, v) - 1e-10:
        tau2 = 0.0
    weights = 1.0 / (v + tau2)
    pooled = float(np.sum(weights * y) / np.sum(weights))
    hk_scale = float(np.sum(weights * np.square(y - pooled)) / (k - 1))
    se = math.sqrt(max(hk_scale, 1.0) / float(np.sum(weights)))
    critical = float(t_dist.ppf(0.975, k - 1))
    p = float(2.0 * t_dist.sf(abs(pooled / se), k - 1)) if se > 0 else 1.0
    fixed_weights = 1.0 / v
    fixed_mean = float(np.sum(fixed_weights * y) / np.sum(fixed_weights))
    q = float(np.sum(fixed_weights * np.square(y - fixed_mean)))
    i2 = max(0.0, (q - (k - 1)) / q * 100.0) if q > 0 else 0.0
    return {"k": k, "pooled_effect": pooled, "tau2_reml": tau2, "ci_low_modified_hk": pooled - critical * se, "ci_high_modified_hk": pooled + critical * se, "p_modified_hk": p, "i2_percent": i2}


def gene_meta_and_leave_one_source(cohorts: dict[str, CohortInput], pathways: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    universe = sorted({gene for members in pathways.values() for gene in members})
    by_cohort = {}
    effect_rows = []
    for name, cohort in cohorts.items():
        metadata = cohort.metadata.reset_index(drop=True)
        samples = metadata["sample_id"].tolist()
        available = [gene for gene in universe if gene in cohort.expression.index]
        values = cohort.expression.loc[available, samples].to_numpy(float)
        g, variance = hedges_g(values, metadata["group"].eq("DKD").to_numpy())
        by_cohort[name] = {gene: (effect, var) for gene, effect, var in zip(available, g, variance)}
        effect_rows.extend({"cohort": name, "gene_symbol": gene, "hedges_g": effect, "variance": var} for gene, effect, var in zip(available, g, variance))
    loo_rows = []
    for omitted in ("NONE", *COHORTS):
        retained = [name for name in COHORTS if name != omitted]
        for gene in universe:
            pairs = [by_cohort[name].get(gene, (np.nan, np.nan)) for name in retained]
            y = np.array([pair[0] for pair in pairs], dtype=float)
            v = np.array([pair[1] for pair in pairs], dtype=float)
            result = reml_hk(y, v)
            loo_rows.append(
                {
                    "gene_symbol": gene,
                    "omitted_source": omitted,
                    "retained_sources": ";".join(retained),
                    **result,
                }
            )
    return pd.DataFrame(effect_rows), pd.DataFrame(loo_rows)


def coverage_table(cohorts: dict[str, CohortInput], pathways: dict[str, list[str]], names: dict[str, str]) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    common = {}
    rows = []
    gene_sets = {name: set(cohort.expression.index) for name, cohort in cohorts.items()}
    for sid, members in pathways.items():
        mapped = {name: set(members) & genes for name, genes in gene_sets.items()}
        intersection = set.intersection(*mapped.values())
        union = set.union(*mapped.values())
        common[sid] = sorted(intersection)
        rows.append(
            {
                "reactome_id": sid,
                "reactome_name": names[sid],
                "reference_members": len(members),
                **{f"{name}_mapped": len(mapped[name]) for name in COHORTS},
                "common_three_cohort_members": len(intersection),
                "mapped_set_jaccard_intersection_over_union": len(intersection) / len(union),
                "common_genes": ";".join(sorted(intersection)),
            }
        )
    return pd.DataFrame(rows), common


def replication_summary(pathway_results: pd.DataFrame) -> pd.DataFrame:
    primary = pathway_results.loc[pathway_results["estimand"].eq("common_three_cohort_intersection")]
    rows = []
    for (sid, name), group in primary.groupby(["reactome_id", "reactome_name"], sort=False):
        significant = group.loc[group["studentized_maxT_fwer_p"] < 0.05]
        directions = np.sign(significant["observed_mean_hedges_g"].to_numpy(float))
        met = len(significant) >= 2 and len(set(directions.tolist())) == 1
        rows.append(
            {
                "reactome_id": sid,
                "reactome_name": name,
                "positive_sources": int((group["observed_mean_hedges_g"] > 0).sum()),
                "negative_sources": int((group["observed_mean_hedges_g"] < 0).sum()),
                "sources_studentized_maxT_lt_0_05": len(significant),
                "operational_two_of_three_call": bool(met),
                "interpretation": "Criterion met under stated cohort-specific exchangeability and common-measurement assumptions" if met else "Criterion not met; this is not evidence of absence",
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepare-inputs", action="store_true")
    args = parser.parse_args()
    if args.prepare_inputs:
        prepare_inputs()
    OUT.mkdir(parents=True, exist_ok=True)
    pathways, names = load_pathways()
    cohorts = load_cohorts()
    coverage, common = coverage_table(cohorts, pathways, names)
    result_frames, precision_frames = [], []
    for name in COHORTS:
        print(f"M20 robust pathway inference: {name}", flush=True)
        result, precision = permutation_analysis(cohorts[name], pathways, names, common)
        result_frames.append(result)
        precision_frames.append(precision)
    pathway_results = pd.concat(result_frames, ignore_index=True)
    precision = pd.concat(precision_frames, ignore_index=True)
    effects, gene_loo = gene_meta_and_leave_one_source(cohorts, pathways)
    replication = replication_summary(pathway_results)

    coverage.to_csv(OUT / "primary_pathway_measurement_coverage.csv", index=False)
    pathway_results.to_csv(OUT / "primary_pathway_studentized_maxT_results.csv", index=False)
    precision.to_csv(OUT / "permutation_monte_carlo_precision.csv", index=False)
    replication.to_csv(OUT / "primary_pathway_operational_replication.csv", index=False)
    effects.to_csv(OUT / "frozen_primary_gene_effects.csv", index=False)
    gene_loo.to_csv(OUT / "gene_meta_leave_one_source_out.csv", index=False)
    summary = {
        "analysis": "M20 design-aware reviewer-response robustness analysis",
        "primary_estimand": "pathway mean Hedges' g over the three-cohort common measurable gene intersection",
        "multiplicity": "two-sided studentized maxT FWER across seven fixed Reactome pathways within each cohort",
        "permutations": {"GSE96804": "100,000 Monte Carlo allocations restricted within sex", "GSE30528": "all 497,420 allocations", "GSE104948_H7": "all 480,700 allocations within H7"},
        "bootstrap_replicates_per_cohort": BOOTSTRAP_REPLICATES,
        "replicated_pathways": replication.loc[replication["operational_two_of_three_call"], "reactome_name"].tolist(),
        "exchangeability_boundary": "Association tests remain conditional on archived-label exchangeability within the stated restriction; unavailable covariates cannot be repaired computationally.",
        "input_manifest_sha256": sha256(INPUT / "input_manifest.csv"),
    }
    (OUT / "m20_robustness_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)
    print("M20_ROBUSTNESS_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
