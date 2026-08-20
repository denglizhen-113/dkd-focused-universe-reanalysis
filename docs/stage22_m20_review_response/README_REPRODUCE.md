# M20 clean reproduction

This archive is self-contained for the reviewer-requested primary robustness analysis. It does not download mutable external resources and does not import project-external helper modules.

From the extracted archive root, run:

```text
python scripts/stage22_m20_review_response/run_m20_robustness.py
python scripts/stage22_m20_review_response/test_m20_reproducibility.py
```

The first command regenerates the common-measurement pathway results, restricted/exact studentized maxT inference, Monte Carlo precision table, bootstrap intervals and gene leave-one-source-out table under `tables/stage22_m20_review_response/`. The second command verifies frozen-input hashes, exact allocation counts, pathway coverage, non-zero exact P values, operational calls and the 783-gene leave-one-source-out family.

Tested environment: Python 3.12.12, NumPy 2.4.6, pandas 3.0.3 and SciPy 1.18.0. The plotting package was Matplotlib 3.10.9. No network access is required.

Expected terminal result:

```text
M20_ROBUSTNESS_COMPLETE
M20_NUMERICAL_REGRESSION=PASS
```

The fixed authoritative pathway input for this reproduction is `data_processed/m20_primary_reproduction/canonical_pathways.gmt`. Historical Reactome/GOA source-file hashes remain available for provenance, but rebuilding this frozen GMT from potentially mutable URLs is deliberately outside the clean-run dependency chain.
