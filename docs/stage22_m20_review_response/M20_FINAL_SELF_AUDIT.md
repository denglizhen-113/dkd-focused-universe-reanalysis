# M20 final self-audit

Date: 21 August 2026  
Decision: **Submission package technically ready; conditional author actions remain.**

## Automated checks

- `M20_NUMERICAL_REGRESSION=PASS`
- `M20_SUBMISSION_VALIDATION=PASS abstract_words=164 ready_files=11 prisma_rows=42`
- Frozen-input SHA-256 checks passed.
- Common pathway coverage matched the fixed expected counts: 51, 52, 125, 46, 255, 47 and 85 genes.
- Allocation counts matched: 100,000 restricted Monte Carlo (GSE96804), 497,420 exact (GSE30528) and 480,700 exact (GSE104948 H7).
- No exact or Monte Carlo P value was zero.
- Operational call set matched exactly: complement, vascular-wall interaction, chemokine-receptor binding and extracellular-matrix organization; coagulation remained 0/3.
- Leave-one-source-out output retained all 783 genes for the full and three omission strata.
- Abstract length is 164 words; title is below 20 words; 28 references are sequential; all 11 GEO records are formally cited.
- PRISMA checklist contains all 42 item/sub-item rows.
- Supplementary workbook contains README, data dictionary and S1–S27.
- Each main figure is a one-page vector PDF and was visually inspected after generation.
- Upload set contains exactly 11 files and every hash matches `final_upload_manifest.csv`.

## Clean-room reproduction

`Source_Code_M20.zip` was extracted to `work/m20_clean_room_v2`. The two documented commands were executed with no network access and no files imported from the parent project. The run printed:

```text
M20_ROBUSTNESS_COMPLETE
M20_NUMERICAL_REGRESSION=PASS
```

The SHA-256 hashes of all seven regenerated outputs were identical to the package-building workspace outputs:

| Output | SHA-256 |
| --- | --- |
| primary_pathway_measurement_coverage.csv | CC6F1B343222CA6D53D19BB0F8B1DE320CD59441E0D5B61521703FA67C3EA82A |
| primary_pathway_studentized_maxT_results.csv | 54DE16B1F3A11CE67496BB04B0D6722D438B20D9485F116B22E6556C18A21B59 |
| permutation_monte_carlo_precision.csv | 0E81D1052C7AB726F9418F25DD95E714CD14F117CD8887245750E7AF09DF0249 |
| primary_pathway_operational_replication.csv | 00D0E3DAF80EE1AC815DCCBA28E19A263B4C05831D6D0E21DD37B329B927D6B1 |
| frozen_primary_gene_effects.csv | 70A4C51A49245842C51DEEF36AD933EF696F01CDC4A6E8F80A2F1ED10A33F74F |
| gene_meta_leave_one_source_out.csv | D589DF295C191DBF6E0B409756EF84CEBB9B2C535C3869099D6178D41C393F73 |
| m20_robustness_summary.json | B4A92896380049448E2712F60BB6A87A8405C3CE1A0641DF87975BD411CE5AF3 |

## Visual review

- Figure 1: PRISMA-style dataset-record flow, explicit reporting unit and search date; counts reconcile 374−52=322 and 322−269=53, with 53−42=11.
- Figure 2: control/DKD sample-size mirror and source/compartment independence map; color palette is color-vision-deficiency compatible.
- Figure 3: descriptive selection is prominent in the title and panel annotation; 0/783 is visible.
- Figure 4: uncertainty intervals, source-specific shapes/colors, filled/open maxT status and operational count are readable without relying on color alone.
- Supplementary Figure S2: direction changes use an outlined X in addition to color.

## Scientific self-review

The analysis is materially stronger than M19, but it remains an observational reanalysis with only three independent primary sources. Restricted permutations do not correct unmeasured confounding; bootstrap intervals are conditional on archived cohorts; exact enumeration does not make the biological design randomized; source risk remains moderate-high; and the two-of-three rule is operational rather than a universal replication standard. These limitations are now explicit in the Abstract, Results, Methods and Discussion.

No evidence of fabricated primary data, image manipulation or undisclosed positive-result substitution was found. The change from three to four operational pathway calls is fully disclosed, including weakening of complement in GSE96804. The audit does not certify unresolved original-study ethics fields and does not substitute AI review for an independent human screener.

## Author-only actions before portal submission

- Confirm local institutional ethics/exemption requirements.
- Perform or arrange independent human screening confirmation if demanded.
- Enter reviewer suggestions/exclusions and editorial-board relationships personally.
- Verify all author/contact details and every source ethics statement against the full original publications.
- Optionally deposit the tagged release in a DOI-bearing repository and update Code availability.
