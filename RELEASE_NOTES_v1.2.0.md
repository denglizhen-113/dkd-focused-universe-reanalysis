# v1.2.0 — M20 reviewer-response release

Date: 21 August 2026

This release implements the complete internal Scientific Reports-style major-revision audit.

## Statistical changes

- Primary pathway estimands now use the fixed intersection of genes measurable across all three primary glomerular cohorts.
- Pathway family-wise inference now uses two-sided studentized maxT.
- GSE96804 uses 100,000 sex-restricted Monte Carlo allocations with plus-one probabilities and simulation-precision intervals.
- GSE30528 and GSE104948 H7 use exhaustive enumeration (497,420 and 480,700 allocations).
- Every pathway/source effect has a 5,000-resample bootstrap 95% interval.
- Gene synthesis now includes leave-one-source-out results for all 783 family members.

## Result changes

- Four pathways meet the operational two-of-three-source rule: complement cascade, vascular-wall interaction, chemokine-receptor binding and extracellular-matrix organization.
- Complement is weak in GSE96804 after common-gene restriction and meets the rule through the other two sources.
- Coagulation remains directionally mixed and meets the rule in zero sources.
- No gene meets BH FDR<0.05 across the 783-gene family.

## Reporting and reproducibility

- Full PRISMA 2020 item/sub-item checklist and an adapted dataset-record flow diagram.
- Formal source-level control, confounding, overlap and ethics/consent provenance table.
- Transparent AI-assisted second-pass screening audit that is not misrepresented as independent human dual screening.
- Formal citation of all 11 GEO Series and direct comparison with related DKD integrations.
- Figures 1–4 revised to journal scale; Figure 4 now displays uncertainty and Figure 3 makes selection/multiplicity status explicit.
- `Source_Code_M20.zip` contains frozen inputs, requirements, README and regression tests. A clean extraction reproduced all seven primary outputs byte-for-byte.

## Residual author actions

The author must personally confirm local ethics requirements, unresolved original-source ethics fields, any independent-human screening requirement, reviewer/exclusion suggestions and editorial-board relationships. No DOI is claimed.
