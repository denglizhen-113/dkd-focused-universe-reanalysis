# M19 post-beautification submission audit

**Audit date:** 21 August 2026  
**Target journal:** Scientific Reports  
**Package:** `submission_ready_scientific_reports_m19`  
**Overall decision:** PASS for technical upload and internal scientific consistency, conditional only on author-supplied portal declarations noted below.

## Main findings

- Figures 1–4 are single-page vector PDFs with embedded fonts and maximum width 182.9 mm.
- Figure 1 is a compact dataset-level PRISMA-style flow; Figure 2 separates sample size, source identity, compartment, and analysis role; Figure 3 separates selected-gene forest display from the complete-gene landscape; Figure 4 separates study-wise maxT evidence from the replication rule.
- All categorical colors are based on the color-vision-accessible Okabe–Ito palette. The Figure 4 signed-effect scale is symmetric about zero, and the complete figure remains interpretable in grayscale through direct values, asterisks, a threshold line, and text calls.
- No search count, cohort definition, numerical result, statistical decision, or scientific conclusion changed from v1.1.1.
- The manuscript renders as 7 pages, the Supplementary Information as 4 pages, and the cover letter and PRISMA checklist as one page each.
- The upload directory contains exactly 11 files and all final-manifest hashes match.
- Automated validation passes in both the working project and release checkout: `M19_SUBMISSION_VALIDATION=PASS abstract_words=176 ready_files=11`.

## Remaining disclosed limitations

The primary analysis has three independent glomerular sources; screening was performed by one reviewer; the pathway family was not prospectively preregistered; archived covariates and preprocessing are heterogeneous; and pathway transcriptional effects do not establish causality, protein activity, independent mechanisms, or cell of origin.

## Author-only confirmations

Before portal submission, the corresponding author must confirm any prior discussion with a Scientific Reports Editorial Board Member, provide suggested/excluded reviewer details if requested, and verify the final portal spelling of author and affiliation data. No personal declaration or reviewer identity was inferred or fabricated.
