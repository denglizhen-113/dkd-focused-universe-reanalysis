# v1.2.1 — M20.1 network-verification and provenance release

Date: 21 August 2026

This release does not change the M20 statistical rules, numerical results or operational pathway calls. It strengthens literature positioning, source-lineage auditing, disclosure and submission administration.

## Added

- Zhang et al. (2025), which used GSE30528 and GSE104948 for discovery and GSE96804 for external validation, and Jiao et al. (2022), which reported glomerular complement/C1q/C3 findings using GSE30528, GSE104948, GSE96804 and GSE99339.
- Supplementary Table S28 with exact GSE99339-to-GSE104948 mappings for H7 donors DN901, DN910, DN914, DN916, DN932, DN941 and DN947.
- GEO-declared series-level lineage warnings linking older-CDF analyses in GSE47183 and GSE32591 with GSE104948, without inventing unverified one-to-one sample mappings.
- Explicit AI accountability language: AI did not autonomously determine the retained inferential conclusions.
- Cover-letter placeholders for suggested reviewers, excluded referees and any prior discussion with a *Scientific Reports* Editorial Board Member.

## Verified

- Current and archived search evidence consistently resolve to 21 August 2026 in Asia/Shanghai; the GEO UTC timestamp converts to that local date.
- The complete PRISMA 2020 checklist retains 42 item/sub-item rows.
- The M20 numerical regression test passes unchanged.
- All seven primary result files are byte-identical when rerun from the extracted code archive under Python 3.12.12, NumPy 2.4.6, pandas 3.0.3 and SciPy 1.18.0.
- Figures 1–4 remain one-page vector PDFs and passed visual review.

## Author action required

`submission_ready_scientific_reports_m20/cover_letter_AUTHOR_COMPLETION_REQUIRED.pdf` must not be uploaded unchanged. The author must replace the three bracketed declarations, confirm the local secondary-analysis ethics determination and rename the file to `cover_letter.pdf`.
