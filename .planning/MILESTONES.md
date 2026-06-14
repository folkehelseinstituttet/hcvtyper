# Milestones

## v3.0 Strain Model Redraw (Shipped: 2026-06-14)

**Phases completed:** 5 phases, 14 plans, 24 tasks

**Key accomplishments:**

- HCVGLUE refactored to a per-BAM hermetic single-container task (custom MySQL 5.7 + GLUE image, runtime-staged params.hcvglue_db) with the host docker.sock mount and bin/run_hcvglue.sh both eliminated.
- workflows/hcvtyper.nf rewired so each BAM feeds one isolated HCVGLUE task — the SQL dump is staged once as a broadcast value channel and combined per-emission, the bulk-collect feed is dropped (GLUE-02/GLUE-03), and versions use .first() to keep the snapshot deterministic.
- 1. [Rule 2 - Missing Critical] Guarded empty/X1-less depth frame (T-06-01 DoS mitigation)
- 1. [Rule 3 - Blocking] Verify-block disk guard tripped at 94-95% (host root fs genuinely near-full)
- 1. [Rule 1 - Bug] `remainder:true` on the optional legacy-emit joins to avoid dropping single-candidate samples
- Task 1 — §4b roll-up in `bin/blast_parse.R` (commit 88e0884)
- Task 1 — sourceable join helper `bin/assembly_support_join.R` (commit 196dfdf)
- Pure, unit-testable bin/classify_roles.R encoding the breadth-evenness-weighted dominance score + D-01..D-14 strain-role logic + verbatim-recovered is_valid_minor() exceptions, with the shipped denovo floor reconciled to the validated 1000/2.0/90 and the false-4g→background / genuine-2b→co-infection evidence-table cases asserted against the real helper.
- summarize.R now runs the Phase-8 N-candidate role classifier over the Phase-7-joined candidate frame — emitting per-candidate role/dominance_score/role_reason + one overall_sample_call per sample, retiring the legacy apply_denovo_layer/minor_denovo_status/coinfection_flag path (D-15), rewiring review_flag onto roles, writing the enriched long candidates.csv, and staging classify_roles.R into the SUMMARIZE module.

---

## v1.0 De novo Minor Confirmation (Shipped: 2026-06-08)

**Phases completed:** 4 phases, 12 plans, 19 tasks
**Git range:** `refactor(01-01)` → `feat(04-03)` (2026-06-06 → 2026-06-07)
**Audit:** tech_debt — 20/20 requirements satisfied, no functional blockers (see milestones/v1.0-MILESTONE-AUDIT.md)

**Delivered:** De novo/BLAST evidence and a major-gate now drive minor-strain reporting end-to-end, behind `--denovo_confirm_minor` (default ON), fixing the two benchmarked failure modes while preserving genuine co-infections.

**Key accomplishments:**

- **Change 1 — Major-gate (GATE-01..06):** minor strains are only evaluated when the major passes `minRead`+`minCov`; the decision moved into `summarize_mapping_to_all_references.R` as `minor_call`/`gate_flag` CSV columns (R-emits, Nextflow-routes), killing the `NA.toInteger()` crash class. Fixes ERR1810469.
- **Selection-script hardening:** removed the dead `strategy='denovo'` branch + undeclared `params.minDenovoLength` (schema + all 3 config profiles); single-sourced the 2k1b-aware `genotype_from_subtype()` helper; fixed the always-truthy `length(minor_ref > 0)` line-144 bug.
- **De novo evidence plumbing (PLUMB-01..04):** BLASTPARSE per-contig CSVs wired into SUMMARIZE via `left_join` (NA-fill, no row loss), 5 typed `denovo_*` params declared everywhere; proven zero-behaviour-change via an additive-only column-subset diff.
- **Change 2 — De novo confirmation (CONF-01..07, REPORT-01):** pure `classify_minor_denovo()` core turns per-contig BLAST evidence into confirmed_by_denovo / refuted / unconfirmed via a calibrated substantial-contig test + asymmetric refute rule at genotype level; downgrade-only `apply_denovo_layer()` surfaces `minor_denovo_status` and never nulls Minor_* columns; 1a/1b + 2k1b exceptions preserved.
- **Regression guard (TEST-01, TEST-02):** `bin/tests/run_all.sh` (4 R test files) covers the major-gate, all three confirm/refute/fall-back branches, flag-OFF golden-baseline reproduction, and non-suppression of genuine co-infections — wired into CI as the `r-regression` job.
- **Security:** all 4 phases verified (`*-SECURITY.md`, `threats_open: 0`), including the CI surface (pinned image, no PR-input interpolation).

**Known deferred items at close: 6** (see STATE.md Deferred Items) — 2 v2 todos (Tanoti removal, HCVGLUE parallel Docker), Change 3 (REFSEL-01), plus accepted tech debt: Phase-4 machine VERIFICATION.md absent (verified via UAT + live suite), Nyquist partial for phases 01–03, latent Phase-1 no_mapping crash (T-04-03), extreme-ratio IVT refute limitation (D-09).

---
