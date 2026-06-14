# HCVTyper — De novo-Informed Strain Selection

## What This Is

HCVTyper is a Nextflow DSL2 pipeline that genotypes Hepatitis C virus (HCV) from NGS data and detects major/minor strain co-infections. As of **v3.0 (Strain Model Redraw)**, the pipeline replaces the overloaded major/minor model with three clearly separated axes: dominance-neutral candidate selection during the run, independent per-genotype assembly support (de novo/BLAST), and summary-step strain-role classification (dominant / co-infection / background) driven by a breadth-evenness-weighted dominance score. Index-hopping/contamination false co-infections are now caught by the "guilty until corroborated" classification rule, while genuine low-abundance co-infections are preserved.

For: the HCVTyper maintainers/analysts at FHI and downstream users running clinical/research HCV genotyping. Tied to manuscript revision Access Microbiology ACMI-D-26-00042.

## Core Value

A reported minor strain (co-infection) must be backed by orthogonal de novo/BLAST evidence — never reported on top of a failed major, and never reported when de novo refutes it as a cross-mapping artefact — while genuine low-abundance co-infections are preserved.

## Paused: v2.0 HCVGLUE Refactor

Phase 05 is paused mid-execution (Plan 03 Task 2, blocked by `run-glue.sh` container defects — see `.continue-here.md`). GLUE-02..05 remain Active. Resume via `/gsd-resume-work` (Phase 5 directory + `.continue-here.md` intact).

## Requirements

### Validated

<!-- Inferred from existing codebase (brownfield) -->

- ✓ QC/trim (fastp|cutadapt) + Kraken2 HCV read selection — existing
- ✓ First-pass mapping (bowtie2|tanoti) of HCV reads vs full reference set — existing
- ✓ Major/minor strain selection from read-count + coverage thresholds (`parsefirstmapping`) — existing
- ✓ De novo assembly (SPAdes) + BLAST vs reference set, parsed to per-sample major/minor genotype (`blastparse` / `blast_parse.R`) — existing, **currently output-only**
- ✓ Targeted re-mapping + iVar consensus + HCV-GLUE genotyping/RAVs — existing
- ✓ Cross-sample contamination detection (all-vs-all BLAST) — existing
- ✓ Result aggregation (SUMMARIZE) + MultiQC report — existing
- ✓ 1a/1b co-infection exception and 2k/1b recombinant suppression in selection logic — existing

- ✓ **Change 1 — Major-gate** (GATE-01..06) — v1.0. A minor is only evaluated when the major passes `minRead`+`minCov`; failed-major reports first-mapping stats + a `gate_flag` reason, no minor/major call. Fixes ERR1810469.
- ✓ **Change 2 — De novo confirmation** (CONF-01..07) — v1.0. Confirm/refute/fall-back of the candidate minor against calibrated substantial-contig thresholds, genotype-level match, asymmetric refute rule.
- ✓ **Configurable + non-breaking** (CONF-07) — v1.0. `--denovo_confirm_minor` default ON; flag-OFF reproduces legacy output (golden-baseline guarded).
- ✓ **Reporting status field** (REPORT-01) — v1.0. `minor_denovo_status` = confirmed_by_denovo / unconfirmed / refuted / not_evaluated in both output CSVs.
- ✓ **Preserved exceptions** (CONF-06) — v1.0. 1a/1b co-infection + 2k1b recombinant handling intact after both changes.
- ✓ **Tests / regression guard** (TEST-01, TEST-02) — v1.0. `bin/tests/run_all.sh` covers major-gate, all three branches, genotype-match, flag-OFF reproduction, and non-suppression of genuine co-infections; wired into CI.

- ✓ **Per-genotype assembly support** (ASUP-01, ASUP-02) — v3.0 Phase 7. De novo/BLAST evidence is reframed as a dominance-neutral per-subtype roll-up (single best contig by `sc_length` → best contig length, BLAST identity, BLAST alignment length, k-mer coverage) emitted as `*.assembly_support.csv` from `blast_parse.R`, then left-joined to the neutral candidates at genotype level (parameterized `denovo_match_level`, default genotype) via the sourceable `assembly_support_join.R`. Candidates are the LEFT side (no row loss); unmatched candidates resolve to explicit `assembly_support = "none"` with NA metrics.
- ✓ **Neutral candidate selection** (REFSEL-01, REFSEL-02, REFSEL-03) — v3.0 Phase 6. `is_valid_minor()` validity filtering deleted; `params.n_candidates` (default 2) drives uniform read-recruitment ranking; single `TARGETED_MAPPING` fan-out via `splitCsv.flatMap`; `confirmation_status` per candidate replaces `gate_flag`/`minor_call`.
- ✓ **Dominance score + strain-role classification** (SCORE-01, SCORE-02, CLASS-01, CLASS-02, CLASS-03, CLASS-04) — v3.0 Phase 8. Breadth-evenness-weighted dominance score (`score_weight_evenness=3.0 > score_weight_reads=1.0`); guilty-until-corroborated classification (dominant/co-infection/background) with `overall_sample_call`; false 4g refuted to background; genuine 2b co-infection preserved; `classify_roles.R` pure sourced helper in SUMMARIZE.
- ✓ **Compatibility, filename migration + regression suite** (COMPAT-01, COMPAT-02, COMPAT-03, COMPAT-04, TEST-01) — v3.0 Phase 9. `.major.`/`.minor.` → `.cand{rank}.` lockstep migration; legacy `Major_*`/`Minor_*` aliases retained one release; `test_compat.R` golden reproduction + exception coverage; wired into CI.

### Active

<!-- v2.0 scope: HCVGLUE Refactor — PAUSED -->
- [ ] **GLUE-02**: Concurrent HCVGLUE tasks do not collide — each task is isolated and parallel-safe
- [ ] **GLUE-03**: Output format preserved — `*.json` + `*.html` per sample, compatible with existing GLUE_PARSER
- [ ] **GLUE-04**: Process works under Docker, Singularity/Apptainer, and Podman profiles
- [ ] **GLUE-05**: `bin/run_hcvglue.sh` simplified or replaced; no apt-get inside containers, no deprecated `--link` flag

<!-- Post-v3.0 candidates for next milestone -->
- [ ] **COMPAT-03-DROP**: Drop the legacy `Major_*`/`Minor_*` column aliases after the one-release deprecation window
- [ ] **NYQUIST-01-03**: Close Nyquist validation for phases 01–03 (nyquist_compliant: false in VALIDATION.md)
- [ ] **CRASH-01**: Fix the latent Phase-1 `no_mapping` FASTA-write crash (T-04-03, documented)

### Out of Scope
- **Full manuscript re-benchmarking as part of DoD** — re-running sim + real datasets and diffing against `combined_analysis.tsv` is a separate manual analyst step. This project ships code + tests; the analyst validates afterward.
- **Deciding co-infection vs contamination biology** — when a minor is confirmed genuine, whether it represents co-infection or contamination is a biology question, not a pipeline one.
- **Changing `minRead`/`minCov` thresholds as the fix** — benchmarking proved thresholds cannot separate artefact minors from genuine ones; the discriminator must be orthogonal (de novo), not a stricter cutoff.

## Context

- **Brownfield.** Codebase already maps cleanly (`.planning/codebase/`). Core selection logic lives in `modules/local/parsefirstmapping/` (with R helpers `summarize_mapping_to_all_references.R`); de novo evidence is parsed in `modules/local/blastparse/` (`bin/blast_parse.R`) producing per-sample `major_ref`/`minor_ref` genotype calls. Final reporting is `modules/local/summarize/` (`bin/summarize.R`).
- **Architecture note / correction:** the handoff guessed the parse scripts were Python; they are **R** (`blast_parse.R`, `summarize_mapping_to_all_references.R`, `summarize.R`). Verify exact selection/reporting touchpoints during planning.
- **The enabling fact:** `blastparse` already reduces de novo contigs to a major/minor genotype call per sample. Change 2 is largely about *consuming* existing output at decision time, not new computation.
- **Benchmark evidence (why this matters):**
  - Failure A: sim 1a single-infection → first-pass flagged a false 4g minor (53,279 reads, 67.7% breadth, 1,212× depth) clearing all thresholds; de novo produced only 1a contigs, no 4g → refute.
  - Failure B: ERR1810469 (3a) → major 3a failed coverage (8.8% breadth, 3.5× depth) yet a 1a minor was still called → major-gate fixes it.
  - Corollary: most flagged minors in single-genotype samples are REAL (de novo recovered full second-genotype genomes, e.g. ERR1810447 full 9,207 bp 2b). Goal is confirm genuine + refute artefacts, NOT suppress all minors.
- **Datasets (read-only mounts, for the analyst's later re-benchmark):** sim `/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2026/HCVTyper_sim_data/`; real `/mnt/N/.../HCVTyper_SRA_data/`; truth+results `/home/jon.brate/HCV_nextflow_pipeline/data/thomson2016/combined_analysis/combined_analysis.tsv`.

## Constraints

- **Tech stack**: Nextflow DSL2 + nf-core conventions; HCV-specific logic in local modules with R/Python helper scripts in `bin/`. New logic should follow existing module/script patterns.
- **Compatibility**: Non-breaking — existing behaviour must be reproducible via the `--denovo_confirm_minor` flag (default ON). Existing 1a/1b and 2k/1b exceptions must be preserved.
- **Discriminator design**: Minor confirmation must match at **genotype** level, not subtype (subtype matching is too strict — cf. 3b→3a). Match level should be a parameter.
- **"Substantial contig" calibration**: A minimum contig length (≈≥1,000 bp or ≥X% of genome) plus possibly minimum BLAST identity/length and a k-mer-coverage floor, tuned to confirm partial genuine contigs (e.g. ERR1810453's 2,949 bp 2b) while ignoring ~300 bp / ~1× noise contigs.
- **Timeline**: Manuscript resubmission deadline 2026-08-11 — implementing + releasing before then lets future-work claims be upgraded to features, but the manuscript stands as-is if not.

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Scope = Changes 1 + 2 only; defer Change 3 | Phased by risk; 1+2 are lower-risk and fix the two benchmarked failures. Change 3 needs its own validation. | ✓ Good — shipped v1.0; Change 3 carried to v2 |
| Refuted minor → report as single infection + QC flag | Don't report an unconfirmed co-infection; surface the candidate in QC status fields for analyst visibility. | ✓ Good — downgrade-only mutation; Minor_* columns preserved |
| De novo-confirmation flag defaults ON | Make the improved behaviour standard; old behaviour reproducible by toggling off. | ✓ Good — flag-OFF golden-baseline reproduction guarded in CI |
| Re-benchmarking NOT in DoD | Ship code + tests; full dataset re-benchmark is a separate manual analyst step. | ✓ Good — shipped with R regression suite; analyst re-benchmark pending |
| Match de novo↔minor at genotype level (parameterized) | Subtype matching too strict; genotype level handles divergent/novel subtypes (3b→3a). | ✓ Good — `denovo_match_level` param, default genotype |
| Asymmetric refute rule | Absence of a minor contig only refutes when de novo otherwise assembled a substantial major contig; protects genuine low-yield co-infections. | ✓ Good — verified across evidence-table fixtures |
| Remove the `strategy='denovo'` branch + `strategy` param entirely (Phase 1) | Dead/non-functional (`minDenovoLength` undeclared), never default, in no config profile, undocumented; it's a crude Change-3 prototype that bypasses all selection safeguards — the #1 pitfall flagged in research. | ✓ Good — removed at all 4 sites; security re-verified |
| Verify Phase 4 via UAT + live R suite (no machine VERIFICATION.md) | gsd-verifier never ran on Phase 4; UAT 5/5 + a live `run_all.sh` pass is stronger evidence than a re-derivation. | ⚠️ Revisit — generate the verifier artifact retroactively in v2 if needed |
| Delete `is_valid_minor()` from selection; port to `classify_roles.R` (v3.0 D-05) | Validity rules (different-genotype / 1a-1b / 2k1b) belong at classification time, not selection time. | ✓ Good — classification is now the single place for HCV exception logic |
| Reconcile denovo floors to validated 1000/2.0/90 (v3.0 Pitfall 1) | Old 10.0 k-mer floor would refute ERR1810453's genuine partial 2b (k-mer ~5); 2.0 passes it correctly. | ✓ Good — test_classify_roles.R Test 2b locks this against regression |
| Breadth-evenness weighting 3× reads in dominance score (v3.0 SCORE-02) | Read counts alone cannot separate index-hopping contamination from genuine co-infection; breadth uniformity is the discriminator. | ✓ Good — 4g artefact loses to genuine even minor even at 12× fewer reads |
| Lockstep `.major.`/`.minor.` → `.cand{rank}.` migration in single phase (v3.0 COMPAT-02) | Producer (6 ext.prefix closures + FASTA loop) and consumer (`summarize.R` join) must change atomically; a split risks silent empty rows. | ✓ Good — one-phase lockstep confirmed by integration checker (0 broken flows) |
| Legacy `Major_*`/`Minor_*` aliases for one release (v3.0 COMPAT-03) | Downstream consumers read old column names; breaking change without a deprecation window is unsafe. | ✓ Good — drop in next milestone after deprecation window |

## Current State

**Shipped: v3.0 — Strain Model Redraw (2026-06-14).** 4 phases (6-9), 11 plans, 45 commits over 3 days. The pipeline now classifies each candidate into dominant / co-infection / background at the summary step, driven by a breadth-evenness-weighted dominance score. The false 4g artefact is refuted to background; genuine 2b co-infections are preserved. Filename slots migrated `.major.`/`.minor.` → `.cand{rank}.`; legacy `Major_*`/`Minor_*` columns aliased for one release; regression suite covers all four COMPAT requirements. 16/16 requirements satisfied.

**Tech debt from v3.0 close:** Phase 8 Docker runtime test pending (run `bash bin/tests/run_all.sh` in the pipeline container); nf-test workflow snapshot needs `--update-snapshot`; `modules_hcv.config` L123-129 has a stale comment.

**Shipped: v1.0 — De novo Minor Confirmation (2026-06-08).** 4 phases, 12 plans. Changes 1 + 2 are live and default-ON.

- **What works now:** de novo/BLAST evidence gates minor-strain reporting end-to-end; major-gate blocks minors behind a failed major; `minor_denovo_status` surfaces the basis for each call; `--denovo_confirm_minor false` reproduces legacy output exactly.
- **Safety net:** `bin/tests/run_all.sh` (4 R test files) covers the major-gate, confirm/refute/fall-back branches, flag-OFF golden reproduction, and genuine-co-infection non-suppression — wired into CI as the `r-regression` job. All 4 phases security-verified (`*-SECURITY.md`, `threats_open: 0`).
- **Accepted tech debt (see `.planning/v1.0-MILESTONE-AUDIT.md`):** Phase 4 lacks a machine VERIFICATION.md (verified via UAT + live suite); Nyquist validation partial for phases 01–03; one documented latent Phase-1 `no_mapping` FASTA-write crash; one known extreme-ratio IVT refute limitation (D-09).
- **Manuscript:** ACMI-D-26-00042; resubmission deadline 2026-08-11. v1.0 lets future-work claims be upgraded to shipped features.

## Post-v2.0 Candidates

- ▶ Change 3 — de novo-informed reference selection (REFSEL-01) — **promoted to active milestone v3.0 (Strain Model Redraw), 2026-06-12.**
- Close Nyquist validation (phases 01–03) + fix the latent Phase-1 `no_mapping` crash.
- Finish v2.0 HCVGLUE refactor (Phase 5, paused) — see `.continue-here.md`.

## Evolution

This document evolves at phase transitions and milestone boundaries.

**After each phase transition** (via `/gsd-transition`):
1. Requirements invalidated? → Move to Out of Scope with reason
2. Requirements validated? → Move to Validated with phase reference
3. New requirements emerged? → Add to Active
4. Decisions to log? → Add to Key Decisions
5. "What This Is" still accurate? → Update if drifted

**After each milestone** (via `/gsd-complete-milestone`):
1. Full review of all sections
2. Core Value check — still the right priority?
3. Audit Out of Scope — reasons still valid?
4. Update Context with current state

---
*Last updated: 2026-06-14 after v3.0 milestone — all REFSEL/ASUP/CLASS/SCORE/COMPAT/TEST requirements moved Active → Validated. v2.0 HCVGLUE Refactor remains paused at Phase 5.*
