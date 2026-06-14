---
gsd_state_version: 1.0
milestone: v3.0
milestone_name: Strain Model Redraw
status: Awaiting next milestone
stopped_at: Completed 09-02-PLAN.md
last_updated: "2026-06-14T21:07:27.007Z"
last_activity: 2026-06-14 — Milestone v3.0 completed and archived
progress:
  total_phases: 5
  completed_phases: 4
  total_plans: 14
  completed_plans: 13
  percent: 80
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-06-12 — v3.0 started)

**Core value:** A reported minor strain must be backed by orthogonal de novo/BLAST evidence — never on top of a failed major/dominant, never when de novo refutes it as a cross-mapping artefact — while genuine low-abundance co-infections are preserved.
**Current focus:** Phase 09 — compatibility-filename-migration-regression-suite

## Current Position

Phase: Milestone v3.0 complete
Plan: —
Status: Awaiting next milestone
Last activity: 2026-06-14 — Milestone v3.0 completed and archived

## Performance Metrics

**Velocity:**

- Total plans completed: 11 (v1.0) + 2 (v2.0 Phase 5, partial)
- Average duration: -
- Total execution time: 0 hours (v3.0)

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01 | 3 | - | - |
| 02 | 3 | - | - |
| 06 | 3 | - | - |
| 07 | 2 | - | - |

**Recent Trend:**

- Last 5 plans: -
- Trend: -

*Updated after each plan completion*
| Phase 01 P01 | 12min | 2 tasks | 5 files |
| Phase 01 P02 | 18min | 3 tasks | 9 files |
| Phase 01 P03 | 30min | 4 tasks | 8 files |
| Phase 02 P01 | 6min | 2 tasks | 3 files |
| Phase 02 P03 | 75min | 2 tasks | 5 files |
| Phase 03 P01 | 9min | 1 tasks | 2 files |
| Phase 03 P02 | 11min | 2 tasks | 2 files |
| Phase 03 P03 | 4min | 3 tasks | 8 files |
| Phase 04 P01 | 4min | 3 tasks | 7 files |
| Phase 04 P02 | 12m | 3 tasks | 12 files |
| Phase 04 P03 | 5min | 1 tasks | 1 files |
| Phase 05 P01 | 8min | 2 tasks | 8 files |
| Phase 06 P01 | 30min | 3 tasks | 5 files |
| Phase 06 P02 | 10min | 2 tasks | 3 files |
| Phase 06 P03 | 35min | 3 tasks | 3 files |
| Phase 07 P01 | 12m | 3 tasks | 3 files |
| Phase 07 P02 | ~4 min | 3 tasks | 4 files |
| Phase 08 P01 | 22min | 3 tasks | 5 files |
| Phase 08 P02 | 35min | 3 tasks | 4 files |
| Phase 09 P01 | 3min | 3 tasks | 4 files |
| Phase 09 P02 | ~7min | 3 tasks | 4 files |
| Phase 09 P03 | 25min | 2 tasks | 1 files |
| Phase 09 P04 | 40min | 3 tasks | 2 files |

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- [v3.0 Roadmap]: Phases continue at 6 (v1.0 Phases 1-4 shipped, v2.0 Phase 5 paused/preserved). v3.0 = Phases 6-9.
- [v3.0 Roadmap]: Dependency order enforced — classification (Phase 8) depends on BOTH neutral candidates (Phase 6) and assembly support (Phase 7) existing first.
- [v3.0 Roadmap]: COMPAT-02 (filename slot migration) and `bin/summarize.R` parsing kept lockstepped in a single phase (Phase 9, last/highest-risk).
- [Phase 3]: Match de novo↔minor at genotype level (parameterized), not subtype; asymmetric refute rule. (carries into ASUP-02 design)
- [Phase 3]: `--denovo_confirm_minor` defaults ON; OFF reproduces legacy behaviour. (carries into COMPAT-01)
- [Phase 01]: Minor-gate decision moved into summarize_mapping_to_all_references.R as minor_call/gate_flag columns (R-emits-decision, Nextflow-routes). (REFSEL-03 replaces this with per-candidate confirmation_status)
- [Phase 02]: Both `*.blastparse.csv` and `*_blast_out.csv` are wired; the latter carries the k-mer-cov floor (needed for ASUP-01).
- [Phase ?]: HCVGLUE is now a hermetic single-image per-BAM task; docker.sock mount removed (T-05-DS). (v2.0, paused)
- [Phase ?]: [Phase 06]: confirmation_status vocabulary = pass / below_threshold (per-candidate reads>minRead && cov>minCov); Phase 8 redefines against dominance scoring.
- [Phase ?]: [Phase 06]: Selection is mechanical read-recruitment ranking (top ref per distinct subtype, top-N by reads); is_valid_minor validity filtering deleted and moved to Phase 8 (D-05).
- [Phase ?]: [Phase 06]: Long-format *.candidates.csv emitted alongside reconstructed legacy wide CSV + _major.fa/_minor.fa shim (D-06); all 10 legacy columns always present, NA-filled.
- [Phase 06]: 06-02: Pinned PARSEFIRSTMAPPING legacy csv/major/minor globs to *.parsefirstmapping.csv so the new *.candidates.csv routes only to the candidates emit (T-06-05 glob-collision mitigation).
- [Phase 06]: 06-02: GATE-03/04 minor_mapping.size()==0 assertions kept unchanged — their fixtures are single-subtype, so neutral ranking yields one candidate regardless of D-05 validity-filter removal.
- [Phase 06]: 06-03: Collapsed MAJOR/MINOR_MAPPING aliases into ONE per-candidate TARGETED_MAPPING fan-out (splitCsv ALL rows); .major./.minor. filename slot derived from meta.candidate_rank in config so summarize.R parses unchanged (T-1).
- [Phase 06]: 06-03: Per-candidate routing on R-emitted confirmation_status=='pass' (no NA.toInteger); per-rank FASTA from legacy emits keeps meta.reference byte-identical; remainder:true joins so single-candidate samples survive.
- [Phase 06]: 06-03: D-03 behavior change confirmed (2nd candidate 3i_JX227955 vs old 4k_EU392173); golden reproduction deferred to Phase 9/COMPAT-01, NOT a Phase-6 gate.
- [Phase ?]: [Phase 08]: 08-01: Pitfall 1 resolved by reconciling nextflow.config to validated denovo floors 1000/2.0 (not re-validating 500/10.0); the old 10.0 k-mer floor would refute ERR1810453's genuine partial 2b.
- [Phase ?]: [Phase 08]: 08-01: Dominance score = wr*log10(reads) + we*breadth + we*cv_evenness + wk*log10(1+min(kmer,cap)); breadth+evenness share score_weight_evenness (default 3.0 > reads 1.0) so spiky high-read 4g loses to a genuine even minor (SCORE-02).
- [Phase ?]: [Phase 08]: 08-01: classify_roles.R is a pure sourced helper (score_candidates + classify_roles + verbatim-recovered is_valid_minor); summarize.R wiring deferred to Plan 02.
- [Phase ?]: [Phase 08]: 08-02: classify_roles() floor wired to min_targeted_read/min_targeted_cov (args[9]/[10] = configured minRead/minCov); fixed b23e021's undefined minRead/minCov call (Rule 1 bug).
- [Phase ?]: [Phase 08]: 08-02: legacy apply_denovo_layer/minor_denovo_status/coinfection_flag consumption retired (D-15); review_flag rewired onto per-sample role rollup + overall_sample_call; denovo_layer.R still sourced + unit-tested directly.
- [Phase ?]: [Phase 08]: 08-02: additive Major_role_*/Minor_role_* + overall_sample_call wide columns; legacy Major_*/Minor_* aliasing deferred to Phase 9 (COMPAT-03); enriched long candidates.csv surfaces every candidate incl. background (CLASS-03).
- [Phase ?]: [Phase 09]: 09-01: Migrated .major./.minor. -> uniform .cand{rank}. slot in all 6 TARGETED_MAPPING ext.prefix closures + N-candidate FASTA write loop in summarize_mapping_to_all_references.R (COMPAT-02 producer side).
- [Phase ?]: [Phase 09]: 09-01: PARSEFIRSTMAPPING major_mapping/minor_mapping collapsed to one N-FASTA candidate_fasta emit (*_cand*.fa); blastparse fasta emits collapsed too (no consumers); workflow fan-out rewire + N>2 guard lift deferred to Plan 02 (same-wave lockstep).
- [Phase ?]: [Phase 09]: 09-02: Workflow fan-out rewired to single PARSEFIRSTMAPPING.out.candidate_fasta N-FASTA emit; FASTA picked by rank-indexed _cand{rank}. basename match (replaces rank==1?major:minor two-slot pick); n_candidates>2 guard lifted (D-01).
- [Phase ?]: [Phase 09]: 09-02: All four anti-regression invariants preserved (remainder:true, rank-as-String, confirmation_status pass and fasta-not-null filter, id==sample assert); both module snapshots tool-regenerated with cand-slot, diff confined to slot rename, all non-filename md5s byte-identical (Pitfall 4).
- [Phase ?]: [Phase 09]: 09-03: summarize.R recovers candidate rank by joining cleaned candidate_ref to candidate_rank (candidates CSV) across all 3 stats loops + cv_by_ref + consensus-distance block; replaces filename position-3 first_major_minor parse (COMPAT-02 consumer side, D-02). Major_reference/Minor_reference cleaned via candidate_ref (no slot suffix) keep L366/coverage join keys byte-identical; flag-OFF golden baseline verified.
- [Phase ?]: test_compat.R drives the REAL summarize.R via system2 for COMPAT-01/02/03 and sources classify_roles.R for COMPAT-04 role_reason assertions
- [Phase ?]: [Rule 1] empty positional summarize.R args [4-8] passed as shQuote("") so system2 does not drop them and shift minRead/minCov off args[9]/[10]

### Pending Todos

- Refactor HCVGLUE module for per-sample parallel Docker execution — `.planning/todos/pending/2026-06-06-refactor-hcvglue-parallel-docker.md` (v2.0 Phase 5, paused)

### Blockers/Concerns

- [v3.0 Phase 8]: Combined dominance score weighting (breadth-evenness vs read count vs k-mer cov) is the highest-design-risk sub-task — must catch the benchmarked false-4g minor (refute to background) while preserving the genuine ERR1810447 2b co-infection. Calibrate against the handoff evidence table.
- [v3.0 Phase 9]: Filename-slot migration `.major.`/`.minor.` → `.cand1.`/`.cand2.` must land lockstep with `bin/summarize.R` parsing — a mismatch silently drops rows. Highest-risk integration; last phase.
- [v3.0 Phase 6]: `MAJOR_MAPPING`/`MINOR_MAPPING` aliases in `workflows/hcvtyper.nf` are consumed downstream by HCVGLUE feed wiring — confirm Phase 5 (paused) wiring is not broken by the per-candidate generalization; coordinate on resume.
- [v2.0 Phase 5]: PAUSED at Plan 03 Task 2 — blocked by `run-glue.sh` container defects (see `.continue-here.md`).

### Quick Tasks Completed

| # | Description | Date | Commit | Directory |
|---|-------------|------|--------|-----------|
| 260608-gef | Remove TANOTI process entirely from the pipeline | 2026-06-08 | 517b9eb | [260608-gef-remove-tanoti-process-entirely-from-the-](./quick/260608-gef-remove-tanoti-process-entirely-from-the-/) |
| 260609-dj8 | Add coinfection_flag QC warning column to Summary.csv | 2026-06-09 | 5feb8f2 | [260609-dj8-add-a-qc-warning-flag-to-summary-csv-for](./quick/260609-dj8-add-a-qc-warning-flag-to-summary-csv-for/) |
| 260609-odh | de novo BLAST subtype vs first-pass mapping subtype comparison + review_flag column | 2026-06-09 | 82df402 | [260609-odh-i-want-to-extract-the-subtype-from-the-n](./quick/260609-odh-i-want-to-extract-the-subtype-from-the-n/) |
| 260610-8t0 | Rewrite review_flag to human-readable sentences with co-infection vs single-infection distinction | 2026-06-10 | ff12009 | [260610-8t0-rewrite-review-flag-in-bin-summarize-r-t](./quick/260610-8t0-rewrite-review-flag-in-bin-summarize-r-t/) |

## Deferred Items

Items acknowledged and deferred at v3.0 milestone close on 2026-06-14:

| Category | Item | Status |
|----------|------|--------|
| verification_gap | Phase 08: 08-VERIFICATION.md [human_needed] — Docker runtime test pending (all 13 static checks pass; tidyverse not installed in host) | accepted — Docker test needed |
| verification_gap | Phase 09: 09-VERIFICATION.md [gaps_found] — both blockers fixed post-verification in commits 8321c18 / ee7f839; document not re-run | accepted — close in next session |
| uat_gap | Phase 09: 09-UAT.md [passed] — 0 pending scenarios (audit false-positive; status is already passed) | accepted |
| quick_task | 260608-gef-remove-tanoti-process-entirely (status: unknown in index — completed in commit 517b9eb) | accepted — already done |
| quick_task | 260609-odh-extract-subtype-from-nextflow (status: unknown in index — completed in commit 82df402) | accepted — already done |
| todo | refactor-hcvglue-parallel-docker — v2.0 Phase 5 paused work | accepted — carry to v4.0 |
| todo | remove-tanoti-mapper-completely — completed as quick task 260608-gef | accepted — already done |

Known deferred items at v3.0 close: 7 (see above)

Items acknowledged and carried forward from previous milestone close:

| Category | Item | Status | Deferred At |
|----------|------|--------|-------------|
| Feature | Change 3 — de novo-informed reference selection (REFSEL-01) | SHIPPED in v3.0 (Strain Model Redraw) | 2026-06-12 |

Items acknowledged and deferred at v1.0 milestone close on 2026-06-08:

| Category | Item | Status |
|----------|------|--------|
| todo | refactor-hcvglue-parallel-docker | in progress (v2.0 Phase 5, paused) |
| todo | remove-tanoti-mapper-completely | done (quick 260608-gef) |
| tech-debt | Phase 04 has no machine VERIFICATION.md — verified via 04-UAT + live R suite (see v1.0-MILESTONE-AUDIT.md) | accepted |
| tech-debt | Nyquist VALIDATION partial for phases 01/02/03 (nyquist_compliant: false) | accepted — run /gsd-validate-phase 1-3 later |
| tech-debt | Latent Phase-1 no_mapping FASTA-write crash (T-04-03) — documented, not fixed (D-03) | accepted — separate milestone |
| known-limitation | Extreme-ratio IVT minor (~300bp/~1x) is refuted (D-09) | accepted |

## Session Continuity

Last session: 2026-06-14T15:10:14.537Z
Stopped at: Completed 09-02-PLAN.md
Resume file: None

## Operator Next Steps

- Start the next milestone with /gsd-new-milestone
