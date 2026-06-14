---
phase: 08-dominance-scoring-strain-role-classification
plan: 01
subsystem: testing
tags: [R, tidyverse, dominance-score, strain-role, classify_roles, is_valid_minor, nextflow-config, nf-schema]

# Dependency graph
requires:
  - phase: 06-neutral-candidate-selection
    provides: long-format *.candidates.csv (candidate_subtype/genotype/reads/cov) — the classifier input
  - phase: 07-per-genotype-assembly-support
    provides: genotype-level assembly_support_* join columns — the corroboration evidence
provides:
  - Pure sourced helper bin/classify_roles.R hosting score_candidates() + classify_roles() + verbatim-recovered is_valid_minor()
  - Calibrated dominance score where breadth-evenness dominates raw read count (the benchmarked false 4g loses to a genuine even minor)
  - D-01..D-14 strain-role logic (gated dominant, asymmetric refute, HCV exceptions, 3-value sample call)
  - Reconciled denovo corroboration floor (1000/2.0/90) + new score_weight_* params in nextflow.config + nextflow_schema.json
  - Two Rscript unit-test files asserting the handoff §2 evidence-table outcomes against the real helper
affects: [Plan 08-02 summarize.R wiring, Phase 9 COMPAT/TEST regression suite]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Pure sourced R helper (functions only; no commandArgs/file I/O; defensive tidyverse guard)"
    - "Typed zero-row / NULL guard generalized to new score/role columns (T-08-01 / CLASS-03)"
    - "Calibration-as-test: evidence-table fixtures drive the score-weight defaults"

key-files:
  created:
    - bin/classify_roles.R
    - bin/tests/test_dominance_score.R
    - bin/tests/test_classify_roles.R
  modified:
    - nextflow.config
    - nextflow_schema.json

key-decisions:
  - "Pitfall 1 resolved by RECONCILING config to validated 1000/2.0 (Open Question #1 option a), not by re-validating 500/10.0"
  - "Score components: wr*log10(reads) + we*breadth_frac + we*cv_evenness + wk*log10(1+min(kmer,cap)); breadth and cv_evenness share score_weight_evenness (the breadth-evenness headline)"
  - "Default weights evenness=3.0 / reads=1.0 / kmercov=0.5 / evenness_k=1.0 — calibrated so 4g (low evenness) loses to a genuine even minor even with no de novo boost"
  - "is_valid_minor() reconstructed as a pure 4-arg function (cand/dom subtype+genotype); 2k1b demotion reported as recombinant_2k1b, generic same-genotype as same_genotype_as_dominant"
  - "role_reason vocabulary: dominant / corroborated / refuted_denovo / uncorroborated_kept / below_floor / same_genotype_as_dominant / recombinant_2k1b"

patterns-established:
  - "Score-then-classify pipeline: score_candidates() emits dominance_score, classify_roles() consumes it (defensive re-score if absent)"
  - "Per-sample grouping (group_split by sampleName) so dominant determination + asymmetric refute are sample-local"

requirements-completed: [SCORE-01, SCORE-02, CLASS-01, CLASS-02, CLASS-03, CLASS-04]

# Metrics
duration: 22min
completed: 2026-06-13
---

# Phase 8 Plan 01: Dominance Scoring + Strain-Role Classification Summary

**Pure, unit-testable bin/classify_roles.R encoding the breadth-evenness-weighted dominance score + D-01..D-14 strain-role logic + verbatim-recovered is_valid_minor() exceptions, with the shipped denovo floor reconciled to the validated 1000/2.0/90 and the false-4g→background / genuine-2b→co-infection evidence-table cases asserted against the real helper.**

## Performance

- **Duration:** ~22 min
- **Started:** 2026-06-13 (execution)
- **Completed:** 2026-06-13
- **Tasks:** 3
- **Files modified:** 5 (2 config, 3 new R)

## Accomplishments
- Reconciled the SHIPPED de novo corroboration floor to the 03-RESEARCH-validated 1000/2.0/90 (Pitfall 1), and declared the four new `score_weight_*` / `score_evenness_k` params in both `nextflow.config` and `nextflow_schema.json`.
- Authored the pure sourced helper `bin/classify_roles.R` — `score_candidates()` (SCORE-01/02), `classify_roles()` (D-01..D-14), and the verbatim-recovered `is_valid_minor()` (D-12) — following the `assembly_support_join.R` / `denovo_confirm.R` house pattern (functions only, defensive tidyverse guard, typed zero-row guard, match_level stopifnot).
- Calibrated the score weights so the benchmarked false 4g (53,279 reads, 67.7% breadth, low cv_evenness) scores BELOW a genuine even minor with far fewer reads, and so a missing k-mer cov never penalises a genuine low-yield dominant (D-05).
- Authored two Rscript unit-test files asserting the headline evidence-table outcomes against the REAL sourced helper: false-4g→background/refuted_denovo, full+partial 2b→co-infection/corroborated (incl. a lock that the old 10.0 floor would have refuted the partial 2b), sim1 1a:1b + sim2 2a:3a preserved, IVT→uncorroborated_kept, D-12 same_genotype + recombinant_2k1b demotions, no-gate→indeterminate, one-role-per-candidate, and zero-row/NULL typed-frame guards.
- `bash bin/tests/run_all.sh` is green (all pre-existing R tests plus the two new files).

## Task Commits

Each task was committed atomically:

1. **Task 1: Reconcile denovo_* floors + declare score_weight_* params** - `ccff05d` (chore)
2. **Task 2: Author bin/classify_roles.R** - `9b76de5` (feat)
3. **Task 3: Author the two unit-test files** - `55ff571` (test)

_Note: although the tasks carry `tdd="true"`, this plan delivers a standalone helper whose RED/GREEN distinction is the helper-vs-test split; the helper was committed (feat) before the tests (test). Calibration converged on the first attempt — no weight re-tuning was needed (the Task-1 defaults already make 4g lose to the even minor)._

## Files Created/Modified
- `bin/classify_roles.R` - Pure sourced helper: `score_candidates()`, `classify_roles()`, `is_valid_minor()`.
- `bin/tests/test_dominance_score.R` - SCORE-01/02 unit coverage incl. the 4g-loses-to-even-minor headline + D-05 bonus-only check.
- `bin/tests/test_classify_roles.R` - CLASS-01..04 + D-11/D-12 unit coverage anchored to the handoff §2 evidence table.
- `nextflow.config` - denovo_min_contig_length 500→1000, denovo_min_kmer_cov 10.0→2.0; declared score_weight_evenness/reads/kmercov + score_evenness_k.
- `nextflow_schema.json` - mirrored denovo defaults; added four typed `number` properties for the new params.

## Decisions Made
- **Pitfall 1 / Open Question #1:** chose option (a) — reconcile the config to the validated 1000/2.0 — over (b) re-validating 500/10.0. The k-mer floor of 10.0 is stricter than the validated 2.0 and would refute ERR1810453's genuine partial 2b (k-mer cov ~5); Test2b in `test_classify_roles.R` explicitly locks this.
- **Score form:** breadth fraction and the cv_evenness factor BOTH carry `score_weight_evenness` (the combined "breadth-evenness" headline of D-02/D-03); reads carry `score_weight_reads`; the capped log10(1+kmer) bonus carries `score_weight_kmercov`. Defaults 3.0 / 1.0 / 0.5 make evenness strictly dominant.
- **`score_candidates()` input flexibility:** prefers an explicit `cand_cov_breadth` column, falls back to `candidate_cov`, coerces 0–100 percent to a 0–1 fraction; accepts a precomputed `cv_evenness` factor (production, computed in the Plan-02 cov loop) or a raw `cv` transformed via `1/(1+k*CV)`. The helper does NOT read depth files (purity contract).
- **D-12 reason coding:** a 2k1b-involving demotion reports `recombinant_2k1b`; any other same-genotype demotion reports `same_genotype_as_dominant`. The exception may only demote a would-be co-infection to background, never promote.

## Deviations from Plan

None - plan executed exactly as written. The Task-1 score-weight defaults made the evidence-table assertions pass on the first calibration attempt, so no Task-3-driven weight re-tuning was required (the plan explicitly permitted adjusting the weights if needed).

## Issues Encountered
None. `Rscript -e 'source(... classify_roles.R)'` emits the tidyverse attach banner to stderr because the standalone Rscript has not pre-loaded tidyverse, so the defensive `if (!exists("group_by"))` guard fires — this is expected (the production consumer `summarize.R` already loads tidyverse) and the source exits 0.

## Known Stubs
None. The helper is complete and pure. The only deliberately deferred wiring is the `summarize.R` consumption (arg parse, cov-loop `cv_evenness` computation, `review_flag` rewire, enriched `*.candidates.csv` emit) and the `conf/modules_hcv.config` positional `ext.args` append — all explicitly Plan 08-02 per the plan scope, NOT stubs in this plan's artifacts.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- `bin/classify_roles.R` is ready to be sourced + wired into `bin/summarize.R` in Plan 08-02 (after `genotype_utils.R`, mirroring the `assembly_support_join.R` staging).
- Plan 08-02 must: append `score_weight_*` to the `SUMMARIZE` `ext.args` (END only — positional coupling), parse them at `args[12+]`, compute `cv_evenness` inside the cov loop (`cov$X3`, with the zero-mean guard), run score→classify over the Phase-7-joined candidate frame, fill the wide Major_/Minor_ slots + the new `overall_sample_call`, rewire `review_flag` onto roles (D-15), and emit the enriched long `*.candidates.csv` (D-16). The SUMMARIZE nf-test snapshot regen remains deferred (MEMORY: Phase 7 snapshot regen deferred; needs docker).

## Self-Check: PASSED
- All three created R files present on disk.
- SUMMARY.md present on disk.
- All three task commits (ccff05d, 9b76de5, 55ff571) present in git history.
