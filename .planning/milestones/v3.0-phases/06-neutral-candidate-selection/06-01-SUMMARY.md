---
phase: 06-neutral-candidate-selection
plan: 01
subsystem: selection
tags: [nextflow, nf-schema, R, tidyverse, idxstats, candidate-ranking, shim]

# Dependency graph
requires:
  - phase: 04-baseline
    provides: minor_call/gate_flag R-emits-decision selection script + R regression guard
  - phase: 03-denovo-confirm
    provides: genotype-level matching precedent (genotype_from_subtype, parameterized match level)
provides:
  - "params.n_candidates typed integer param (default 2, minimum 1) across schema + nextflow.config + modules_hcv.config"
  - "Neutral read-recruitment candidate ranking (top ref per distinct subtype, top-N subtypes by reads) with NO validity filtering"
  - "Long-format candidate table (<sample>.candidates.csv) with per-candidate confirmation_status — the Plan 02/03 downstream contract"
  - "Legacy 10-column wide CSV + _major.fa/_minor.fa reconstructed behind the D-06 compatibility shim"
  - "Empty-depth / no-mapping crash guard (T-06-01) so the no_mapping default is actually reachable"
affects: [phase-07-assembly-support, phase-08-classification, phase-09-compat-migration, blastparse-emit, hcvtyper-workflow-fanout]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Typed-param triple: declare new params in nextflow_schema.json + nextflow.config + conf/modules_hcv.config, mirroring the denovo_* precedent"
    - "Long-format candidate table emitted alongside a reconstructed legacy wide CSV (interface-first; downstream locks the long schema)"
    - "Distinct-subtype dedup via genotype_from_subtype() instead of is_valid_minor() validity rules"

key-files:
  created: []
  modified:
    - nextflow_schema.json
    - nextflow.config
    - conf/modules_hcv.config
    - bin/summarize_mapping_to_all_references.R
    - bin/tests/test_candidate_selection.R

key-decisions:
  - "confirmation_status vocabulary: pass / below_threshold, derived from the generalized reads>minRead && cov>minCov per-candidate comparison (Open Q2 resolution)"
  - "Long-format CSV emitted to a distinct *.candidates.csv glob so it never collides with the legacy *.parsefirstmapping.csv wide CSV"
  - "n_candidates coerced with as.integer; NA or <1 falls back to the safe default 2 (T-06-02 hardening)"

patterns-established:
  - "Selection is now mechanical read-recruitment ranking only; validity rules (different-genotype / 1a-1b / 2k1b) move to Phase 8 classification (D-05)"
  - "All 10 legacy columns are ALWAYS present (NA-filled) so bin/summarize.R can pull them unconditionally"

requirements-completed: [REFSEL-01, REFSEL-02]

# Metrics
duration: 30min
completed: 2026-06-12
---

# Phase 6 Plan 01: Neutral Candidate Selection Summary

**Rewrote the selection script to rank candidates neutrally by read recruitment (top ref per distinct subtype, top-N subtypes), emit a long-format `*.candidates.csv` with per-candidate `confirmation_status`, and reconstruct the full legacy 10-column CSV + `_major.fa`/`_minor.fa` behind the D-06 shim — with `is_valid_minor()` validity filtering deleted.**

## Performance

- **Duration:** ~30 min productive (across two sessions; a usage-limit pause separated Task 2 from Task 3)
- **Started:** 2026-06-12T13:35:01Z
- **Completed:** 2026-06-12T16:51:01Z (Task 3 commit)
- **Tasks:** 3
- **Files modified:** 5

## Accomplishments
- `params.n_candidates` declared as a typed integer (default 2, minimum 1) across the schema/config/modules_hcv triple, mirroring the `denovo_*` precedent (REFSEL-02).
- Selection script generalized to uniform read-recruitment ranking: total reads per subtype orders the subtypes, the top reference per distinct subtype is the per-rank pick, `head(n = n_candidates)` truncates. Distinct-subtype dedup (D-02) preserves the 1a/1b co-infection without an explicit allow-rule (REFSEL-01).
- `is_valid_minor()` and all validity filtering deleted from selection — the different-genotype / 1a-1b-allow / 2k1b-block rules now belong to Phase 8 classification (D-05).
- New long-format `*.candidates.csv` (one row per candidate; `sample, candidate_rank, candidate_ref, candidate_subtype, candidate_genotype, candidate_reads, candidate_cov, confirmation_status`) is the locked downstream contract for Plans 02/03.
- Legacy 10-column wide CSV (rank1→major_*, rank2→minor_*) plus `_major.fa`/`_minor.fa` reconstructed behind the shim, all columns always present with NA-fill on single-candidate samples (D-06).
- New auto-globbed R unit test `bin/tests/test_candidate_selection.R` encodes the ranking, dedup, no-filter, long-format, and shim contract; full `run_all.sh` suite is green (0 failures, no regressions).

## Task Commits

Each task was committed atomically:

1. **Task 1: Declare params.n_candidates across the typed-param triple** - `52fda8d` (feat)
2. **Task 2: R unit test for neutral ranking + shim reconstruction (RED)** - `ea71811` (test)
3. **Task 3: Rewrite selection script — neutral ranking + long-format + shim (GREEN)** - `43904de` (feat)

_Note: This plan is TDD — Task 2 was the RED test, Task 3 made it GREEN._

## Files Created/Modified
- `nextflow_schema.json` - Typed `n_candidates` property (integer, default 2, minimum 1) in the denovo_* group.
- `nextflow.config` - `n_candidates = 2` in the params block.
- `conf/modules_hcv.config` - `n_candidates = 2` in the params{} block.
- `bin/summarize_mapping_to_all_references.R` - PRIMARY REWRITE: neutral ranking, long-format emit, confirmation_status, legacy shim reconstruction, empty-depth guard; `is_valid_minor()` deleted.
- `bin/tests/test_candidate_selection.R` - New subprocess-contract test (ranking, dedup, no-filter, long-format, shim, NA-fill, empty no_mapping).

## Decisions Made
- **confirmation_status vocabulary = `pass` / `below_threshold`** (per-candidate `reads>minRead && cov>minCov`). Keeps the legacy major-pass behaviour observable; Phase 8 redefines the vocabulary against dominance scoring (Open Q2).
- **Distinct glob for the long-format CSV** (`*.candidates.csv`) to avoid colliding with the legacy `*.parsefirstmapping.csv`.
- **`as.integer` coercion with NA/<1 → default 2** for `n_candidates`, hardening T-06-02 at the R boundary in addition to the nf-schema `minimum: 1`.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 2 - Missing Critical] Guarded empty/X1-less depth frame (T-06-01 DoS mitigation)**
- **Found during:** Task 3 (running the Task-2 RED test to GREEN — the `empty` no-mapping case)
- **Issue:** When the depth file is empty (no reads mapped), `read_tsv(depth, col_names=FALSE)` returns a 0-column tibble with no `X1`, so `group_by(X1)` aborts BEFORE the candidate logic — the script crashed and never reached its `gate_flag="no_mapping"` default. The original committed script had this same latent crash (the inline comment acknowledged "the reading fails" but never guarded it). The Task-2 test asserts no crash + `no_mapping` default, which is exactly the T-06-01 `mitigate` disposition in the plan's threat register.
- **Fix:** Read depth into `depth_raw`; when it is 0-row or lacks `X1`, build a well-typed 0-row `cov` tibble so the no-mapping branch reaches the `gate_flag="no_mapping"` default and still writes both output CSVs.
- **Files modified:** bin/summarize_mapping_to_all_references.R
- **Verification:** Direct empty-input run now exits 0 and writes `EMPTY.parsefirstmapping.csv` with `gate_flag=no_mapping`; the test `empty` case PASSes.
- **Committed in:** 43904de (Task 3 commit)

**2. [Rule 1 - Bug] Fixed test-fixture `writeLines(NULL)` crash on empty refs**
- **Found during:** Task 3 (GREEN run)
- **Issue:** In `bin/tests/test_candidate_selection.R` (committed in Task 2), the `empty` case passes `refs = list()`, making `fa_lines <- unlist(lapply(...))` return `NULL`; `writeLines(NULL, ...)` errors ("can only write character objects"), aborting the test before the empty assertions ran.
- **Fix:** Coerce `fa_lines` to `character(0)` when NULL (mirroring the existing `depth_lines` guard) before `writeLines`.
- **Files modified:** bin/tests/test_candidate_selection.R
- **Verification:** The `empty` case now reaches and passes its `gate_flag=no_mapping` assertion.
- **Committed in:** 43904de (Task 3 commit)

**3. [Rule 3 - Blocking] Reworded an `is_valid_minor` mention out of a comment**
- **Found during:** Task 3 (acceptance grep `! grep -q 'is_valid_minor'`)
- **Issue:** The function was genuinely deleted, but a code comment still contained the literal token `is_valid_minor()`, which would trip the verify-block grep that asserts the symbol is absent.
- **Fix:** Reworded the comment to "The old validity-filtering function is DELETED here" without the literal token.
- **Files modified:** bin/summarize_mapping_to_all_references.R
- **Verification:** `! grep -q 'is_valid_minor'` now passes; only intent is documented.
- **Committed in:** 43904de (Task 3 commit)

---

**Total deviations:** 3 auto-fixed (1 missing-critical T-06-01 guard, 1 test bug, 1 blocking grep-comment)
**Impact on plan:** The empty-depth guard (Rule 2) is the substantive one — it implements the T-06-01 mitigation the plan's threat register requires and that the Task-2 test asserts; without it the no_mapping path was unreachable. The other two are mechanical. No scope creep, no architectural change.

## Issues Encountered
- A usage-limit pause interrupted the prior executor mid-Task-3 (script edit uncommitted). On resume the uncommitted diff already satisfied the Task 3 spec; verification surfaced the empty-input crash (above), which was the only real gap. Tasks 1-2 were confirmed present and not redone.

## Threat Flags
None - no new security surface introduced beyond the plan's threat register. T-06-01 (DoS on empty frame) is now mitigated; T-06-02 hardened at the R boundary; T-06-03 unchanged (controlled reference names).

## Next Phase Readiness
- The long-format `*.candidates.csv` schema is locked; Plan 02 (module emit) and Plan 03 (workflow fan-out) can build on it.
- The legacy shim (`*.parsefirstmapping.csv` + `_major.fa`/`_minor.fa`) is fully reconstructed, so `bin/summarize.R` and the existing Nextflow routing are undisturbed.
- Concern carried forward: `MAJOR_MAPPING`/`MINOR_MAPPING` aliases in `workflows/hcvtyper.nf` feed HCVGLUE downstream — confirm the per-candidate generalization (Plan 03) does not break the paused-v2.0 Phase-5 wiring.

## Self-Check: PASSED

All 5 plan files + SUMMARY.md exist on disk; all 3 task commits (52fda8d, ea71811, 43904de) present in git history.

---
*Phase: 06-neutral-candidate-selection*
*Completed: 2026-06-12*
