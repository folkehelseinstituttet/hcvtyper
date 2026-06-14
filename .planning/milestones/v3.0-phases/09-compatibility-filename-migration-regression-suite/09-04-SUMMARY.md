---
phase: 09-compatibility-filename-migration-regression-suite
plan: 04
subsystem: regression-suite
tags: [compat-01, compat-02, compat-03, compat-04, test-01, regression, summarize]
requires:
  - "09-01 producer-side .cand{rank}. slot rename (config + parsefirstmapping/blastparse/consensus)"
  - "09-02 N-FASTA workflow fan-out"
  - "09-03 summarize.R candidate_rank-by-join consumer cutover"
provides:
  - "bin/tests/test_compat.R — regression suite covering COMPAT-01/02/03/04, auto-discovered by run_all.sh"
  - "human-verified bin/tests/fixtures/compat_golden.csv anchor for the two COMPAT-01 cases"
  - "executable guard against RESEARCH Pitfall 1 (silent Major_*/Minor_* emptying under the cand-slot rename)"
affects:
  - "the r-regression CI job (runs run_all.sh in the pinned container) now also exercises summarize.R end-to-end"
tech-stack:
  added: []
  patterns:
    - "subprocess test driving the REAL summarize.R via system2 on a staged tempdir input tree"
    - "dual-harness test: subprocess Summary.csv signal + sourced-helper exact role/role_reason assertions"
    - "shQuote(\"\") for empty positional slots so system2 does not drop and shift later args"
key-files:
  created:
    - bin/tests/test_compat.R
  modified:
    - bin/tests/fixtures/compat_golden.csv
decisions:
  - "Reused the COMPAT-01 co-infection (1a/1b) subprocess run as the COMPAT-04 end-to-end signal — it already exercises the production score_candidates -> classify_roles call site, no second 1a/1b run needed"
  - "COMPAT-04 helper-level role_reason assertions source the REAL classify_roles.R rather than re-run summarize.R, because the wide Summary.csv does not surface per-candidate role_reason"
  - "[Rule 1] empty positional args [4-8] passed as shQuote(\"\") — a bare \"\" is dropped by system2, silently shifting minRead/minCov off args[9]/[10]"
metrics:
  duration: ~40min
  completed: 2026-06-14
---

# Phase 9 Plan 04: Compatibility + Filename-Migration Regression Suite Summary

Added `bin/tests/test_compat.R` — the validation half of Phase 9 — a single auto-discovered regression file that drives the REAL `bin/summarize.R` end-to-end on synthetic `.cand{rank}.` fixtures and asserts the four COMPAT requirements (golden strain-call reproduction, cand-slot lockstep, legacy+role column co-presence, and the 1a/1b-allowed / 2k1b-suppressed exceptions), green in the pinned Seqera container with no CI YAML change.

## What Was Built

**Task 1 — compat_golden.csv (commit 01f7c0b, human-APPROVED):**
- Two-row human-verified golden anchor (`MONO` monoinfection, `COINF` 1a/1b co-infection) holding the five core strain-call columns. Reference accessions `1a_M62321` / `1b_D90208` confirmed real genotype-1 entries; `Major_/Minor_genotype_mapping` confirmed as the subtype token (`1a`/`1b`) matching `summarize.R`'s `separate(candidate_ref, sep="_")` actual output. Locked as-is by the human "Yes to both".

**Task 2 — COMPAT-01/02/03 subprocess cases (commit ae4fce5):**
- New `bin/tests/test_compat.R` with the test_candidate_selection.R entrypoint boilerplate (`args_all`/`file_arg`/`this_dir`/`bin_dir`), `fail()`/`ok()` helpers, and a `run_summarize()` harness modelled on `run_select`.
- `run_summarize()` stages a one-sample summarize.R input tree in a tempdir — `trimmed/`, `kraken_classified/`, `id/` (each with one file so their `1:length(files)` loops do not crash on `1:0`), plus `parsefirst_mapping/` (`.candidates.csv` 8-col long + `.parsefirstmapping.csv` 10-col shim), `stats_withdup/`, `stats_markdup/`, `depth/` with `.cand1.`/`.cand2.` slot fixtures — copies the five cwd-relative sourced helpers in, invokes the REAL `summarize.R` via `system2`, and reads back `Summary.csv`.
- COMPAT-01: compares the five core columns vs `compat_golden.csv` for both MONO and COINF (NA-tolerant minor cells).
- COMPAT-02 smoke: asserts the co-infection `Major_reference`/`Minor_reference` are non-NA and carry NO `cand{rank}` suffix (the RESEARCH Pitfall 1 silent-empty guard).
- COMPAT-03: asserts the legacy `Major_*`/`Minor_*` columns and the new `Major_role_*`/`Minor_role_*`/`overall_sample_call` columns are simultaneously present (validate-only).

**Task 3 — COMPAT-04 end-to-end (commit d117a34):**
- Part A (end-to-end): reuses the COMPAT-01 co-infection (1a/1b) run to assert `overall_sample_call == "co-infection"` with `Minor_role_reference == "1b_D90208"`; adds a 2a-dominant + 2k1b-second `run_summarize()` case asserting the recombinant is suppressed (sample not co-infection, 2k1b ref absent from `Minor_role_reference`).
- Part B (helper): sources the REAL `classify_roles.R` + `genotype_utils.R` and mirrors the test_classify_roles.R sim1/2k1b cases to lock the exact vocabulary — 1b candidate `role == "co-infection"`; 2k1b candidate `role == "background"` && `role_reason == "recombinant_2k1b"`.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Empty positional args dropped by system2 (arg-position shift)**
- **Found during:** Task 2 (first run yielded `overall_sample_call == "indeterminate"` instead of `monoinfection`).
- **Issue:** The plan's prescribed `system2("Rscript", c(..., "", "", "", "", "", minRead, minCov, n_candidates))` passes the five empty denovo-param slots as bare `""`. `system2` DROPS zero-length string arguments, so `minRead`/`minCov`/`n_candidates` landed at args[4]/[5]/[6] instead of [9]/[10]/[11]. `min_targeted_read`/`min_targeted_cov` then parsed as NA, `classify_roles()`'s `reads > minRead & cov > minCov` floor never passed (NA guard), no candidate became dominant, and every sample collapsed to `overall_sample_call == "indeterminate"`. This is a test-harness bug (the production Nextflow ext.args string never has empty positional gaps), not a summarize.R bug.
- **Fix:** Pass each empty slot as `shQuote("")` so the empty token survives as a distinct positional argument; `nchar(args[i]) > 0` is then 0 and summarize.R's index-guarded defaults apply correctly. Verified `commandArgs` receives 11 args with [4-8] empty and [9-11] = 500/30/2.
- **Files modified:** bin/tests/test_compat.R
- **Commit:** ae4fce5

## Verification

- Task 2 verify (`COMPAT_GREEN`): parse-check + `Rscript bin/tests/test_compat.R` prints `ALL PASS` in the pinned `community.wave.seqera.io/library/r-seqinr_r-tidyverse:5358395134867368` container.
- Task 3 verify (`COMPAT4_GREEN`): same, plus the log contains the COMPAT-04 / 2k1b / 1a/1b / co-infection markers.
- Auto-discovery + no regression: `docker run ... bash bin/tests/run_all.sh` → `>>> ALL R TESTS PASSED`, with `test_compat.R` listed among the 10 discovered `test_*.R` files (TEST-01, no CI YAML edit).

## TDD Gate Compliance

Task 3 was marked `tdd="true"`. The behaviour it asserts (the D-12 1a/1b allowance + 2k1b recombinant suppression in `is_valid_minor`/`classify_roles`) was already SHIPPED in Phase 8 and is unit-covered in `test_classify_roles.R`; COMPAT-04 is an additive regression lock at the integration (summarize.R) level, not new behaviour. A RED-then-GREEN cycle is therefore not applicable — there is no production code to implement, only the regression assertion. The block was authored and verified green directly; the commit uses `feat(...)` per the plan's `<files>`/intent. No new strain-call logic was added to any source file.

## Threat Surface

- T-09-04 (silent data loss at the summarize.R candidate_rank join under the cand-slot rename): MITIGATED as planned — the COMPAT-02 smoke assertion fails loudly if `.cand1.`/`.cand2.` filenames empty `Major_*`/`Minor_*`. This test is the executable guard for the highest-risk integration site (RESEARCH Pitfall 1).
- T-09-04-FP (false-green golden baseline): MITIGATED — the Task 1 blocking human-verify checkpoint confirmed the golden values before assertions locked.
- No new security-relevant surface introduced (test-only file; no network/auth/schema changes). No package installs (T-09-SC accept; runs in the already-pinned container).

## Self-Check: PASSED

- `bin/tests/test_compat.R` exists (created this plan).
- `bin/tests/fixtures/compat_golden.csv` exists (Task 1, human-approved).
- Commits present in `git log`: 01f7c0b (Task 1), ae4fce5 (Task 2), d117a34 (Task 3).
