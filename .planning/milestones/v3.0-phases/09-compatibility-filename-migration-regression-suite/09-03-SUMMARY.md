---
phase: 09-compatibility-filename-migration-regression-suite
plan: 03
subsystem: summarize-r-consumer
tags: [compat-02, summarize, candidate-rank, filename-migration]
requires:
  - "09-01 producer-side .cand{rank}. slot rename (config + parsefirstmapping/blastparse/consensus)"
provides:
  - "summarize.R recovers candidate rank by joining candidate_ref to candidate_rank (no filename position-3 parse)"
  - "Major_*/Minor_* logic driven by candidate_rank == 1 / == 2"
  - "byte-identical cleaned Major_reference/Minor_reference join keys (no slot suffix)"
  - "cv_by_ref + consensus-distance blocks migrated to _cand{rank} parsing"
affects:
  - "09-04 regression suite (test_compat.R) exercises this consumer-side parse"
tech-stack:
  added: []
  patterns:
    - "candidate_rank-by-join replaces fragile filename-slot string extraction"
    - "byte-identical cleaned-ref join key preserved across slot rename"
key-files:
  created: []
  modified:
    - bin/summarize.R
decisions:
  - "Hoisted candidates_long load (+ candidate_rank_lookup) before the three stats loops so each loop can join on candidate_rank"
  - "first_mapping rows kept via reference == \"first_mapping\" guard (no candidate_rank assigned)"
  - "[Rule 1] migrated a fifth, plan-unlisted slot-coupled site (consensus distance block) that Plan 01's CONSENSUS_DISTANCE prefix rename had silently broken"
metrics:
  duration: ~25min
  completed: 2026-06-14
---

# Phase 9 Plan 03: summarize.R candidate_rank Consumer Cutover Summary

Refactored `bin/summarize.R` so it recovers each stats file's candidate rank by joining the cleaned reference token to `candidate_rank` from the candidates CSV, instead of parsing a `.major.`/`.minor.` slot from filename position 3 — completing the consumer half of the COMPAT-02 lockstep with Plan 01's producer-side `.cand{rank}.` rename.

## What Was Built

**Task 1 — candidate_rank-by-join across the three stats loops (commit bbe124e):**
- Hoisted the `candidates_long` load block (and a new `candidate_rank_lookup` of distinct `(sampleName, candidate_ref, candidate_rank)`) from its original position near the Phase-7 assembly-support join to *before* the first stats loop, since the load depends only on `path_3` (set near the top). The Phase-7 join below reuses the same hoisted frame.
- Block 1 (withdup stats), Block 2 (nodup stats), Block 3 (coverage): removed the `str_split(basename(...), "\\.")[[1]][3]` `first_major_minor` extraction; after each loop, derive `candidate_ref = str_remove(reference, "_cand[0-9]+$")` and `left_join(candidate_rank_lookup, by = c("sampleName", "candidate_ref"))`.
- Replaced every `case_when(first_major_minor == "major" ~ X)` with `case_when(candidate_rank == 1 ~ X)` and `"minor"` with `candidate_rank == 2`.
- `Major_reference`/`Minor_reference` now populated from the cleaned `candidate_ref` (no slot suffix), keeping them byte-identical to the legacy value so the downstream join keys at L366 (`full_join`) and the coverage `left_join` still match.
- `first_mapping` row in Block 2 preserved via `reference == "first_mapping"` guard (these rows have no `candidate_rank`).
- Updated the per-loop `select(-first_major_minor)` collapses to drop `candidate_rank`/`candidate_ref` so the per-sample output schema is unchanged.

**Task 2 — cv_by_ref slot strip migration (commit b49b1cf):**
- Changed the `cv_by_ref` derivation from `str_remove(reference, "_(major|minor)$")` to `str_remove(reference, "_cand[0-9]+$")`, identical to the stats-loop regex so per-candidate `cv_evenness` still joins to the Phase-6 `candidate_ref`.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Migrated a fifth, plan-unlisted slot-coupled site (consensus distance block)**
- **Found during:** Task 1 (the `NO_SLOT_PARSE` verification surfaced remaining `first_major_minor == "major"/"minor"` at L1046-1049).
- **Issue:** The plan enumerated four slot-coupled sites (three stats loops + cv_by_ref) but `bin/summarize.R` has a fifth: the consensus-distance block (`df_consensus_distance` / `df_distance_wide`). It parsed `major_minor <- fname_parts[2]` and drove `case_when(first_major_minor == "major"/"minor" ~ ...)`. Plan 01 had already renamed the `CONSENSUS_DISTANCE` `ext.prefix` (conf/modules_hcv.config L263) from `.major.`/`.minor.` to `.cand{rank}.`, so this block would have silently emptied the `Major_*/Minor_consensus_*` columns (RESEARCH Pitfall 1, the documented silent-NA failure mode). This is squarely in COMPAT-02 scope (the consumer-side filename migration in summarize.R).
- **Fix:** Parse `slot <- fname_parts[2]`, derive `cand_rank <- as.integer(str_remove(slot, "^cand"))`, store as `candidate_rank`, and drive the wide pivot off `candidate_rank == 1 / == 2`. The internal `sample`-column validation regex was likewise changed from `str_extract(sample_col, "major|minor")` to `str_extract(sample_col, "cand[0-9]+")` (the iVar consensus header carries the `${meta.id}.cand{rank}` prefix from IVAR_CONSENSUS).
- **Files modified:** bin/summarize.R
- **Commit:** bbe124e (folded into Task 1 since it shares the same file and the same verification gate).

## Verification

- Task 1 automated: `NO_SLOT_PARSE`, `HAS_RANK_LOGIC`, `NO_POS3_PARSE` all pass.
- Task 2 automated: `HAS_CAND_STRIP`, `NO_LEGACY_STRIP` all pass.
- Remaining `major`/`minor` non-comment strings are out-of-scope first-mapping percent columns (L175-178, upstream `parsefirstmapping` data) and the legacy COMPAT-03 output column names (`Reads_*_mapped_major`, etc.) that must stay — none are reference-cleanup slot strips.
- Syntax: `Rscript -e 'parse("bin/summarize.R")'` → `PARSE_OK` in the pinned `r-seqinr_r-tidyverse` container.
- Regression: full `bin/tests/run_all.sh` → `ALL R TESTS PASSED`, including `Block4 (D-05/D-06): flag-OFF reproduces committed legacy golden baseline` — confirms the byte-identical cleaned-ref join keys did not regress the COMPAT golden output.

## Threat Surface

T-09-03 (silent data loss at the `Major_reference`/`Minor_reference` join keys) mitigated as planned: cleaned-ref values use `_cand[0-9]+$` strip and stay byte-identical. The discovered consensus-distance site was the concrete instance of this threat the plan's threat model implied; it is now closed. No new security-relevant surface introduced.

## Self-Check: PASSED

- bin/summarize.R exists and was modified (commits bbe124e, b49b1cf present in `git log`).
