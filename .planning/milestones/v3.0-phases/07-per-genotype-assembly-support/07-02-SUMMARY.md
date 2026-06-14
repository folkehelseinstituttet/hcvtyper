---
phase: 07-per-genotype-assembly-support
plan: 02
subsystem: denovo-evidence
tags: [summarize, assembly-support, ASUP-02, genotype-join, r-utility, nextflow-staging]
requires:
  - "*.assembly_support.csv per sample (subtype grain, six Plan-01 columns) from BLASTPARSE.out.support"
  - "*.candidates.csv long-format candidate table (8 Phase-6 columns incl. precomputed candidate_genotype/candidate_subtype) from PARSEFIRSTMAPPING.out.candidates"
  - "genotype_from_subtype() (bin/genotype_utils.R), 2k1b-aware"
provides:
  - "bin/assembly_support_join.R::join_assembly_support(candidates_df, support_df, match_level) — genotype-level candidate<->support join + NA-fill"
  - "Per-candidate assembly-support columns in Summary.csv (wide cand_<rank>_assembly_support* slots)"
  - "bin/tests/test_assembly_support_join.R function-level test (criteria #2/#3/#4)"
affects:
  - "bin/summarize.R (reads both CSVs, joins, pivots to wide, left_joins onto final)"
  - "workflows/hcvtyper.nf (stages candidates + assembly_support into SUMMARIZE)"
  - "Phase 8 (consumes the per-candidate assembly_support signal for the substantiality verdict)"
tech-stack:
  added: []
  patterns:
    - "Sourceable pure join helper mirroring denovo_layer.R (no commandArgs, no source())"
    - "PLUMB-02 typed-empty-tibble guards on the new candidate + support reads"
    - "left_join candidates as LEFT side -> explicit assembly_support='none' + NA, no row loss"
    - "Collapse subtype->genotype via slice_max(best_contig_length, with_ties=FALSE) (D-03)"
    - "pivot_wider per-candidate support to cand_<rank>_* slots before joining one-row/sample final"
key-files:
  created:
    - "bin/assembly_support_join.R"
    - "bin/tests/test_assembly_support_join.R"
  modified:
    - "bin/summarize.R"
    - "workflows/hcvtyper.nf"
decisions:
  - "Candidate-side match key reuses precomputed candidate_genotype/candidate_subtype, not recomputed (07-PATTERNS line 118)"
  - "Support-side match key derived with the exact denovo_confirm.R line-64 branch on subtype"
  - "Per-candidate support pivoted to wide cand_<rank>_* slots so it left_joins onto the one-row/sample final without exploding rows"
  - "Legacy minor-coupled path (apply_denovo_layer/review_flag/denovo_minor_*) byte-unchanged; intentional column duplication this phase (D-04)"
metrics:
  duration: "~4 min"
  completed: 2026-06-13
---

# Phase 7 Plan 02: Per-Genotype Assembly Support (join side) Summary

Genotype-level join of the Plan-01 per-subtype assembly-support table onto the Phase-6 long-format candidate set. A new sourceable `join_assembly_support()` left-joins candidates (the LEFT side, no row loss) to the per-subtype support at the parameterized `denovo_match_level` (default genotype), collapsing subtypes up to genotype by the single best contig and NA-filling unmatched candidates as an explicit `assembly_support = "none"`. `summarize.R` reads both CSVs with typed-empty guards, calls the helper, pivots the per-candidate support to wide `cand_<rank>_*` slots, and carries them into `Summary.csv` alongside the unchanged legacy `denovo_minor_*` columns (D-04). The workflow stages both CSVs into the existing SUMMARIZE dirs.

## What Was Built

**Task 1 — sourceable join helper `bin/assembly_support_join.R` (commit 196dfdf)**
A pure sourced helper (no `commandArgs`, no `source()`, mirroring `denovo_layer.R`) defining `join_assembly_support(candidates_df, support_df, match_level = "genotype")`. Candidate-side match key reuses the precomputed `candidate_genotype` / `candidate_subtype` columns; support-side key is derived with the exact `denovo_confirm.R` line-64 branch `if (match_level == "subtype") subtype else genotype_from_subtype(subtype)`. Support rows are collapsed to one row per `(sampleName, match_key)` via `slice_max(best_contig_length, n = 1, with_ties = FALSE)` (single-best-contig, D-03), then `left_join`ed onto candidates so no candidate row is ever dropped. Unmatched candidates get `assembly_support = "none"` + NA metrics. Zero-row candidate or empty support inputs return typed frames, never abort (T-07-03).

**Task 3 — function test `bin/tests/test_assembly_support_join.R` (commit a9668e5)**
Modelled on `test_summarize_denovo.R`, sources the REAL helper and asserts: (1) criterion #2 — a 3b candidate is corroborated by 3a support at `match_level="genotype"`, inheriting the 3a contig length 2949; (2) the subtype-level non-match leaves the 3b candidate at `"none"`; (3) criterion #3 — output row count == input candidate row count, explicit `"none"` + NA metrics for the unmatched candidate; (3b/3c) empty support / empty candidates guards never abort; (4) criterion #4 — on an ERR1810453-style fixture the `cand_2` (minor-slot) `assembly_support_best_contig_length` equals the best different-genotype contig length (2949) that today's `blast_parse.R` §7 minor logic puts in `denovo_minor_contig_length`, and the 300 bp noise 2b hit does NOT win the collapse.

**Task 2 — wire into `summarize.R` + workflow staging (commit d80b3ed)**
`summarize.R`: sources `assembly_support_join.R` after `genotype_utils.R`; reads `*.candidates.csv` from `parsefirst_mapping/` (8-column typed-empty fallback, `sample`->`sampleName`) and `*.assembly_support.csv` from `denovo/` (6-column typed-empty fallback); calls `join_assembly_support(candidates_long, support_df, denovo_match_level)`; pivots the per-candidate support to wide `cand_<rank>_assembly_support*` slots and `left_join`s them onto `final` after the existing `df_denovo` join (samplesheet still anchors the left side). The select reorder ends in `everything()`, so the new columns carry through; `distinct()` still yields one row/sample. The legacy `apply_denovo_layer()` / `coinfection_flag` / `review_flag` / `denovo_minor_*` path is byte-unchanged (D-04). `workflows/hcvtyper.nf`: mixes `PARSEFIRSTMAPPING.out.candidates` into `ch_summarize_first_mapping` (lands in `parsefirst_mapping/`) and `BLASTPARSE.out.support` into `ch_denovo` inside the `!params.skip_assembly` guard (lands in `denovo/`); SUMMARIZE module input arity unchanged.

## Verification

- `Rscript -e 'parse("bin/summarize.R")'` -> exit 0 (PARSE_OK).
- `Rscript bin/tests/test_assembly_support_join.R` -> ALL PASS (criteria #2/#3/#4 + empty guards).
- `bash bin/tests/run_all.sh` -> ALL R TESTS PASSED (7 test files; new test auto-discovered; no regression).
- Task 2 `grep` gate -> WIRED_OK (summarize.R sources helper + reads both CSVs; workflow mixes both channels).
- `git diff bin/summarize.R` over legacy `apply_denovo_layer`/`review_flag`/`denovo_minor_*`/`coinfection_flag` lines -> NO_LEGACY_LINES_REMOVED (D-04 byte-unchanged).
- Integration sanity: the join+pivot produces clean `cand_1_*`/`cand_2_*` slots and `cand_2_assembly_support_best_contig_length == 2949` end-to-end (criterion #4 through the pivot).

## Deviations from Plan

None - plan executed exactly as written.

The plan's task ordering lists Task 1 (helper) before Task 3 (test). Both are `tdd="true"`; TDD discipline is preserved because the test sources the real `join_assembly_support()` and asserts on its output (the test cannot run without the Task-1 helper, so it is load-bearing rather than vacuous). The `feat` implementation commit (196dfdf) and the `test` commit (a9668e5) both exist.

## TDD Gate Compliance

Plan type is `execute` (not plan-level `tdd`), so the per-plan RED/GREEN/REFACTOR gate sequence does not apply. The three `tdd="true"` tasks are covered by `feat` implementation commits (196dfdf, d80b3ed) and a `test` commit (a9668e5); the test exercises the real helper and passes. No REFACTOR commit was needed.

## Known Stubs

None. The new `cand_<rank>_assembly_support*` columns are wired to real data via `join_assembly_support()`; they NA-fill (explicit `"none"`) only when a sample genuinely has no genotype-level de novo support, which is the intended criterion-#3 behaviour, not a stub.

## Notes for Phase 8

- The per-candidate `assembly_support` signal is RAW (length / pident / aln_length / kmer_cov + `"supported"`/`"none"` status) — NO substantiality floor is applied (D-01). Phase 8 (CLASS-02/SCORE) applies the `denovo_min_*` floors as the corroboration verdict.
- Both the legacy `denovo_minor_*` columns and the new `cand_<rank>_assembly_support*` columns appear in `Summary.csv` this phase (intentional duplication, D-04) — Phase 8/9 retires the legacy minor-coupled path.
- The wide pivot is keyed by `candidate_rank`; at the default `n_candidates=2` two-slot shim only `cand_1`/`cand_2` slots appear. Criterion #4 anchors `cand_2` to today's `denovo_minor_contig_length`.

## Self-Check: PASSED

- Created files verified on disk: `bin/assembly_support_join.R`, `bin/tests/test_assembly_support_join.R`, `07-02-SUMMARY.md`.
- Commits verified in git log: 196dfdf (Task 1), a9668e5 (Task 3), d80b3ed (Task 2).
