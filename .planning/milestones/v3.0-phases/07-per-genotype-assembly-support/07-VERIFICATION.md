---
phase: 07-per-genotype-assembly-support
verified: 2026-06-13T12:00:00Z
status: passed
score: 4/4 must-haves verified
overrides_applied: 0
re_verification:
  previous_status: gaps_found
  previous_score: 3/4
  gaps_closed:
    - "Assembly support is joined to candidates at genotype level using a parameterized match level (default genotype) and the result reaches Summary.csv in a real pipeline run"
  gaps_remaining: []
  regressions: []
---

# Phase 7: Per-Genotype Assembly Support — Verification Report

**Phase Goal:** De novo/BLAST evidence stops being an output-only or minor-only QC artefact and becomes an independent, per-genotype "assembly support" signal computed without reference to which candidate is dominant — ready to corroborate any candidate at genotype level.

**Verified:** 2026-06-13T12:00:00Z
**Status:** passed
**Re-verification:** Yes — after gap closure (commit 4fe3c85)

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | For every genotype seen in the de novo contigs, assembly support is summarized as best contig length, BLAST identity, BLAST alignment length, and k-mer coverage, computed independently of mapping/candidate result | VERIFIED | `blast_parse.R` lines 166-199 (section 4b): `group_by(subtype) %>% slice_max(sc_length, n=1, with_ties=FALSE) %>% distinct()` with all four metrics; no threshold floor applied; raw subtype emitted; empty guard writes header-only CSV on zero hits |
| 2 | Assembly support is joined to candidates at genotype level using a parameterized match level (default genotype, not subtype), so 3b is corroborated by 3a | VERIFIED | `join_assembly_support()` in `bin/assembly_support_join.R` implements the correct genotype-level collapse and `left_join`. Fix confirmed in commit 4fe3c85: `path(assembly_support_join)` added as input 18 in `modules/local/summarize/main.nf` line 33; `file("${projectDir}/bin/assembly_support_join.R")` passed as 18th argument at `workflows/hcvtyper.nf` line 531; `input[17]` added to both test cases in `modules/local/summarize/tests/main.nf.test` (lines 40, 79). The `source("assembly_support_join.R")` at `summarize.R:16` will now succeed because the file is staged into the task workdir. Module input arity (18) matches call arity (18). |
| 3 | A candidate with no corresponding genotype-level assembly evidence resolves to an explicit "no support" state rather than silently dropping the row | VERIFIED | `join_assembly_support()` uses candidates as the LEFT side; `assembly_support = if_else(is.na(assembly_support_best_contig_length), "none", "supported")`; empty support returns all "none" with NA metrics; Test3, Test3b, Test3c all pass. |
| 4 | With default flags and candidate count 2, the assembly-support join reproduces the de novo evidence currently attached to the minor slot for the regression fixtures | VERIFIED | Test4 in `test_assembly_support_join.R` asserts `cand_2 assembly_support_best_contig_length == denovo_minor_contig_length` (both 2949 for the ERR1810453 anchor); passes. |

**Score:** 4/4 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `bin/blast_parse.R` | Per-subtype support roll-up emitting `*.assembly_support.csv` | VERIFIED | Section 4b at lines 158-199; correct group_by/slice_max/distinct/transmute/write_csv pattern; no threshold floor; legacy §6/§7 unchanged |
| `modules/local/blastparse/main.nf` | BLASTPARSE declares and stubs `*.assembly_support.csv` | VERIFIED | Line 25: `tuple val(meta), path("*assembly_support.csv"), emit: support` (non-optional); stub line 63 writes correct 6-column header |
| `bin/assembly_support_join.R` | Sourceable `join_assembly_support()` helper | VERIFIED | Function exists, implements correct genotype/subtype match via `genotype_from_subtype()`, CR-01 `as.character()` coercion present (lines 90, 110), CR-02 typed empty support frame present, WR-04 `stopifnot(match_level %in% ...)` present |
| `bin/summarize.R` | Sources helper, reads both CSVs, calls join, carries columns | VERIFIED | Lines 16, 519-579: sources helper; reads `*.candidates.csv` with pinned col_types (CR-01/CR-02 fix); reads `*.assembly_support.csv` with pinned col_types; calls `join_assembly_support()`; WR-01/WR-02 fixed schema via `expected_wide_cols` |
| `modules/local/summarize/main.nf` | Declares `path(assembly_support_join)` as 18th input | VERIFIED | Line 33: `path(assembly_support_join)` present after `path(denovo_layer)` — the file is now staged into the SUMMARIZE task workdir on a real run |
| `workflows/hcvtyper.nf` | Stages `*.candidates.csv` + `*.assembly_support.csv` + `assembly_support_join.R` into SUMMARIZE | VERIFIED | Line 483 stages `PARSEFIRSTMAPPING.out.candidates` into `parsefirst_mapping/`; line 501 mixes `BLASTPARSE.out.support` into `ch_denovo`; line 531 passes `file("${projectDir}/bin/assembly_support_join.R")` as the 18th argument to the SUMMARIZE() call |
| `modules/local/summarize/tests/main.nf.test` | nf-test extended with `input[17]` for the helper | VERIFIED | Lines 40 and 79: `input[17] = file("${projectDir}/bin/assembly_support_join.R", checkIfExists: true)` in both the comprehensive and stub test cases |
| `bin/tests/test_assembly_support.R` | Subprocess-contract test for per-subtype roll-up | VERIFIED | 3 cases (single-best-by-length, multi-hit de-dup, empty no-abort); all pass |
| `bin/tests/test_assembly_support_join.R` | Function test incl. CR-01/CR-02 locking and criterion #4 | VERIFIED | 9 tests (Tests 1-6, 3b, 3c, 5b) covering genotype match, subtype non-match, no-row-loss, criterion #4 reproduction, CSV round-trip (CR-01/CR-02), numeric key coercion, match_level validation; all pass |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `bin/blast_parse.R` | `scaf_top` | `group_by(subtype) %>% slice_max(sc_length` | WIRED | Lines 168-172: exactly the pattern declared in PLAN |
| `modules/local/blastparse/main.nf` | `*.assembly_support.csv` | `emit: support` + stub printf | WIRED | Lines 25, 63: output declared non-optional; stub header matches real header |
| `bin/summarize.R` | `bin/assembly_support_join.R` | `source() + join_assembly_support() call` | WIRED | `source("assembly_support_join.R")` present at line 16; `join_assembly_support(candidates_long, support_df, denovo_match_level)` called at line 579; file is now staged via `path(assembly_support_join)` in the module input block |
| `workflows/hcvtyper.nf` | SUMMARIZE staging | `PARSEFIRSTMAPPING.out.candidates + BLASTPARSE.out.support collected + assembly_support_join.R file()` | WIRED | All three channels/file arguments staged; `assembly_support_join.R` passed at line 531 |
| `bin/assembly_support_join.R` | `genotype_from_subtype()` | match-level key on both candidate and support sides | WIRED | Line 111: `if (match_level == "subtype") subtype else genotype_from_subtype(subtype)` |

### Data-Flow Trace (Level 4)

| Artifact | Data Variable | Source | Produces Real Data | Status |
|----------|---------------|--------|--------------------|--------|
| `bin/blast_parse.R §4b` | `support_tbl` | `scaf_top` (BLAST hits) | Yes — `group_by/slice_max/transmute` from real BLAST output | FLOWING |
| `bin/summarize.R` | `candidate_support_wide` | `join_assembly_support(candidates_long, support_df)` | Yes — helper now staged; data flows from BLASTPARSE assembly_support.csv through the join into Summary.csv | FLOWING |

### Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| Per-subtype roll-up emits correct CSV | `Rscript bin/tests/test_assembly_support.R` | ALL PASS (3/3 cases) | PASS |
| Genotype-level join helper function | `Rscript bin/tests/test_assembly_support_join.R` | ALL PASS (9 tests) | PASS |
| Full R test suite (7 test files) | `bash bin/tests/run_all.sh` | ALL R TESTS PASSED | PASS |
| blast_parse.R parses without error | `Rscript -e 'invisible(parse("bin/blast_parse.R"))'` | OK | PASS |
| assembly_support_join.R defines function | `Rscript -e 'source("bin/genotype_utils.R"); source("bin/assembly_support_join.R"); stopifnot(is.function(join_assembly_support))'` | OK | PASS |

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|-------------|-------------|--------|----------|
| ASUP-01 | 07-01-PLAN.md | De novo/BLAST evidence expressed as per-genotype assembly support (best contig length, BLAST identity/length, k-mer coverage), computed independently of mapping | SATISFIED | `blast_parse.R` §4b: 4-metric roll-up per subtype, dominance-neutral, no threshold floor |
| ASUP-02 | 07-02-PLAN.md | Assembly support joined to candidates at genotype level (parameterized match level), not subtype | SATISFIED | `assembly_support_join.R` implements the join; `summarize.R` calls it; staging gap closed by commit 4fe3c85 — all three artifacts (module input, workflow call, nf-test) now consistent |

REQUIREMENTS.md Traceability table marks both ASUP-01 and ASUP-02 as "Complete" (lines 76-77). Both are confirmed correct.

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| None found | — | No TBD/FIXME/XXX markers in phase-modified files | — | — |

The SUMMARIZE module stub (`modules/local/summarize/main.nf` lines 70-77) does not include the new Phase-7 `cand_*_assembly_support*` columns in the stub Summary.csv header or data row. This is a pre-existing stub-vs-reality gap for all new columns, consistent with how the stub was handled in earlier phases; it does not affect the functional path.

The summarize module nf-test `.snap` regeneration is intentionally deferred because the host disk is near-full and docker nf-test runs risk ENOSPC. Summary.csv now carries the `cand_*_assembly_support*` schema, so the snapshot will need a `--update-snapshot` docker run later. This is a CI/test-infra follow-up, NOT a production-path gap.

### Human Verification Required

None — all verification is automatable for this phase.

## CR-01/CR-02 Fix Verification

The code review (07-REVIEW.md) identified two default-path blockers fixed in commits `ac30a5d` and `bd2461e`. Both fixes are confirmed present after the re-verification commit:

**CR-01 (numeric/character join-key mismatch):** Fixed at two layers:
1. `bin/assembly_support_join.R` lines 90, 110: `as.character()` coerces both sides of `.match_key` before the join.
2. `bin/summarize.R` lines 519-528: `col_types = cols(candidate_genotype = col_character(), ...)` pins the type at read time.

**CR-02 (mixed header-only/populated CSV bind):** Fixed:
1. `bin/summarize.R` lines 519-528: pinned `col_types` on `candidates_long` read.
2. `bin/summarize.R` lines 553-561: pinned `col_types` on `support_df` read.

**Strengthened test (IN-01):** Test5 and Test5b in `test_assembly_support_join.R` write candidates/support to real CSVs, read them back with the same `col_types` as `summarize.R`, and assert the join succeeds with a numeric-inferred `candidate_genotype` — locking in both fixes against regression. Both pass.

## Gap Closure Confirmation

The single blocker from the initial verification has been resolved:

**Closed:** `assembly_support_join.R` was not staged into the SUMMARIZE Nextflow task workdir. Fixed in commit `4fe3c85` (message: "fix(07): stage assembly_support_join.R into SUMMARIZE workdir (V-01)"):
- `modules/local/summarize/main.nf` line 33: `path(assembly_support_join)` added as the 18th declared input.
- `workflows/hcvtyper.nf` line 531: `file("${projectDir}/bin/assembly_support_join.R")` added as the 18th argument in the `SUMMARIZE()` call.
- `modules/local/summarize/tests/main.nf.test` lines 40, 79: `input[17] = file("${projectDir}/bin/assembly_support_join.R", checkIfExists: true)` added to both test cases.

The module input arity (18) now matches the workflow call arity (18). The `source("assembly_support_join.R")` at `summarize.R:16` will succeed on a real pipeline run.

---

_Verified: 2026-06-13T12:00:00Z_
_Verifier: Claude (gsd-verifier)_
