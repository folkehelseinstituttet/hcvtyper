---
phase: 09-compatibility-filename-migration-regression-suite
plan: 02
subsystem: workflow-fanout-module-test-migration
tags: [COMPAT-02, D-01, D-03, D-04, nextflow, nf-test, candidate_fasta]
requires:
  - "Plan 01 candidate_fasta N-FASTA emit + cand-slot FASTAs (*_cand{rank}.fa)"
provides:
  - "Rank-indexed N-FASTA fan-out in workflows/hcvtyper.nf consuming PARSEFIRSTMAPPING.out.candidate_fasta"
  - "n_candidates > 2 guard lifted (D-01) — N>2 candidates now route end-to-end"
  - "parsefirstmapping nf-test asserting cand-slot filenames against the candidate_fasta emit"
  - "Regenerated parsefirstmapping + blastparse .snap golden files (cand-slot, no major/minor)"
affects:
  - "bin/summarize.R (Plan 03 must parse the .cand{rank}. slot in lockstep — consumer side)"
tech-stack:
  added: []
  patterns:
    - "Rank-indexed FASTA pick: fastas?.find { it.toString().contains(\"_cand${rank}.\") }"
    - "Single candidate_fasta join with remainder: true (replaces two major/minor joins)"
    - "List-normalize collected per-sample FASTA emit before rank lookup"
key-files:
  created: []
  modified:
    - workflows/hcvtyper.nf
    - modules/local/parsefirstmapping/tests/main.nf.test
    - modules/local/parsefirstmapping/tests/main.nf.test.snap
    - modules/local/blastparse/tests/main.nf.test.snap
decisions:
  - "FASTA picked by rank-indexed basename match (_cand${rank}.) over the collected candidate_fasta list, replacing the rank==1?major:minor two-slot branch"
  - "candidate_fasta join uses remainder: true (single optional emit) to preserve monoinfection-sample survival"
  - "Snapshot version/timestamp metadata churn (nf-test 0.9.2 vs 0.9.3) accepted as benign tooling drift; all non-filename md5s byte-identical"
metrics:
  duration: ~7min
  completed: 2026-06-14
  tasks: 3
  files: 4
---

# Phase 9 Plan 02: Workflow Fan-out + Module-Test Migration (consumer side) Summary

Completed the COMPAT-02 producer migration on the workflow-routing + module-test side:
deleted the `n_candidates > 2` guard, generalized the two-slot FASTA pick to a rank-indexed
lookup against the single `PARSEFIRSTMAPPING.out.candidate_fasta` N-FASTA emit, migrated the
parsefirstmapping nf-test assertions to cand-slot filenames, and tool-regenerated both module
`.snap` golden files — all seven module nf-tests green, all four anti-regression invariants intact.

## What Was Built

### Task 1 — N-FASTA fan-out + guard lift (`workflows/hcvtyper.nf`)
- Deleted the entire `if (params.n_candidates > 2) { error ... }` guard block + leading
  comment (D-01; the error message pointed at this very phase).
- Replaced the two legacy joins (`.join(PARSEFIRSTMAPPING.out.major_mapping, remainder: true)`
  + `.join(PARSEFIRSTMAPPING.out.minor_mapping, remainder: true)`) with a single
  `.join(PARSEFIRSTMAPPING.out.candidate_fasta, remainder: true)`.
- The candidate_fasta emit is `tuple(meta, parsefirstmapping_csv, [*_cand{rank}.fa, ...])`.
  Inside the `flatMap`, the collected FASTA list is normalized (bare path / list / null → list),
  then each candidate's FASTA is picked by rank-indexed basename match:
  `def fasta = fastas?.find { it.toString().contains("_cand${rank}.") }` — replacing the hard
  `def fasta = (rank == '1') ? major_fasta : minor_fasta` two-slot branch.
- Updated the explanatory comment blocks to describe the single candidate_fasta emit + the
  rank-indexed lookup.
- **All four anti-regression invariants preserved verbatim** (RESEARCH Pitfall 3 / T-09-02):
  (1) `remainder: true` on the candidate_fasta join; (2) `def rank = new_meta.candidate_rank.toString()`
  (rank stays a String, never `.toInteger()` in the workflow); (3) the final
  `.filter { entry -> entry[0]['confirmation_status'] == 'pass' && entry[1] != null }`;
  (4) the `assert new_meta.id == new_meta.sample` check.

### Task 2 — parsefirstmapping nf-test → cand-slot (`modules/local/parsefirstmapping/tests/main.nf.test`)
- Replaced `major_mapping: process.out.major_mapping` + `minor_mapping: process.out.minor_mapping`
  in all five `snapshot([...]).match()` maps with a single `candidate_fasta: process.out.candidate_fasta`;
  kept `csv`/`candidates`/`versions`.
- Migrated the GATE-03 hardcoded assertion from `...major_mapping.get(0).get(2)...contains("2k1b_AB031663_major.fa")`
  to `...candidate_fasta.get(0).get(2)...contains("2k1b_AB031663_cand1.fa")`.
- Replaced both `minor_mapping.size() == 0` single-candidate assertions (GATE-03 + GATE-04)
  with the equivalent statement on the candidate_fasta tuple's FASTA list (position `.get(2)`):
  asserts `_cand1.fa` is present and `_cand2.fa` is absent.
- Updated stale `*_minor.fa` / `minor_mapping` mentions in surrounding comments for accuracy.
- Left untouched: idxstats/tsv fixture inputs, candidate_rank header + row-count assertions,
  and the GATE-02 `cell()` assertions on the legacy wide CSV.

### Task 3 — regenerate both module snapshots (`*/tests/main.nf.test.snap`)
- Tool-regenerated both `.snap` files via `conda run -n NEXTFLOW nf-test ... --update-snapshot`
  (nf-test 0.9.2; host `nextflow` unusable per RESEARCH).
- **Diff-guarded against Pitfall 4 (T-09-02-SNAP):** the only content changes are the slot
  rename (`*_major.fa`/`*_minor.fa` → `*_cand1.fa`/`*_cand2.fa`, `toy.major.fa`/`toy.minor.fa`
  → `toy.cand1.fa`/`toy.cand2.fa`) plus the snapshot-key rename
  (`major_mapping`/`minor_mapping`/`major_fasta`/`minor_fasta` → `candidate_fasta`).
  **Every non-filename md5 is byte-identical** to the pre-rename snapshot — no contig/depth/
  stats/png/candidates.csv/assembly_support.csv md5 churn. The only other diff is benign tooling
  metadata (nf-test 0.9.2 vs prior 0.9.3, nextflow 24.10.1 vs prior 24.10.5, and timestamps).
- Re-ran both module tests WITHOUT `--update-snapshot`: 7/7 green against the regenerated snapshots.

## Verification

| Check | Result |
|-------|--------|
| `params.n_candidates > 2` guard removed | PASS (GUARD_GONE) |
| Fan-out joins `PARSEFIRSTMAPPING.out.candidate_fasta` | PASS (USES_NEW_EMIT) |
| No `major_mapping`/`minor_mapping` in workflow | PASS (NO_LEGACY_EMIT) |
| `remainder: true` + `candidate_rank.toString()` preserved | PASS (GUARDS_PRESERVED) |
| `confirmation_status == 'pass' && entry[1] != null` filter preserved | PASS (FILTER_OK) |
| `assert new_meta.id == new_meta.sample` preserved | PASS (ASSERT_OK) |
| parsefirstmapping test: no `major_mapping`/`minor_mapping`; uses candidate_fasta; asserts `_cand1.fa` | PASS |
| GATE-03 hardcoded assert targets `2k1b_AB031663_cand1.fa` | PASS |
| Neither `.snap` contains `major.fa`/`minor.fa`/`major_fasta`/`minor_fasta`/`major_mapping`/`minor_mapping` | PASS (0/0) |
| blastparse `.snap` carries cand-slot / candidate_fasta entries | PASS (5 matches) |
| Snap diff confined to slot rename; no non-filename md5 churn (Pitfall 4) | PASS |
| All 7 module nf-tests green without `--update-snapshot` | PASS |

## Deviations from Plan

None — plan executed exactly as written. (Comment-only cleanup of stale `*_minor.fa` /
`minor_mapping` references in test docstrings is within the Task 2 scope of removing legacy
emit names; not a behavioral deviation.)

Note on accepted snapshot churn (not a deviation): the regenerated `.snap` files carry a
tooling-metadata diff (nf-test 0.9.2/nextflow 24.10.1 in the local conda env vs the 0.9.3/24.10.5
recorded in the committed snapshots) plus fresh timestamps. Per Pitfall 4 the load-bearing check
is md5 content — every non-filename md5 is byte-identical, so this is benign environment drift,
not a regression.

## Known Stubs

None introduced.

## Self-Check: PASSED

- FOUND: workflows/hcvtyper.nf
- FOUND: modules/local/parsefirstmapping/tests/main.nf.test
- FOUND: modules/local/parsefirstmapping/tests/main.nf.test.snap
- FOUND: modules/local/blastparse/tests/main.nf.test.snap
- FOUND commit 652a343 (Task 1)
- FOUND commit eea5c20 (Task 2)
- FOUND commit 4106f7e (Task 3)
