---
phase: 07-per-genotype-assembly-support
plan: 01
subsystem: denovo-evidence
tags: [blast_parse, assembly-support, ASUP-01, nextflow-module, r-utility]
requires:
  - "scaf_top (best-hit-per-scaffold tibble) in bin/blast_parse.R §4"
  - "sc_length / kmer_cov / subtype / pident / length columns already extracted per BLAST hit"
provides:
  - "*.assembly_support.csv per sample (subtype grain, raw four ASUP-01 metrics, single-best-by-sc_length)"
  - "BLASTPARSE emit channel `support` + lockstep stub"
  - "bin/tests/test_assembly_support.R subprocess-contract test"
affects:
  - "bin/summarize.R (Plan 02 consumes *.assembly_support.csv for the genotype-level candidate join)"
tech-stack:
  added: []
  patterns:
    - "Neutral per-subtype roll-up: group_by(subtype) %>% slice_max(sc_length) %>% distinct()"
    - "Typed header-only CSV on empty input (T-07-01 DoS guard)"
    - "Subprocess-contract test invoking the real script via system2 on inline fixtures"
key-files:
  created:
    - "bin/tests/test_assembly_support.R"
  modified:
    - "bin/blast_parse.R"
    - "modules/local/blastparse/main.nf"
decisions:
  - "Emit subtype grain + raw subtype; genotype key derived in summarize.R (D-02) — no genotype_utils.R staged into BLASTPARSE"
  - "Raw metrics only, no denovo_min_* floor (D-01) — substantiality verdict is Phase 8"
  - "Single best contig by sc_length backs all four metrics (D-03); distinct() de-dup is load-bearing"
  - "New dedicated *.assembly_support.csv (non-optional) rather than folding into *_blast_out.csv — cleaner join target for Plan 02"
metrics:
  duration: "~12 min"
  completed: 2026-06-13
---

# Phase 7 Plan 01: Per-Genotype Assembly Support (emit side) Summary

Dominance-neutral per-subtype assembly-support roll-up: for every subtype in the de novo contigs, `bin/blast_parse.R` now emits one `*.assembly_support.csv` row carrying the single best contig (by `sc_length`, D-03) and that contig's four ASUP-01 metrics — full contig length, BLAST % identity, BLAST alignment length, and k-mer coverage — computed with no reference to the mapping/candidate result. The BLASTPARSE module declares and stubs the new output, and a subprocess-contract test pins the roll-up behaviour.

## What Was Built

**Task 1 — §4b roll-up in `bin/blast_parse.R` (commit 88e0884)**
Inserted a new section after §4 `scaf_top` that groups `scaf_top` by `subtype`, picks the single best contig per subtype via `slice_max(sc_length, n = 1, with_ties = FALSE)`, and carries that winning row's `pident`, BLAST `length`, `sc_length`, and `kmer_cov` into a six-column tibble (`sample, subtype, best_contig_length, best_contig_pident, best_contig_aln_length, best_contig_kmer_cov`). Raw metrics only — no `denovo_min_*` floor (D-01). Raw subtype token carried — no genotype derivation here (D-02). Empty/skip-assembly input writes a typed header-only CSV and exits 0 (T-07-01 DoS guard). Legacy §6/§7 (`major_ref`/`minor_ref`, `*.blastparse.csv`, `_blast_out.csv`) left byte-unchanged (D-04).

**Task 2 — BLASTPARSE declare + stub (commit f9489da)**
Added the non-optional output `tuple val(meta), path("*assembly_support.csv"), emit: support` and a lockstep `stub:` printf emitting the real six-column header so a `-stub-run` yields a parseable file. No `genotype_utils.R` input staged (genotype key is a summarize.R concern, D-02).

**Task 3 — subprocess-contract test (commit 3287a6c)**
`bin/tests/test_assembly_support.R`, modelled on `test_candidate_selection.R`: builds inline outfmt6 blast_out / contigs FASTA / references FASTA fixtures in per-case tempdirs, invokes the real `blast_parse.R` via `system2`, and asserts the emitted CSV. Case A (single-best-by-`sc_length`: 3a row carries the 2949 contig's four metrics, not the 300 one; exactly 2 subtype rows). Case B (multi-hit de-dup: one contig with two same-ref hits → one row, length not duplicated). Case C (empty input → header-only CSV, exit 0, six-column contract). Auto-discovered by `bin/tests/run_all.sh`.

## Verification

- `Rscript -e 'parse("bin/blast_parse.R")'` → exit 0 (PARSE_OK).
- `Rscript bin/tests/test_assembly_support.R` → ALL PASS (exit 0).
- `bash bin/tests/run_all.sh` → ALL R TESTS PASSED (no regression; new test auto-discovered).
- `grep` confirms BLASTPARSE declares `emit: support` and stubs `assembly_support.csv` (count 2).
- `git diff HEAD~3 HEAD -- bin/blast_parse.R` over §6/§7 legacy lines → no changes (D-04 preserved, byte-unchanged).
- TDD load-bearing check: the test was confirmed RED against the pre-Task-1 `blast_parse.R` (old script writes no `*.assembly_support.csv`), then GREEN after Task 1.

## Deviations from Plan

None - plan executed exactly as written.

The plan tagged Task 1 and Task 3 as `tdd="true"`. Task 1 (implementation) was committed before Task 3 (test) per the plan's task ordering; TDD discipline was preserved by separately confirming the Task 3 test is RED against the pre-Task-1 script and GREEN after, so the test is load-bearing rather than vacuous.

## TDD Gate Compliance

Plan type is `execute` (not plan-level `tdd`), so the per-plan RED/GREEN/REFACTOR gate sequence does not apply. The two `tdd="true"` tasks are covered: a `feat(...)` implementation commit (88e0884) and a `test(...)` commit (3287a6c) both exist, and the test was verified to fail against the pre-implementation script (RED) and pass after (GREEN). No REFACTOR commit was needed.

## Known Stubs

None. The BLASTPARSE `stub:` printf is an intentional Nextflow stub-run artifact (every declared non-optional output needs a deterministic stub), not a data stub — the real `script:` path writes the populated CSV via `blast_parse.R`.

## Notes for Plan 02

- `*.assembly_support.csv` is emitted at **subtype grain** with a **raw `subtype`** token; `bin/summarize.R` derives the genotype key via `genotype_from_subtype()` and collapses to `denovo_match_level` at join time (D-02).
- The four metrics on each row all come from one coherent best-by-`sc_length` contig (D-03); the genotype-collapse step in summarize.R must re-apply the same single-best discipline when rolling subtype rows up to genotype.
- No threshold floor is applied here (D-01); the substantiality verdict stays a Phase 8 concern.
- The legacy `*.blastparse.csv` / `denovo_minor_*` path runs in parallel this phase (D-04) — both appear in `Summary.csv` until Phase 8/9.

## Self-Check: PASSED

- Created files verified on disk: `bin/tests/test_assembly_support.R`, `07-01-SUMMARY.md`.
- Commits verified in git log: 88e0884 (Task 1), f9489da (Task 2), 3287a6c (Task 3).
