---
status: passed
phase: 09-compatibility-filename-migration-regression-suite
source:
  - 09-01-SUMMARY.md
  - 09-02-SUMMARY.md
  - 09-03-SUMMARY.md
  - 09-04-SUMMARY.md
started: 2026-06-14T18:30:00Z
updated: 2026-06-14T18:30:00Z
---

## Current Test

number: 5
name: n_candidates > 2 guard is gone from the workflow
expected: no output from grep
awaiting: complete

## Tests

### 1. COMPAT regression suite passes end-to-end
expected: |
  `bash bin/tests/run_all.sh` in the pinned container prints `>>> ALL R TESTS PASSED`
  with test_compat.R auto-discovered alongside the other test_*.R files.
result: pass

### 2. No `.major.`/`.minor.` slot in ext.prefix closures
expected: |
  grep -c "'major'\|'minor'" conf/modules_hcv.config
  Returns 0 (no legacy slot literals in any ext.prefix closure).
  grep -c "cand\${" conf/modules_hcv.config
  Returns 6 (all six TARGETED_MAPPING closures emit cand{rank}).
result: pass
notes: 1 hit is a comment on L123 describing the old scheme — not functional code. 6 cand${ closures confirmed.

### 3. summarize.R no longer parses filename position 3 for major/minor
expected: |
  grep "first_major_minor" bin/summarize.R
  Returns no output (the position-3 filename parse is gone, replaced by candidate_rank join).
  grep "candidate_rank" bin/summarize.R | wc -l
  Returns a non-zero count showing the rank-join is in place.
result: pass
notes: 0 first_major_minor hits, 51 candidate_rank references confirmed.

### 4. parsefirstmapping nf-tests pass with cand-slot snapshots
expected: |
  conda run -n NEXTFLOW nf-test test \
    modules/local/parsefirstmapping/tests/main.nf.test \
    modules/local/blastparse/tests/main.nf.test \
    --profile test,docker
  All tests PASS; no major.fa or minor.fa references in the .snap files.
result: pass

### 5. n_candidates > 2 guard is gone from the workflow
expected: |
  grep "n_candidates > 2" workflows/hcvtyper.nf
  Returns no output — the guard has been deleted and N>2 candidates now route end-to-end.
result: pass

## Summary

total: 5
passed: 5
issues: 0
pending: 0
skipped: 0

## Gaps

[none yet]
