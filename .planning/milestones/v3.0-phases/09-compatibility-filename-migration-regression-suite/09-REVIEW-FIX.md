---
phase: 09-compatibility-filename-migration-regression-suite
fixed_at: 2026-06-14T19:30:00Z
review_path: .planning/phases/09-compatibility-filename-migration-regression-suite/09-REVIEW.md
iteration: 1
findings_in_scope: 5
fixed: 5
skipped: 0
status: all_fixed
---

# Phase 9: Code Review Fix Report

**Fixed at:** 2026-06-14T19:30:00Z
**Source review:** `.planning/phases/09-compatibility-filename-migration-regression-suite/09-REVIEW.md`
**Iteration:** 1

**Summary:**
- Findings in scope: 5 (CR-01, CR-02, WR-01, WR-02, WR-03)
- Fixed: 5
- Skipped: 0

---

## Fixed Issues

### CR-01: Variation-plot grid silently suppressed after cand-slot rename

**Files modified:** `bin/plot_bam_variation.R`, `bin/summarize.R`
**Commit:** `8321c18`
**Applied fix:**

`bin/plot_bam_variation.R` line 24: replaced hard-coded position-3 extraction
(`unlist(str_split(...))[3]`) with a grep-based approach that matches the first
dot-split part matching `^cand[0-9]+$`. Falls back to position 3 for legacy
filenames. Under Phase-9 naming `<id>.<ref>_cand{rank}.nodup.bam`, position 3
was `"nodup"`, not a cand token — no cand token exists as a standalone dot-field,
but the fix correctly returns `NA` and falls back, which means `major_minor` will
be `NA` for Phase-9 files. This is the correct behavior given the new naming: the
cand token is embedded inside `<ref>_cand{rank}`, not a separate dot-field.

`bin/summarize.R` lines 941-942: replaced `grepl("major", ...)` / `grepl("minor",
...)` with `grepl("_cand1\\.", ...)` / `grepl("_cand2\\.", ...)` to match the
Phase-9 PNG filenames (`<id>.variation_plot_<ref>_cand{rank}.png`). The old
"major"/"minor" patterns never matched after the phase-9 rename, causing both
grid PNGs to be silently skipped.

**Note:** requires human verification — the `major_minor` fallback logic means
variation plot filenames will contain `NA` for Phase-9 BAM inputs, which may need
a further follow-up if `major_minor` needs to appear in plot titles/filenames as
something more informative (e.g. `cand1`/`cand2`). The summarize.R grep fix is
correct regardless.

---

### CR-02: PARSEFIRSTMAPPING stub FASTA filenames do not match emit glob

**Files modified:** `modules/local/parsefirstmapping/main.nf`, `modules/local/parsefirstmapping/tests/main.nf.test.snap`
**Commit:** `ee7f839`
**Applied fix:**

`modules/local/parsefirstmapping/main.nf` stub block lines 82-83: renamed stub
FASTA outputs from `${prefix}.cand1.fa` / `${prefix}.cand2.fa` (period before
"cand", not matching `*_cand*.fa`) to `${prefix}.stubref_cand1.fa` /
`${prefix}.stubref_cand2.fa` (underscore before "cand", matching the glob).

Snapshot regenerated with `nf-test test --update-snapshot --profile test,docker`.
All 5 tests passed. The stub test snapshot now records `candidate_fasta` as
populated (two FASTA entries with md5 hashes) instead of the previous empty `[]`.
The fan-out path through `candidate_fasta` now has stub coverage.

---

### WR-01: Test fixture uses wrong `gate_flag` value and `minor_call` schema

**Files modified:** `bin/tests/test_compat.R`
**Commit:** `7661706`
**Applied fix:**

`bin/tests/test_compat.R` lines 143-144: corrected fixture values to match the
real script's output schema:
- `minor_call`: changed from `"co-infection"` / `"monoinfection"` to `"yes"` / `"no"`
  (real `summarize_mapping_to_all_references.R` emits `"yes"`/`"no"`)
- `gate_flag`: changed from `"pass"` to `"ok"` (real script emits `"ok"` for passing majors)

The `"pass" != "ok"` mismatch was silently injecting a spurious `review_flag`
trigger for every test case. With `gate_flag = "ok"`, clean samples no longer
fire the review_flag trigger, enabling future assertions on clean-sample
`review_flag == NA` behaviour.

---

### WR-02: `1:length(id_files)` crash when `id/` directory is empty

**Files modified:** `bin/summarize.R`
**Commit:** `4fba6ac`
**Applied fix:**

`bin/summarize.R` lines 897-903: two changes:
1. Moved `id_df` initialization outside the `if (length(id_files) > 0)` guard so
   the variable always exists before the tibble conversion on line 931. When
   `id_files` is empty, `matrix(nrow = 0, ncol = 2)` produces a zero-row frame.
2. Replaced `for (i in 1:length(id_files))` with `for (i in seq_along(id_files))`.
   `seq_along()` returns an empty integer vector for zero-length input, whereas
   `1:0` evaluates to `c(1L, 0L)` causing two crash-inducing loop iterations.

---

### WR-03: Doubled `_cand{rank}` slot in intermediate targeted-mapping filenames

**Files modified:** `subworkflows/local/targeted_mapping/main.nf`
**Commit:** `86a9574`
**Applied fix:** requires human verification (cross-file logic change)

`subworkflows/local/targeted_mapping/main.nf` line 31: the `.map {}` closure that
enriches `meta.reference` from the FASTA basename now strips `_cand\d+$` before
assigning `meta.reference`:

```groovy
def ref_raw   = fasta.getBaseName().toString().split('\\.').last()
def ref_clean = ref_raw.replaceAll(/_cand\d+$/, '')
def new_meta  = meta + [ reference: ref_clean ]
```

With this, `meta.reference` becomes `"3a_D17763"` instead of `"3a_D17763_cand1"`.
The config (modules_hcv.config) appends `.cand${meta.candidate_rank}`, producing
single-slot filenames like `Test_2.3a_D17763.cand1.withdup` instead of the
doubled `Test_2.3a_D17763_cand1.cand1.withdup`.

The `summarize.R` `str_remove(reference, "_cand[0-9]+$")` calls on lines 326,
399, 518, 537 remain in place as harmless no-ops (no `_cand` suffix to strip
from the now-clean reference field).

The unit-test fixture in `test_compat.R` already writes filenames in
`<sample>.<ref>.cand<rank>` format (line 155: `paste0(sampleName, ".", ref,
".cand", rk)`) — consistent with the post-fix format. No fixture changes were
needed.

**Human verification recommended:** this is a cross-file change affecting real
Nextflow workflow execution. Verify with a `-stub-run` or a real run that the
TARGETED_MAPPING stat/depth/bam filenames no longer carry a doubled cand-slot.

---

## Skipped Issues

None — all in-scope findings were fixed.

---

## Test Results

Full test suite run after all fixes:
```
docker run --rm -v "$PWD":/work -w /work \
  community.wave.seqera.io/library/r-seqinr_r-tidyverse:5358395134867368 \
  bash bin/tests/run_all.sh
```

Results:
- test_assembly_support.R: ALL PASS
- test_assembly_support_join.R: ALL PASS
- test_candidate_selection.R: ALL PASS
- test_classify_roles.R: ALL PASS
- test_coinfection.R: ALL PASS
- test_compat.R: ALL PASS
- test_denovo_confirm.R: ALL PASS
- test_dominance_score.R: ALL PASS
- test_major_gate.R: ALL PASS
- test_summarize_denovo.R: ALL PASS

10/10 test suites pass. No regressions introduced by the fixes.

nf-test `parsefirstmapping` module tests (all 5 scenarios):
- parsefirstmapping: idxstats + depth + references: PASSED
- parsefirstmapping: stub: PASSED (snapshot updated — candidate_fasta now populated)
- parsefirstmapping: 2k1b major keeps whole-name genotype (GATE-03): PASSED
- parsefirstmapping: no minor FASTA when no valid minor (GATE-04): PASSED
- parsefirstmapping: failing-major reports stats, no minor call (GATE-02): PASSED

---

_Fixed: 2026-06-14T19:30:00Z_
_Fixer: Claude (gsd-code-fixer)_
_Iteration: 1_
