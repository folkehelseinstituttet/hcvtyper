---
phase: 09-compatibility-filename-migration-regression-suite
reviewed: 2026-06-14T18:00:00Z
depth: standard
files_reviewed: 11
files_reviewed_list:
  - bin/summarize.R
  - bin/summarize_mapping_to_all_references.R
  - bin/tests/fixtures/compat_golden.csv
  - bin/tests/test_compat.R
  - conf/modules_hcv.config
  - modules/local/blastparse/main.nf
  - modules/local/blastparse/tests/main.nf.test.snap
  - modules/local/parsefirstmapping/main.nf
  - modules/local/parsefirstmapping/tests/main.nf.test
  - modules/local/parsefirstmapping/tests/main.nf.test.snap
  - workflows/hcvtyper.nf
findings:
  critical: 2
  warning: 3
  info: 2
  total: 7
status: issues_found
---

# Phase 9: Code Review Report

**Reviewed:** 2026-06-14T18:00:00Z
**Depth:** standard
**Files Reviewed:** 11
**Status:** issues_found

## Summary

Phase 9 renamed intermediate filename slots from `.major.`/`.minor.` to `.cand{rank}.` across
the pipeline. The candidate-rank join approach in `summarize.R` is mechanically sound: all three
stats loops (withdup, nodups, coverage) correctly strip `_cand[0-9]+$` from the reference token
before joining to `candidate_rank_lookup`. The PARSEFIRSTMAPPING and BLASTPARSE module outputs,
the Nextflow fan-out logic in `hcvtyper.nf`, and the regression test structure (COMPAT-01..04)
are well-considered.

Two blockers were found:

1. **The variation-plot grid is silently suppressed** after the rename. `plot_bam_variation.R`
   encodes the slot in output filenames using position 3 from the BAM basename — which now
   yields `nodup` (a suffix word) instead of `major`/`minor`. The `grepl("major", ...)` /
   `grepl("minor", ...)` split in `summarize.R` therefore always returns empty lists, so
   `Variation_plot_major.png` and `Variation_plot_minor.png` are never written.

2. **The PARSEFIRSTMAPPING stub produces FASTA files that never match the declared emit glob**,
   making `candidate_fasta` perpetually empty in `-stub-run` mode. A full-workflow stub run
   consequently never exercises the fan-out path.

Three warnings were found: a fixture `gate_flag` value mismatch that suppresses the intended
review-flag assertion coverage; the well-known `1:length(id_files)` crash path when no id files
exist; and an undocumented doubled `_cand{rank}` slot in intermediate targeted-mapping filenames.

---

## Critical Issues

### CR-01: Variation-plot grid silently suppressed after cand-slot rename

**File:** `bin/summarize.R:941-942` and `bin/plot_bam_variation.R:24,110,114,118`

**Issue:** `summarize.R` splits variation-plot PNGs into major/minor grids using
`grepl("major", ...)` / `grepl("minor", ...)` on filenames (lines 941-942). The PNG filenames
are produced by `plot_bam_variation.R`, which builds the slot token from position 3 of the
BAM basename (line 24: `str_split(basename(bam_file), "\\.")[[1]][3]`).

Under the Phase-9 naming, the BAM fed to `PLOT_BAMVARIATION` is named
`<id>.<ref>_cand{rank}.nodup.bam` (TARGETED_MAPPING:SAMTOOLS_SORMADUP prefix from
`conf/modules_hcv.config:131`). Position 3 of that dot-split is `nodup`, not `major`/`minor`.
The resulting PNG filename pattern becomes `<id>.variation_plot_<ref>_nodup.png`, which matches
neither `grepl("major", ...)` nor `grepl("minor", ...)`. Both `major_plots` and `minor_plots`
are empty for every sample, so neither grid PNG is ever written.

This is a silent regression: `summarize.R` writes no error; the grid PNGs simply never appear
in the output, nor in the MultiQC report.

**Fix:** Two coordinated changes are needed.

Option A (recommended) — change `plot_bam_variation.R` to derive the slot from the `cand{rank}`
field rather than position 3:
```r
# Replace line 24 in plot_bam_variation.R:
# OLD:
major_minor <- unlist(str_split(basename(bam_file), pattern = "\\."))[3]
# NEW (extract cand-slot; falls back to position 3 for legacy filenames):
parts <- unlist(str_split(basename(bam_file), pattern = "\\."))
major_minor <- parts[grep("^cand[0-9]+$", parts)][1]
if (is.na(major_minor)) major_minor <- parts[3]
```
Then update `summarize.R` lines 941-942 to match on `cand`:
```r
major_plots <- variation_plot_files[grepl("_cand1\\.", variation_plot_files)]
minor_plots <- variation_plot_files[grepl("_cand2\\.", variation_plot_files)]
```
Note: `plot_bam_variation.R` is not in the Phase-9 review scope but is the upstream root cause.
The summarize.R fix at lines 941-942 is in-scope and must be applied regardless.

---

### CR-02: PARSEFIRSTMAPPING stub FASTA filenames do not match the declared emit glob

**File:** `modules/local/parsefirstmapping/main.nf:82-83` vs `main.nf:27`

**Issue:** The `candidate_fasta` output is declared with glob `"*_cand*.fa"` (line 27),
which requires an underscore immediately before `cand`. The stub block creates:
```bash
: > ${prefix}.cand1.fa
: > ${prefix}.cand2.fa
```
For `prefix = "Test_2"` this produces `Test_2.cand1.fa` and `Test_2.cand2.fa`. The character
before `cand` is `.` (a period), not `_`, so neither file matches `*_cand*.fa`. The
`candidate_fasta` channel is therefore always empty in a `-stub-run`.

This is confirmed by the snapshot (`main.nf.test.snap:175-177`): the stub test records
`"candidate_fasta": []` while the real-run tests record actual FASTA paths. The stub-run of the
full workflow consequently never exercises the fan-out join (`hcvtyper.nf:400`) through
`candidate_fasta`, meaning the fan-out code path has no stub coverage.

**Fix:** Rename the stub outputs to match the glob pattern:
```bash
# Replace lines 82-83 in the stub block:
: > ${prefix}.ref_cand1.fa
: > ${prefix}.ref_cand2.fa
```
Or use a literal underscore-bearing fake reference name consistent with the real script's
output format (`<sample>.<ref>_cand{rank}.fa`):
```bash
: > ${prefix}.stubref_cand1.fa
: > ${prefix}.stubref_cand2.fa
```
After this change, regenerate the stub snapshot with `nf-test test --update-snapshot`.

---

## Warnings

### WR-01: Test fixture uses wrong `gate_flag` value, weakening review_flag assertion coverage

**File:** `bin/tests/test_compat.R:144`

**Issue:** The `run_summarize` harness writes `gate_flag = "pass"` into the synthetic
`parsefirstmapping.csv` fixture (line 144). The real `summarize_mapping_to_all_references.R`
emits `"ok"` for a passing major (line 196 of that script), never `"pass"`. In `summarize.R`
the review_flag trigger at line 1344 is `if (!is.na(gflag) && gflag != "ok")` — since
`"pass" != "ok"`, EVERY test case built by `run_summarize()` fires the review_flag trigger
"Major strain failed mapping quality thresholds" even for healthy samples.

The COMPAT-01/02/03/04 assertions check only the five core strain-call columns and do not
assert on `review_flag`, so no test currently fails. But the fixture silently injects a
spurious review flag, meaning the test does not verify that clean samples produce `NA`
review_flag. Any future assertion on `review_flag` in these cases would see a wrong value.

Additionally, line 143 writes `minor_call = if (nrow(minor) > 0) "co-infection" else "monoinfection"`.
The real script emits `"yes"` / `"no"` for this column. `minor_call` is not read by
`summarize.R` currently, but the schema mismatch is confusing and could cause failures if
`minor_call` consumption is added.

**Fix:**
```r
# Line 143-144 in test_compat.R:
minor_call = if (nrow(minor) > 0) "yes" else "no",
gate_flag  = "ok"
```

---

### WR-02: `1:length(id_files)` loop crashes when `id/` directory is empty

**File:** `bin/summarize.R:903`

**Issue:** The `id_df` data frame is created only when `length(id_files) > 0` (inside the
`if` block at lines 897-901). The `for` loop at line 903 is `for (i in 1:length(id_files))`.
When `length(id_files) == 0`, `1:0` in R evaluates to `c(1L, 0L)` — a two-element vector,
not an empty sequence. The loop body then references `id_df` (which does not exist) and
crashes with "object 'id_df' not found".

The workflow always stages at least one id file (INSTRUMENTID runs for every sample), so this
path is not reached in normal production runs. However it is reachable in unit tests and on
truncated/partial runs.

Note: the Phase 8 commit (71d6a24) corrected a related `if (length(id_files > 0))` bug
(misplaced parenthesis), but the `1:length()` anti-pattern in the loop body was not fixed.

**Fix:** Replace the anti-pattern with `seq_along` and guard `id_df` creation outside the `if`:
```r
# Lines 895-931 — replace the id_df block:
id_files <- list.files(path = path_9, pattern = "sequencerID.tsv$", full.names = TRUE)

id_df <- tibble(sampleName = character(), sequencer_id = character())

for (i in seq_along(id_files)) {
  # ... existing loop body, with id_df populated by bind_rows() or row assignment ...
}
```

---

### WR-03: Intermediate targeted-mapping filenames carry a doubled `_cand{rank}` slot

**File:** `subworkflows/local/targeted_mapping/main.nf:31` and `conf/modules_hcv.config:131,240,248,228`

**Issue:** The PARSEFIRSTMAPPING script produces FASTAs named `<sample>.<ref>_cand{rank}.fa`
(e.g. `Test_2.3a_D17763_cand1.fa`). Inside `TARGETED_MAPPING`, `meta.reference` is populated
by:
```groovy
fasta.getBaseName().toString().split('\\.').last()
// → "3a_D17763_cand1"  (the last dot-field of the FASTA basename)
```
The config then builds stat/depth prefixes as:
```
"${meta.id}.${meta.reference}.cand${meta.candidate_rank.toInteger()}.withdup"
// → "Test_2.3a_D17763_cand1.cand1.withdup"
```
This produces a doubled slot (`_cand1.cand1`) in intermediate filenames. `summarize.R` extracts
the reference at field position 2 (`3a_D17763_cand1`) and strips `_cand[0-9]+$` to recover
`3a_D17763`, so the join works correctly. The doubled slot is therefore not a correctness bug
in the current implementation — `summarize.R` handles it — but:

1. The doubling is undocumented and will confuse future maintainers.
2. Any tool that reads the intermediate filenames outside summarize.R would see `cand1.cand1`.
3. Published files in `${params.outdir}/samtools/` carry the doubled slot, making the output
   directory harder to interpret.

**Fix:** Extract only the `<ref>_cand{rank}` part (without the sample prefix) at the reference
enrichment step:
```groovy
// In targeted_mapping/main.nf line 31 — split on '.' and take the second-to-last field,
// which is the reference without the sample prefix:
def parts  = fasta.getBaseName().toString().split('\\.')
def new_meta = meta + [ reference: parts.size() > 1 ? parts[-1] : parts[0] ]
```
This keeps `meta.reference = "3a_D17763_cand1"` unchanged (`.last()` already does this), but
the config prefixes should then NOT append `cand${rank}` again. Alternatively, strip the
`_cand{rank}` suffix from `meta.reference` in the enrichment:
```groovy
def ref_raw  = fasta.getBaseName().toString().split('\\.').last()
def ref_clean = ref_raw.replaceAll(/_cand\d+$/, '')
def new_meta  = meta + [ reference: ref_clean ]
```
Then remove the `.cand${meta.candidate_rank.toInteger()}` suffix from all six config prefixes
(lines 131, 228, 240, 248, 257, 263 in `modules_hcv.config`). This would make `summarize.R`'s
strip no-op (correct). The `summarize.R` `str_remove(reference, "_cand[0-9]+$")` strip would
then be a no-op but harmless.

---

## Info

### IN-01: `BLASTPARSE.out.candidate_fasta` declared but never consumed in workflow

**File:** `modules/local/blastparse/main.nf:22` and `workflows/hcvtyper.nf` (no reference)

**Issue:** BLASTPARSE declares `emit: candidate_fasta` (line 22 of main.nf) and the stub
creates `${prefix}.cand1.fa` and `${prefix}.cand2.fa`, but `BLASTPARSE.out.candidate_fasta`
is never referenced in `hcvtyper.nf`. The canonical source of candidate FASTAs is
`PARSEFIRSTMAPPING.out.candidate_fasta`. The dead emit produces confusing stub output and
suggests an incomplete or abandoned migration step.

**Fix:** Either remove the `candidate_fasta` emit from BLASTPARSE entirely (if not planned for
use), or add a comment explaining that it is reserved for a future assembly-guided FASTA
selection path.

---

### IN-02: `versions` channel drop for GET_MAPPING_STATS subworkflows (pre-existing)

**File:** `workflows/hcvtyper.nf:336,360`

**Issue:** Lines 336 and 360 write:
```groovy
versions = GET_MAPPING_STATS_WITHDUP.out.versions
versions = GET_MAPPING_STATS_MARKDUP.out.versions
```
These assign to a new undeclared local variable `versions` instead of mixing into `ch_versions`.
The versions from both GET_MAPPING_STATS subworkflows are silently dropped from the
`software_mqc_versions.yml` report. This is a pre-existing issue not introduced in Phase 9, but
it affects the correctness of the published software versions report.

**Fix:**
```groovy
// Line 336:
ch_versions = ch_versions.mix(GET_MAPPING_STATS_WITHDUP.out.versions)
// Line 360:
ch_versions = ch_versions.mix(GET_MAPPING_STATS_MARKDUP.out.versions)
```

---

_Reviewed: 2026-06-14T18:00:00Z_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
