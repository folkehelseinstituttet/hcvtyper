---
phase: 08-dominance-scoring-strain-role-classification
reviewed: 2026-06-14T12:00:00Z
depth: standard
files_reviewed: 7
files_reviewed_list:
  - bin/classify_roles.R
  - bin/summarize.R
  - bin/tests/test_classify_roles.R
  - bin/tests/test_dominance_score.R
  - bin/tests/test_summarize_denovo.R
  - modules/local/summarize/main.nf
  - workflows/hcvtyper.nf
findings:
  critical: 2
  warning: 4
  info: 2
  total: 8
status: issues_found
---

# Phase 08: Code Review Report

**Reviewed:** 2026-06-14T12:00:00Z
**Depth:** standard
**Files Reviewed:** 7
**Status:** issues_found

## Summary

Phase 8 introduces `bin/classify_roles.R` as a pure sourced helper (score_candidates, classify_roles, is_valid_minor), wires it into `bin/summarize.R`, stages it through the Nextflow module, and adds test coverage for dominance scoring and role classification.

The core classifier logic in `classify_roles.R` is sound: the dominance scoring formula, the asymmetric refute rule (D-11), the HCV exceptions (D-12), the overall-call derivation (D-14), and the zero-row/NULL guard (T-08-01) all hold up under trace. The nine test cases in `test_classify_roles.R` correctly exercise the decision boundaries including the calibration-anchor 2.0 k-mer floor.

Two blockers are present. First, the enriched `candidates.csv` produced by `summarize.R` (D-16) is not declared in the Nextflow module's `output:` block and will therefore never be published or available to downstream processes. Second, the pre-existing `length(id_files > 0)` guard in `summarize.R` is incorrect — it will crash with "object 'id_df' not found" if no sequencerID.tsv files are present. Four warnings cover a stub schema gap (five Phase-8 columns missing), a dead `header` vector that leaves the MultiQC TSV without identification comments, a column-name inconsistency for a drug column in the GLUE-absent path, and a breadth-fraction edge case in score_candidates that inflates scores for candidates with ≤1% coverage (though these candidates always fail the gate so final calls are unaffected).

---

## Critical Issues

### CR-01: `candidates.csv` not declared as a Nextflow `output:` emit

**File:** `modules/local/summarize/main.nf:36-40`
**Issue:** `summarize.R` writes `candidates.csv` to the work directory (line 678) — the Phase-8 enriched long output carrying role / dominance_score / role_reason / overall_sample_call for every candidate including background (D-16). The file is never declared in the `output:` block of the SUMMARIZE process. Without an output declaration, Nextflow does not capture the file: it cannot be published to `params.outdir`, referenced by a downstream `path(candidates)` input, or picked up by any `publishDir` rule. The D-16 requirement to "surface every candidate, never drop backgrounds" is silently not met in an actual pipeline run.

**Fix:** Add an output declaration and a publishDir stanza:

```groovy
// In modules/local/summarize/main.nf output: block
path 'candidates.csv'  , emit: candidates

// In conf/modules_hcv.config SUMMARIZE publishDir array
[
    path: { "${params.outdir}/summary" },
    mode: params.publish_dir_mode,
    pattern: 'candidates.csv'
],
```

Also add the file to the stub block so `-stub` runs produce a header-only `candidates.csv` that matches the real schema.

---

### CR-02: `length(id_files > 0)` — incorrect guard causes crash when no sequencer-ID files are staged

**File:** `bin/summarize.R:857`
**Issue:** The parentheses are wrong. `length(id_files > 0)` evaluates `id_files > 0` (a character vector compared to a number, which returns NA values with a warning) and then takes the length of that result vector. When `id_files` is non-empty the expression accidentally returns a truthy value and `id_df` is created; the loop on line 863 runs correctly. However, when `id_files` is empty (`character(0)`), the expression returns `0` (falsy), so `id_df` is never created. The loop on line 863 then executes `1:length(id_files)` = `1:0` = `c(1, 0)` — two iterations — and the first iteration (`i=1`) attempts `id_df$sampleName[1] <- ...`, immediately crashing with `Error: object 'id_df' not found`. Any run where INSTRUMENTID emits no files (e.g. all reads filtered, empty FASTQ) will abort SUMMARIZE.

```r
# Line 857 — current (broken):
if (length(id_files > 0)) {

# Fix:
if (length(id_files) > 0) {
```

---

## Warnings

### WR-01: Stub schema missing five Phase-8 output columns

**File:** `modules/local/summarize/main.nf:73,78`
**Issue:** The stub's hardcoded `Summary.csv` header (lines 73 and 78) is missing the following columns that `summarize.R` writes unconditionally in Phase 8:
- `review_flag` (the rewired review sentinel, added at line 1274 of summarize.R)
- `denovo_major_subtype` (line 1141)
- `denovo_minor_subtype` (line 1142)
- `denovo_major_subtype_match` (line 1143)
- `denovo_minor_subtype_match` (line 1148)

Any downstream process or integration test that runs with `-stub` and validates the `Summary.csv` schema will see a different column set than a real run produces. This silently breaks stub-based schema checks and any MultiQC or reporting step that expects these columns from the stub output.

**Fix:** Extend the stub header on both lines 73 (CSV) and 78 (TSV) to include the five missing columns after `Minor_role_reason`. Also add corresponding placeholder values in the sample data rows on lines 74 and 79.

---

### WR-02: Dead `header` vector — MultiQC TSV written without identification comment lines

**File:** `bin/summarize.R:1363-1376`
**Issue:** The `header` character vector (lines 1363–1367) defines the MultiQC comment block (`# id: 'summary'`, `# section_name: 'Summary'`, etc.) but is never written to the file. Lines 1376 and 1379 write only the column-name row and the data rows. The `summary_mqc.tsv` is therefore missing the `#`-prefixed header that MultiQC uses to identify the module and assign a section name. Depending on the MultiQC configuration (`fn:` pattern matching), the file may be parsed with a generic or missing section label. The variable is fully dead code.

**Fix:** Write the header lines before the column row:
```r
# Replace lines 1375-1379 with:
write_lines(header, file)                                          # comment block
tt %>% colnames() %>% paste0(collapse = "\t") %>%
  write_lines(file, append = TRUE)                                 # column names
write_tsv(tt, file, append = TRUE)                                 # data rows
```
Also update `# format: 'csv'` in the header to `# format: 'tsv'` since the file format changed in commit d8690d2.

---

### WR-03: Column name inconsistency — `daclasvir` vs `daclatasvir` in GLUE-absent path

**File:** `bin/summarize.R:1175-1177`
**Issue:** When the GLUE report is absent, `summarize.R` adds the drug-resistance column via `add_column("daclasvir" = NA_character_, ...)` (line 1175). The stub schema in `modules/local/summarize/main.nf:73` uses `daclatasvir` (the correct drug name and presumably what the real GLUE report supplies). This means the column name in the Summary.csv differs depending on whether GLUE ran:
- GLUE present: column name is `daclatasvir` (from GLUE output)
- GLUE absent: column name is `daclasvir` (misspelling in add_column block)

Any downstream consumer, MultiQC config, or report template that references this column by name will fail silently in one of the two paths.

**Fix:**
```r
# Line 1175 — fix the column name:
"daclatasvir"       = NA_character_,
"daclatasvir_mut"   = NA_character_,
"daclatasvir_mut_short" = NA_character_,
```

---

### WR-04: Breadth-fraction boundary error in `score_candidates()` for coverage ≤ 1%

**File:** `bin/classify_roles.R:135-137`
**Issue:** `score_candidates()` coerces the breadth source from percent to fraction using `ifelse(breadth_src > 1, breadth_src / 100, breadth_src)`. The `candidate_cov` column (the fallback breadth source when `cand_cov_breadth` is absent) holds values on a 0–100 percent scale (confirmed in `bin/summarize_mapping_to_all_references.R:55-56`). For any candidate with `candidate_cov` between 0 and 1 exclusive (i.e., fractional-percent coverage like 0.5%), the condition `breadth_src > 1` is FALSE, so the value is left as-is and treated as a 0–1 fraction. A candidate with 0.5% coverage is therefore scored as if it has 50% breadth coverage — a 100-fold overestimate.

In the current production flow this does not affect final call accuracy because such candidates always fail the `cov > minCov` gate (minCov = 30 in production) and are classified as `background/below_floor` before the dominance score matters. However, the inflated `dominance_score` value is written to `candidates.csv` (CR-01) and could mislead any analyst reviewing that file.

**Fix:** Use a threshold of `>= 2` or an explicit known-columns flag, or document that `candidate_cov` will always be ≥ 1 for gated candidates so this boundary is safe:
```r
# Option A: raise threshold slightly (benign for percent-scale data, still ambiguous):
breadth_frac <- ifelse(is.na(breadth_src), 0,
                       ifelse(breadth_src >= 2, breadth_src / 100, breadth_src))

# Option B (preferred): accept an explicit `breadth_is_pct` parameter:
score_candidates <- function(df, ..., breadth_is_pct = TRUE) {
  breadth_frac <- if (breadth_is_pct) breadth_src / 100 else breadth_src
  ...
}
```

---

## Info

### IN-01: Multiple `1:length(x)` loops not guarded against empty vectors

**File:** `bin/summarize.R:174,252,316,383`
**Issue:** Four loops use `for (i in 1:length(x))`. When `x` is empty, `1:length(x)` evaluates to `c(1, 0)` — two iterations over an empty vector — instead of zero iterations. For these four specific loops (kraken_files, stats_files×2, cov_files), the preceding pre-allocation `matrix(nrow = length(x), ...)` produces a zero-row frame, so the body attempts to write to a valid (but empty) data frame and may silently skip rather than crash. Still, this is an unsafe pattern that should be `seq_along(x)` for clarity and correctness. (Not introduced by Phase 8; pre-existing.)

**Fix:**
```r
for (i in seq_along(kraken_files)) { ... }
for (i in seq_along(stats_files))  { ... }  # both withdup and markdup blocks
for (i in seq_along(cov_files))    { ... }
```

---

### IN-02: `test_summarize_denovo.R` exercises a retired code path without exercising the replacement

**File:** `bin/tests/test_summarize_denovo.R:119-213`
**Issue:** Block 3 (Tests 3–7) and Block 4 directly test `apply_denovo_layer()`, which the Phase-8 D-15 retirement removed from the production call path in `summarize.R`. The tests correctly verify that `summarize.R` does NOT call `apply_denovo_layer()` (line 81), and the legacy function is still sourced so tests pass, but the test file has no integration tests for the new production path: `score_candidates()` + `classify_roles()` invoked from `summarize.R`'s wiring (the cv_by_ref join, the role-to-wide mapping, or the review_flag rewiring). This means a regression in the wiring code would not be caught by the existing test suite. The dedicated `test_classify_roles.R` covers the pure functions, but the `summarize.R` integration wiring (lines 642–739) has no coverage.

**Fix:** Add an integration test fixture in `test_summarize_denovo.R` (or a new `test_summarize_classify_integration.R`) that exercises the full `candidate_support -> cv_by_ref join -> score_candidates -> classify_roles -> candidate_roles_wide` pipeline with a minimal in-memory input, asserting on `candidate_roles_wide` columns and on `overall_sample_call` values in the wide frame.

---

_Reviewed: 2026-06-14T12:00:00Z_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
