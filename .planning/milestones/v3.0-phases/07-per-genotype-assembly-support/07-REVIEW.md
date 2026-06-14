---
phase: 07-per-genotype-assembly-support
reviewed: 2026-06-13T00:00:00Z
depth: standard
files_reviewed: 7
files_reviewed_list:
  - bin/assembly_support_join.R
  - bin/blast_parse.R
  - bin/summarize.R
  - bin/tests/test_assembly_support.R
  - bin/tests/test_assembly_support_join.R
  - modules/local/blastparse/main.nf
  - workflows/hcvtyper.nf
findings:
  critical: 2
  warning: 4
  info: 2
  total: 8
status: issues_found
---

# Phase 7: Code Review Report

**Reviewed:** 2026-06-13
**Depth:** standard
**Files Reviewed:** 7
**Status:** issues_found

## Summary

Phase 7 adds a per-subtype "assembly support" roll-up in `blast_parse.R` (§4b, emitting
`*.assembly_support.csv`) and a genotype-level left-join of that support onto Phase-6
candidates via the sourceable `join_assembly_support()` helper, wired into `summarize.R`
and staged through `workflows/hcvtyper.nf`.

The single-sample / in-memory logic is correct and the helper's design intent (no row loss,
NA-fill, parameterized match level) is sound. However, two reproducible BLOCKER-class defects
sit on the **default multi-sample pipeline path** and were both masked by the test fixtures:

1. The join key type mismatches between the candidate side (`candidate_genotype` read as
   numeric by readr) and the support side (`genotype_from_subtype()` returns character),
   so `left_join` aborts on essentially every real run at the **default** `match_level="genotype"`.
2. `map_dfr()` over a mix of header-only and populated `*.assembly_support.csv` /
   `*.candidates.csv` files aborts SUMMARIZE — exactly the skip-assembly / no-mapping mixed
   batch the empty-guards claim to support.

Both were missed because the unit tests build in-memory tibbles (character genotype, no
header-only files) and never exercise the `read_csv` → `map_dfr` → `join` chain that the
real `summarize.R` uses.

## Critical Issues

### CR-01: Join key type mismatch aborts `join_assembly_support()` at the default `match_level="genotype"`

**File:** `bin/assembly_support_join.R:78` (candidate side) and `:96` (support side); triggered from `bin/summarize.R:545`

**Issue:**
The candidate-side match key is `candidate_genotype` (line 78). In `summarize.R` this column
comes from `read_csv(*.candidates.csv)` (line 506), where readr **infers numeric** for genotype
values like `3`, `1`, `6` (a purely-digit column). The support-side match key is
`genotype_from_subtype(subtype)` (line 96), which always returns **character** (`substr(...,1,1)`
/ `if_else`). `dplyr::left_join(by = c("sampleName", ".match_key"))` aborts on incompatible
key types:

```
Error: Can't join `x$.match_key` with `y$.match_key` due to incompatible types.
ℹ `x$.match_key` is a <double>.
ℹ `y$.match_key` is a <character>.
```

Reproduced directly: a `candidates.csv` with genotype column `3,1` reads as `<double>`, and the
join against a character `.match_key` errors. Since the documented/parameter default is
`match_level="genotype"` and HCV genotypes 1–7 are numeric, **this crashes SUMMARIZE on
virtually every real run.** The only non-crashing genotype is the `2k1b` recombinant (forces the
column to character). The `match_level="subtype"` path happens to be safe only because
`candidate_subtype` (e.g. `3a`) reads as character.

The unit test `test_assembly_support_join.R` does not catch this: `mk_cand()` (line 42) builds
`candidate_genotype` via `genotype_from_subtype()` in-memory (character), never via `read_csv`,
so the production readr-coercion path is never exercised.

**Fix:** Coerce both sides of the key to character before the join so the key type is
deterministic regardless of how the candidates CSV was typed:

```r
cand <- candidates_df %>%
  mutate(.match_key = as.character(
    if (match_level == "subtype") candidate_subtype else candidate_genotype
  ))
# ... and in the support branch:
mutate(.match_key = as.character(
  if (match_level == "subtype") subtype else genotype_from_subtype(subtype)
))
```

Alternatively pin the column type at read time in `summarize.R:506`
(`read_csv(.x, col_types = cols(candidate_genotype = col_character(), candidate_rank = col_integer()))`).
Prefer fixing it inside the helper so it is robust to any caller. Add a join-helper test that
feeds a numeric `candidate_genotype` (round-tripped through `read_csv`/`write_csv`) to lock it in.

### CR-02: `map_dfr()` over mixed header-only + populated CSVs aborts SUMMARIZE (defeats the skip-assembly empty-guard)

**File:** `bin/summarize.R:526` (`support_df`) and `bin/summarize.R:506` (`candidates_long`)

**Issue:**
`blast_parse.R` **always** writes `*.assembly_support.csv` — header-only when a sample has no
contigs (line 187–196). `summarize_mapping_to_all_references.R` likewise **always** writes
`*.candidates.csv` — header-only (0-row) on the no-mapping branch (line 201, `candidates_long`
left at 0 rows). readr types every column of a header-only CSV as `<character>`, while a
populated CSV from another sample types the metric columns as `<double>`. `map_dfr()` then calls
`bind_rows()` across these and aborts:

```
Error: Can't combine `..1$best_contig_length` <character> and `..2$best_contig_length` <double>.
```

Reproduced directly with a header-only + a 1-row support CSV. This is precisely the
multi-sample batch (one skip-assembly / no-mapping sample alongside a normal one) that the
PLUMB-02 / T-07-03 "typed-empty-tibble, never abort" comments at lines 499–539 claim to handle.
The typed-empty fallback only fires when the file LIST is empty (`length(files) == 0`); it does
nothing when files exist but one is header-only — which is the common case. Both new reads
(`support_df` and `candidates_long`) are affected; the pre-existing `df_blast_out` read at
line 483 has the same shape but is out of this phase's scope.

`df_denovo` (line 448) is NOT affected because `*.blastparse.csv` always carries exactly one
data row (line 330–345 of `blast_parse.R`), so its types are consistent.

**Fix:** Pin column types at read time so header-only and populated files combine cleanly:

```r
support_df <- map_dfr(support_files, ~ read_csv(.x, col_types = cols(
  sample                 = col_character(),
  subtype                = col_character(),
  best_contig_length     = col_double(),
  best_contig_pident     = col_double(),
  best_contig_aln_length = col_double(),
  best_contig_kmer_cov   = col_double()
))) %>% rename(sampleName = sample)
```

Apply the analogous `col_types = cols(...)` to the `candidates_long` read at line 506
(at minimum `candidate_genotype = col_character()`, `candidate_rank = col_integer()`, metrics
`col_double()`), which also resolves CR-01's root cause at the source. Add a summarize-level or
join-level test that mixes a header-only CSV with a populated one.

## Warnings

### WR-01: No-candidate run silently drops the `cand_*` assembly-support columns from `Summary.csv`

**File:** `bin/summarize.R:572-574`

**Issue:**
When `candidate_support` has zero rows (no candidates anywhere in the batch), the else branch
sets `candidate_support_wide <- tibble(sampleName = character())`. This frame contributes **no**
`cand_<rank>_assembly_support*` columns to the `left_join` at line 822, so `Summary.csv` silently
loses the entire assembly-support column block on that run. This is the same schema-drift failure
mode that the denovo_* PLUMB-02 guard (lines 457–469) was explicitly written to prevent — but the
fix was not mirrored here. Downstream consumers (MultiQC config, manual review) that expect a
stable column set will break or mis-render.

**Fix:** Declare the wide columns explicitly in the empty branch so the schema is stable, e.g.
build a zero-row tibble with `sampleName` plus `cand_1_assembly_support`, `cand_2_assembly_support`,
and the per-rank metric columns for `seq_len(params.n_candidates)`. Alternatively gate the
downstream `left_join` and add the columns via `add_column()` when absent (mirror the GLUE-absent
guard at lines 936–987).

### WR-02: `pivot_wider` column count is data-dependent, so `Summary.csv` schema varies run to run

**File:** `bin/summarize.R:559-571`

**Issue:**
`pivot_wider(names_from = candidate_rank, ...)` produces `cand_<rank>_*` columns only for the
ranks that actually appear in the batch. A batch where every sample is single-candidate yields
only `cand_1_*` columns; a batch with a passing rank-2 yields `cand_1_*` and `cand_2_*`. The
resulting `Summary.csv` column set therefore changes depending on the input batch, which is
fragile for any automated downstream parser and for cross-run diffing. Combined with WR-01 this
makes the assembly-support block non-deterministic in shape.

**Fix:** Complete the wide frame to a fixed `1..params.n_candidates` rank set (e.g.
`tidyr::complete()` on `candidate_rank` before pivot, or reindex columns against the known rank
range) so the same columns are always emitted regardless of which ranks are populated.

### WR-03: `blast_parse.R` empty-support `best_contig_aln_length` typed `integer`, populated path typed from `length` — inconsistent with the join helper's `double`

**File:** `bin/blast_parse.R:193` (vs populated `:181`) and `bin/assembly_support_join.R:91`

**Issue:**
The empty `support_tbl` declares `best_contig_aln_length = integer(0)` (line 193) while the
populated path carries `length` (from `read_tsv`, inferred numeric/double). The join helper's
typed-empty support frame (line 91) declares this column `double()`. These three differing type
declarations for the same logical column are a latent type-conflict surface (the same class that
produced CR-02). It does not currently crash because the affected reads route through the
`nrow == 0` typed-empty branch, but it is an accident waiting to be tripped by a future refactor.

**Fix:** Use one consistent numeric type (`double`) for `best_contig_aln_length` across
`blast_parse.R` (both empty and populated paths) and `assembly_support_join.R`. Centralize the
six-column support schema in one place if practical.

### WR-04: `match_level` is unvalidated — any value other than `"subtype"` silently means "genotype"

**File:** `bin/assembly_support_join.R:78,96`

**Issue:**
The branch `if (match_level == "subtype") ... else ...` treats every non-`"subtype"` value
(including a typo like `"geneotype"`, `NA`, or an empty string from a mis-parsed arg) as the
genotype path with no warning. `summarize.R:32` parses `denovo_match_level` from `args[7]` with
only a "non-empty" guard, so a malformed value flows straight in and is silently misinterpreted
rather than rejected. For a clinical genotyping tool, a silent misread of the match-level
parameter is a correctness risk.

**Fix:** Validate up front, e.g.
`stopifnot(match_level %in% c("genotype", "subtype"))` at the top of `join_assembly_support()`,
and mirror the check where `denovo_match_level` is parsed in `summarize.R`.

## Info

### IN-01: Test fixtures bypass the real `read_csv` typing path, hiding CR-01/CR-02

**File:** `bin/tests/test_assembly_support_join.R:42-67`

**Issue:**
`mk_cand()` / `mk_support()` build in-memory tibbles with hand-typed columns (character genotype,
numeric metrics) and the tests never write/read a CSV. As a result the readr type-inference that
breaks the production path (numeric `candidate_genotype`, character header-only columns) is never
exercised, so a green test suite gives false confidence. The `test_assembly_support.R` subprocess
test does write/read CSVs but only for a single prefix, so it never hits the multi-file
`map_dfr` combine.

**Fix:** Add at least one test that round-trips candidates and support through
`write_csv`/`read_csv` (numeric genotype) and one that `map_dfr`-combines a header-only file with
a populated one, then runs the join — locking in CR-01 and CR-02 fixes.

### IN-02: Stale CLI-compat comment in `blast_parse.R` header

**File:** `bin/blast_parse.R:7-9`

**Issue:**
The usage banner still states `<references> and <agens> are kept for CLI compatibility but no
longer used by this script`, but `references` IS now used (read at lines 31–37 and consumed by
`write_ref_fasta`). Only `agens` is unused. The comment is misleading for maintainers.

**Fix:** Update the banner to note that `references` is used and only `agens` is retained for
CLI compatibility.

---

_Reviewed: 2026-06-13_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
