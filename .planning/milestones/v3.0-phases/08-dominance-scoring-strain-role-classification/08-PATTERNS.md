# Phase 8: Dominance Scoring + Strain-Role Classification - Pattern Map

**Mapped:** 2026-06-13
**Files analyzed:** 8 (2 new helpers, 2 new tests, 4 modified)
**Analogs found:** 8 / 8 (every new/modified file has an in-repo analog)

This is a brownfield R-logic-only phase. Every new file has a close, recently-written
analog already in the repo (the Phase-3 `denovo_confirm.R` and Phase-7
`assembly_support_join.R` helpers + their tests are the canonical templates). The
planner should copy structure verbatim from these analogs, not invent shape.

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|-------------------|------|-----------|----------------|---------------|
| `bin/classify_roles.R` (NEW; may also host `score_candidates()`) | utility (pure sourced R helper) | transform | `bin/assembly_support_join.R` | exact (same house pattern, same domain) |
| `bin/dominance_score.R` (NEW, optional — research recommends folding into `classify_roles.R`) | utility (pure sourced R helper) | transform | `bin/denovo_confirm.R` | exact |
| `bin/tests/test_classify_roles.R` (NEW) | test | batch/assert | `bin/tests/test_assembly_support_join.R` | exact |
| `bin/tests/test_dominance_score.R` (NEW) | test | batch/assert | `bin/tests/test_assembly_support_join.R` | exact |
| `bin/tests/fixtures/` rows (NEW) | test fixture | n/a | existing `bin/tests/fixtures/` + inline `mk_*` builders | exact |
| `bin/summarize.R` (MODIFY) | service (terminal aggregator entrypoint) | transform / file-I/O | self (existing cov loop, confirmation layer, review_flag, joins) | self |
| `modules/local/summarize/main.nf` (MODIFY) | config (Nextflow module) | request-response | self + `path()` staging precedent in same file | self |
| `conf/modules_hcv.config` (MODIFY, ~L328) | config (ext.args wiring) | n/a | the `SUMMARIZE` `ext.args` line itself (denovo_* precedent) | self |
| `nextflow.config` (MODIFY, ~L45-53) | config (param defaults) | n/a | the `denovo_*` / `n_candidates` param block | self |
| `nextflow_schema.json` (MODIFY, ~L191-225) | config (nf-schema) | n/a | the `denovo_min_contig_length` / `n_candidates` schema entries | self |
| `is_valid_minor()` (RECOVER from git, port into `classify_roles.R`) | utility (logic fragment) | transform | recovered verbatim from `43904de~1:bin/summarize_mapping_to_all_references.R` | exact (recovered below) |

## Pattern Assignments

### `bin/classify_roles.R` (NEW pure sourced R helper)

**Analog:** `bin/assembly_support_join.R` (closest — same phase-family, same join-over-candidates
domain, same typed-empty discipline). Secondary: `bin/denovo_confirm.R` for the
floors/asymmetric-return idiom.

**House-pattern header + defensive tidyverse guard** (`bin/assembly_support_join.R:37-42`):
```r
# Defensive: consumers already load tidyverse, which provides the pipe,
# group_by()/slice_max()/mutate()/left_join()/if_else(). This guard only fires
# if sourced into a session that has not.
if (!exists("group_by")) {
  library(tidyverse)
}
```
Copy this verbatim. A pure sourced helper defines functions ONLY — no `commandArgs`,
no top-level file I/O, no global mutation. The doc-comment block at the top of
`assembly_support_join.R:3-35` is the template for documenting purpose + purity contract.

**Match-level guard** (`bin/assembly_support_join.R:50`) — copy for any new
`match_level` arg the classifier exposes:
```r
stopifnot(match_level %in% c("genotype", "subtype"))
```

**Floors-as-named-args with calibration-validated defaults** (`bin/denovo_confirm.R:36-38`):
```r
classify_minor_denovo <- function(blast_out_df, major_geno, minor_geno,
                                  min_len = 1000, min_kmer = 2.0, min_pid = 90,
                                  match_level = "genotype") {
```
Mirror this signature shape for `classify_roles(candidates_df, minRead, minCov,
denovo_min_contig_length = 1000, denovo_min_kmer_cov = 2.0, denovo_min_blast_identity = 90,
match_level = "genotype", score_weights = list(...))`. NOTE D-10 cites validated
1000/2.0/90, but the runtime args override the defaults (see Shared Patterns → Calibration
mismatch).

**Substantiality (ANDed floors) — D-10 corroboration verdict** (`bin/denovo_confirm.R:62-66`):
```r
top <- top %>%
  mutate(
    match_key = if (match_level == "subtype") subtype else genotype_from_subtype(subtype),
    substantial = sc_length >= min_len & kmer_cov >= min_kmer & pident >= min_pid
  )
```
For Phase 8 the candidate's joined columns are
`cand_<rank>_assembly_support_best_contig_length` / `_kmer_cov` / `_pident` (already
present after the Phase-7 join). The corroboration test is the same ANDed-floor form
against those columns.

**Asymmetric refute ordering — D-11** (`bin/denovo_confirm.R:68-80`). This is the
template for "did de novo work for the dominant?" Copy the `any(... & match_key == X,
na.rm = TRUE)` presence idiom and the asymmetric return order:
```r
has_minor_sub <- any(top$substantial & top$match_key == minor_geno, na.rm = TRUE)
has_major_sub <- any(top$substantial & top$match_key == major_geno, na.rm = TRUE)
if (has_minor_sub) return("confirmed_by_denovo")
if (has_major_sub) return("refuted")          # major present, candidate absent
"unconfirmed"                                  # de novo failed overall -> never refute
```
Map to Phase 8: candidate-supported → `co-infection`/`corroborated`; candidate-unsupported
BUT dominant-supported → `background`/`refuted_denovo`; neither supported → keep candidate as
`co-infection`/`uncorroborated_kept` (D-11 — never suppress when de novo was inconclusive).
⚠️ Anti-pattern (research): the new classifier assigns roles from scratch, so the old
downgrade-only invariant no longer holds automatically — the exception/refute rules may only
demote to `background`, never promote a below-floor candidate.

**Typed zero-row guard — DoS / skip-assembly** (`bin/assembly_support_join.R:63-80`).
Copy this exact "if null/zero-row, return a typed frame carrying every output column, never
abort" block for the no-candidate / no-mapping path. Generalize it to the new
`role` / `dominance_score` / `role_reason` columns:
```r
if (is.null(candidates_df) || nrow(candidates_df) == 0) {
  out <- candidates_df
  if (is.null(out)) out <- tibble()
  out <- out %>% mutate(<new typed columns: character()/double()>)
  return(out)
}
```

---

### `is_valid_minor()` port — D-12 (recovered from git history; embed in `classify_roles.R`)

**Source:** `git show 43904de~1:bin/summarize_mapping_to_all_references.R` (the function was
removed in commit `43904de` "feat(06-01)"; research A5 confirmed it is no longer in any live
`bin/` file). **Recovered verbatim** (reconstruct as a pure function taking explicit
`cand_subtype`/`cand_genotype`/`dom_subtype`/`dom_genotype` rather than the old `major_*`
closure vars):
```r
is_valid_minor <- function(minor_row) {
    minor_subtype <- minor_row$Subtype
    minor_genotype <- minor_row$Genotype

    # Rule: allow 1a and 1b co-infection
    if ((major_subtype %in% c("1a", "1b")) & (minor_subtype %in% c("1a", "1b")) & (major_subtype != minor_subtype)) {
      return(TRUE)
    }

    # Rule: block 2k1b co-infections with any genotype 1 or 2 (and itself)
    if ((major_genotype == "2k1b" & minor_genotype %in% c("1", "2", "2k1b")) |
        (minor_genotype == "2k1b" & major_genotype %in% c("1", "2", "2k1b"))) {
      return(FALSE)
    }

    # Rule: allow only different genotypes
    return(major_genotype != minor_genotype)
  }
```
The three rules to preserve EXACTLY (COMPAT-04): (1) allow 1a/1b cross-subtype; (2) block
2k1b paired with genotype 1/2/2k1b; (3) otherwise require different genotype. Use
`genotype_from_subtype()` (already sourced) for the genotype comparisons — do NOT hand-roll
`substr`. Note the historical 2k1b block is broader than D-12's wording ("block 2k1b pairs"):
it blocks 2k1b vs {1,2,2k1b}. Port the historical predicate verbatim per D-12 ("port
`is_valid_minor()` verbatim").

**`genotype_from_subtype()` — the canonical 2k1b-aware helper** (`bin/genotype_utils.R:24-26`),
already sourced in `summarize.R:10`:
```r
genotype_from_subtype <- function(subtype) {
  if_else(subtype == "2k1b", subtype, substr(subtype, 1, 1))
}
```

---

### `bin/summarize.R` (MODIFY — entrypoint)

**Analog:** self. Five concrete edit sites, all read in-session:

**1. Source the new helper** (`bin/summarize.R:10-16`) — add `source("classify_roles.R")`
after the existing sources (genotype_utils.R first so `genotype_from_subtype()` is in scope):
```r
source("genotype_utils.R")
source("denovo_confirm.R")   # D-15: retire CONSUMPTION (file may stay or be deleted per planner)
source("denovo_layer.R")     # D-15: retire CONSUMPTION
source("assembly_support_join.R")
# NEW: source("classify_roles.R")
```

**2. Arg parse — append `score_weight_*` at the END** (`bin/summarize.R:29-45`). The existing
defensive index-guarded parse is the template; new args MUST go after `args[11]`
(`n_candidates`) — see Shared Patterns → Positional ext.args coupling:
```r
min_targeted_read <- if (length(args) >= 9  && nchar(args[9])  > 0) as.numeric(args[9])  else NA_real_
min_targeted_cov  <- if (length(args) >= 10 && nchar(args[10]) > 0) as.numeric(args[10]) else NA_real_
n_candidates      <- if (length(args) >= 11 && nchar(args[11]) > 0) as.integer(args[11]) else 2L
# NEW score_weight_* parse continues at args[12], args[13], ... with the same guard form.
```

**3. CV-evenness inside the cov loop — D-03** (`bin/summarize.R:357-413`). The per-position
depth vector `cov$X3` exists ONLY inside this loop and is discarded after L413 — the CV factor
MUST be computed here. Add a column to `tmp_df` (widen `ncol`/`colnames` at L357-358) and
compute it next to `avg_depth` (L379):
```r
cov <- read_tsv(cov_files[i], col_names = FALSE)
ref_length <- nrow(cov)
tmp_df$avg_depth[i] <- mean(cov$X3)
# NEW (D-03): CV over per-position depth incl. zeros (SAMTOOLS_DEPTH runs -aa,
# conf/modules_hcv.config:259), mapped to a 0–1 evenness factor. Guard zero-mean.
cv_raw <- if (ref_length > 0 && mean(cov$X3) > 0) sd(cov$X3) / mean(cov$X3) else NA_real_
tmp_df$cv_evenness[i] <- if (!is.na(cv_raw)) 1 / (1 + cv_raw) else 0
```
Then carry `cv_evenness` through the `df_coverage` pivot at L418-441 (the
`Major_/Minor_*` `case_when` block) so the score function can read it per candidate.
⚠️ Anti-pattern: computing CV after L413 is impossible — `cov$X3` is gone.

**4. Retire the legacy confirmation layer — D-15** (`bin/summarize.R:965-981`). The
`apply_denovo_layer()` call + `coinfection_flag` mutate are REPLACED by the new role classifier:
```r
final <- apply_denovo_layer(            # <- REMOVE (D-15): replaced by classify_roles()
  ...,
  denovo_min_contig_length, denovo_min_kmer_cov, ...
) %>%
  mutate(coinfection_flag = if_else(    # <- REMOVE (D-15)
    minor_typable == "NO" & !is.na(minor_denovo_status) & minor_denovo_status == "confirmed_by_denovo",
    "possible_multiple_strains", ...
  ))
```
The new classifier runs over the long `candidate_support` frame
(`bin/summarize.R:579`, the Phase-7 join output — this is the natural classifier input, one
row per candidate) and emits `role` / `dominance_score` / `role_reason` per candidate; the
dominant + corroborated co-infection then fill the wide `Major_*`/`Minor_*` slots (D-16).

**5. Rewire `review_flag` onto roles — D-13/D-15** (`bin/summarize.R:1096-1124`). Keep the
`pmap_chr` full-sentence idiom; swap the trigger columns from
`minor_denovo_status`/`coinfection_flag`/`gate_flag` onto `role`/`role_reason`/sample-call:
```r
final <- final %>%
  mutate(review_flag = {
    pmap_chr(
      list(<role / role_reason / overall_sample_call columns>),
      function(...) {
        msgs <- character(0)
        if (<condition>) msgs <- c(msgs, "<full human sentence derived from role_reason>")
        ...
        if (length(msgs) == 0) NA_character_ else paste(msgs, collapse = " | ")
      }
    )
  })
```
Also update the column-reorder `select()` at L1127-1164 (drop `minor_denovo_status` /
`coinfection_flag`, add `overall_sample_call` + the new role columns), and the
`write_csv(final, "Summary.csv")` (L1167) + `summary_mqc.tsv` (L1181-1187) schema follow
automatically. The enriched long `*.candidates.csv` (D-16, CLASS-03 — every candidate incl.
background) is a NEW output written from `candidate_support` after classification (the file is
currently READ-only at L511-529; add a `write_csv(candidate_support_with_roles, ...)`).

**Pinned `col_types` precedent** for any new candidate columns read back
(`bin/summarize.R:519-528`) — pin `candidate_genotype = col_character()` etc. to avoid
readr inferring digit-genotypes as double (CR-01/CR-02).

---

### `bin/tests/test_classify_roles.R` and `test_dominance_score.R` (NEW)

**Analog:** `bin/tests/test_assembly_support_join.R` (exact template).

**Self-locating source block** (`bin/tests/test_assembly_support_join.R:24-39`) — copy verbatim,
adding `source(file.path(bin_dir, "classify_roles.R"))`:
```r
suppressPackageStartupMessages(library(tidyverse))
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
this_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
bin_dir  <- normalizePath(file.path(this_dir, ".."))
source(file.path(bin_dir, "genotype_utils.R"))
source(file.path(bin_dir, "classify_roles.R"))
fail <- function(msg) { cat("FAIL:", msg, "\n"); quit(status = 1) }
ok   <- function(msg) cat("PASS:", msg, "\n")
```

**In-memory fixture builders** (`bin/tests/test_assembly_support_join.R:42-67`) — copy the
`mk_cand()` / `mk_support()` tibble-builder idiom; assert on the REAL sourced function
(no inline re-implementation — research anti-pattern). The runner
(`bin/tests/run_all.sh:19`) auto-globs `test_*.R`, so no registration needed.

**Headline calibration assertions** (handoff §2 evidence table): false 4g
(53279 reads, breadth 0.677, spiky/high-CV, no 4g contig, dominant 1a HAS a contig)
→ `role == "background"`, `role_reason == "refuted_denovo"`; ERR1810447 (full 9207 bp 2b)
and ERR1810453 (2949 bp partial 2b, k-mer cov ~5) → `co-infection`/`corroborated`;
sim1 1a:1b + sim2 2a:3a true co-infections preserved; IVT extreme-ratio genuine minors →
`uncorroborated_kept` when de novo failed.

---

### `modules/local/summarize/main.nf` (MODIFY)

**Analog:** self. Two edits:

**Stage the new helper as a `path()` input** — mirror the existing helper staging
(`modules/local/summarize/main.nf:30-33`):
```r
    path(genotype_utils)
    path(denovo_confirm)
    path(denovo_layer)
    path(assembly_support_join)
    // NEW: path(classify_roles)
```
(Also add the corresponding `file()` arg at the `SUMMARIZE(...)` call site in
`workflows/hcvtyper.nf` — mirror how `assembly_support_join` is passed.)

**Update the stub `Summary.csv` header — D-15 + Pitfall 4** (`modules/local/summarize/main.nf:72-78`).
The stub hard-codes the full wide header ending in `...denovo_minor_contig_length,minor_denovo_status`.
Remove `minor_denovo_status` (D-15) and append `overall_sample_call` + the new role columns, in
BOTH the `Summary.csv` heredoc (L72) AND the `summary_mqc.tsv` printf (L77). The
`tests/default.nf.test.snap` then needs regen.
⚠️ MEMORY "Phase 7 SUMMARIZE snapshot regen deferred" — fold this into the same pending
regen (needs docker; MEMORY "Host disk near-full"; nf-test via `PATH=~/.nf-test:$PATH ...
--profile docker`).

---

### `conf/modules_hcv.config` (MODIFY ~L328) — D-04 param wiring

**Analog:** self — the `SUMMARIZE` `ext.args` line, which already wires the `denovo_*` precedent
(`conf/modules_hcv.config:327-328`):
```groovy
withName: 'FOLKEHELSEINSTITUTTET_HCVTYPER:HCVTYPER:SUMMARIZE' {
    ext.args   = { "${params.denovo_min_contig_length} ${params.denovo_min_kmer_cov} ${params.denovo_min_blast_identity} ${params.denovo_match_level} ${params.denovo_confirm_minor} ${params.minRead} ${params.minCov} ${params.n_candidates}" }
```
Append `${params.score_weight_*}` (and the evenness-transform constant) AFTER
`${params.n_candidates}` — never insert mid-string (Pitfall 2: positional coupling; the R side
reads `args[n]` by index).

---

### `nextflow.config` (MODIFY ~L45-53) + `nextflow_schema.json` (MODIFY ~L191-225)

**Analog:** the `denovo_*` param block (`nextflow.config:45-53`) and its schema entries
(`nextflow_schema.json:191,206,217`). Declare each `score_weight_*` default with an inline
comment in the same block, and add matching `nextflow_schema.json` properties (type/default/
description) mirroring `denovo_min_contig_length` / `n_candidates`.
⚠️ Also resolve the calibration-default mismatch here (Shared Patterns below).

## Shared Patterns

### Pure sourced R helper (the house pattern)
**Source:** `bin/assembly_support_join.R:37-50`, `bin/denovo_confirm.R:29-38`, `bin/genotype_utils.R:20-26`
**Apply to:** `bin/classify_roles.R` (+ any `dominance_score.R`)
- Functions only; no `commandArgs`, no top-level file I/O, no global mutation.
- Defensive `if (!exists("group_by")) library(tidyverse)` guard.
- `stopifnot(match_level %in% c("genotype","subtype"))` for any match-level arg.
- Calibration-validated numeric defaults as named function args (1000/2.0/90).
- Doc-comment block stating purpose + purity contract + design decisions.

### Typed-empty / DoS guard (PLUMB-02 / T-03-01 / CR-01,02)
**Source:** `bin/assembly_support_join.R:63-80,97-107`; `bin/denovo_confirm.R:47-49`; `bin/summarize.R:519-543,565-572`
**Apply to:** every new column the classifier emits, and every new CSV read.
- Null/zero-row input → return a typed zero-row frame, NEVER `stop()`.
- Pin `col_types` on read (esp. `candidate_genotype = col_character()`); `as.character` coerce
  join keys (digit-genotype-as-double otherwise aborts the join).

### Positional `ext.args` coupling (Pitfall 2)
**Source:** `conf/modules_hcv.config:328` ↔ `bin/summarize.R:29-45`
**Apply to:** `score_weight_*` wiring.
- The `ext.args` string is space-joined positional; `summarize.R` reads `args[n]` by index.
- Append new args ONLY at the END (after `n_candidates`/`args[11]`), in both the Groovy string
  and the R parse, with the defensive `if (length(args) >= N && nchar(args[N]) > 0) ... else <default>` form.

### Calibration-default mismatch (Pitfall 1 — MUST resolve before scoring)
**Source:** `nextflow.config:46-47` (ships `denovo_min_contig_length = 500`,
`denovo_min_kmer_cov = 10.0`) vs `bin/denovo_confirm.R:37` defaults (1000 / 2.0, validated).
**Apply to:** `nextflow.config` + `nextflow_schema.json`.
- The R-helper 1000/2.0 defaults are overridden by the args, which carry the config's 500/10.0
  → the EFFECTIVE corroboration floor is 500/10.0, not the validated 1000/2.0.
- 10.0 k-mer is STRICTER than validated 2.0 and could wrongly refute ERR1810453's partial 2b
  (k-mer cov ~5) — a CLASS-02 correctness risk.
- Planner must explicitly (a) correct `nextflow.config`+schema to 1000/2.0, or (b) re-validate
  500/10.0 against the evidence table. Resolve BEFORE score calibration (the corroboration
  verdict is a score-classifier input).

### CV-evenness must live inside the cov loop (D-03)
**Source:** `bin/summarize.R:357-413` (the cov loop); `conf/modules_hcv.config:259` (`-aa` verifies
zeros present)
**Apply to:** the new `cv_evenness` column.
- `cov$X3` exists only inside L360-413; compute `1/(1+sd/mean)` there, guard the zero-mean edge.

### pmap_chr human-sentence builder (D-13)
**Source:** `bin/summarize.R:1096-1124`
**Apply to:** the rewired `review_flag` and the derived `role_reason` sentence.
- Keep the `pmap_chr` over per-row trigger columns → `character(0)` accumulator → `paste(...,
  collapse = " | ")` → `NA_character_` when empty. MultiQC orange-highlight is already wired for
  `review_flag` (`assets/multiqc_config.yml`), so keep the column name.

### R-emits decisions, Nextflow-routes
**Source:** project architecture (CLAUDE.md) + this being the terminal SUMMARIZE step.
**Apply to:** whole phase — NO new channels/processes/containers. All classification logic stays
in R; the only Nextflow edits are param wiring + staging one new helper file.

## No Analog Found

None. Every new/modified file has a close in-repo analog (the Phase-3/Phase-7 helper +
test pair, and the existing `summarize.R` / config wiring). The single recovered artifact
(`is_valid_minor()`) is reproduced verbatim above from git history (commit `43904de~1`).

## Metadata

**Analog search scope:** `bin/` (R helpers + tests), `modules/local/summarize/`, `conf/`,
`nextflow.config`, `nextflow_schema.json`, git history of
`bin/summarize_mapping_to_all_references.R`.
**Files scanned:** `assembly_support_join.R`, `denovo_confirm.R`, `genotype_utils.R`,
`summarize.R` (L1-60, 355-444, 508-646, 790-829, 1075-1188), `test_assembly_support_join.R`,
`run_all.sh`, `modules/local/summarize/main.nf`, `conf/modules_hcv.config` (L320-344),
`nextflow.config` (L44-57), plus git-recovered `is_valid_minor()`.
**Pattern extraction date:** 2026-06-13
