# Phase 9: Compatibility, Filename Migration + Regression Suite - Pattern Map

**Mapped:** 2026-06-14
**Files analyzed:** 8 (6 modified, 1 new, 2 regenerated snapshots)
**Analogs found:** 8 / 8

> This is a refactor/migration phase, not a greenfield build. Every file already
> exists or has a direct in-repo analog. "Pattern to copy from" here means: the
> existing two-slot code being replaced (modify-in-place) plus the sibling file
> that models the convention for the one genuinely new file (`test_compat.R`).

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|-------------------|------|-----------|----------------|---------------|
| `conf/modules_hcv.config` (6 `ext.prefix` closures) | config | transform (filename slot) | itself (the 6 closures already carry the N-arm ternary) | exact (self) |
| `workflows/hcvtyper.nf` (guard L85-87 + fan-out L406-444) | workflow/route | event-driven (channel fan-out) | itself + `subworkflows/local/targeted_mapping/main.nf` | exact (self) |
| `modules/local/parsefirstmapping/main.nf` (emit + stub) | module | file-I/O (emit N FASTAs) | itself (current 2-slot emit) | exact (self) |
| `modules/local/blastparse/main.nf` (stub) | module | file-I/O (stub FASTAs) | `modules/local/parsefirstmapping/main.nf` | role-match |
| `bin/summarize.R` (3 parse blocks + cv_by_ref) | utility (R script) | transform (CSV join) | itself (candidates_long load L570-580) | exact (self) |
| `bin/tests/test_compat.R` (NEW) | test | request-response (subprocess) + transform (function-source) | `bin/tests/test_candidate_selection.R` + `bin/tests/test_classify_roles.R` | exact |
| `modules/local/parsefirstmapping/tests/main.nf.test.snap` | test (build artifact) | snapshot | `nf-test --update-snapshot` (tool-generated) | exact |
| `modules/local/blastparse/tests/main.nf.test.snap` | test (build artifact) | snapshot | `nf-test --update-snapshot` (tool-generated) | exact |

## Pattern Assignments

### `conf/modules_hcv.config` (config, filename-slot transform) — COMPAT-02 / D-04

**Analog:** the six closures themselves; they already contain the full N-arm ternary.

**Six closures to change** (verified line numbers):
- L131 `TARGETED_MAPPING:SAMTOOLS_SORMADUP` — `meta.` namespace
- L228 `TARGETED_MAPPING:SAMTOOLS_DEPTH` — **`meta1.` namespace** (joined input — preserve this)
- L240 `TARGETED_MAPPING:STATS_WITHDUP` — `meta.`
- L248 `TARGETED_MAPPING:STATS_MARKDUP` — `meta.`
- L257 `TARGETED_MAPPING:IVAR_CONSENSUS` — `meta.`, no `reference` field (`${meta.id}.<slot>`)
- L263 `TARGETED_MAPPING:CONSENSUS_DISTANCE` — `meta.`, no `reference` field

**Current pattern** (the two-slot special-case to delete — L131):
```groovy
ext.prefix = { "${meta.id}.${meta.reference}.${meta.candidate_rank.toInteger() == 1 ? 'major' : (meta.candidate_rank.toInteger() == 2 ? 'minor' : "cand${meta.candidate_rank}")}.nodup" }
```

**Migration target** (uniform cand-slot — drop the major/minor arms):
```groovy
ext.prefix = { "${meta.id}.${meta.reference}.cand${meta.candidate_rank.toInteger()}.nodup" }
```

**Namespace caution** (L228 SAMTOOLS_DEPTH uses `meta1.`, NOT `meta.`):
```groovy
ext.prefix = { "${meta1.id}.${meta1.reference}.cand${meta1.candidate_rank.toInteger()}.nodup" }
```

**Caution — V5 input validation (Pitfall: NA.toInteger()):** `.toInteger()` is only
safe inside these closures because `candidate_rank` is guaranteed present for a mapped
candidate. Do NOT coerce rank to Integer in the workflow fan-out (see hcvtyper.nf below).

---

### `workflows/hcvtyper.nf` (workflow/route, channel fan-out) — COMPAT-02 / D-01 / D-04

**Analog:** the existing fan-out (L406-444) + the subworkflow it feeds (`subworkflows/local/targeted_mapping/main.nf`).

**Sub-change A — delete the N>2 guard** (L85-87):
```groovy
// DELETE entirely — D-01 lifts this; the message points at this very phase.
if (params.n_candidates > 2) {
    error "params.n_candidates = ${params.n_candidates} is not yet supported: ... N>2 support arrives with the Phase 9 filename-slot migration (COMPAT-02)."
}
```

**Sub-change B — N-FASTA fan-out** (L406-444). Current two-slot pick (L430-431):
```groovy
// Pick the per-rank FASTA from the legacy shim emits (rank 1 -> _major.fa, else _minor.fa).
def fasta = (rank == '1') ? major_fasta : minor_fasta
```
This is the hard two-slot ceiling. Generalize to a rank-indexed lookup against an
N-FASTA emit from PARSEFIRSTMAPPING (see module below). Per RESEARCH Open Q2, lowest-churn
option is a single `candidate_fasta` tuple emit `[meta(with candidate_rank), fasta]` joined
to the candidates fan-out by full meta, letting the existing `flatMap` keep iterating rows.

**Patterns that MUST be preserved verbatim** (anti-regression, RESEARCH Pitfall 3):
```groovy
// 1. remainder: true semantics on the FASTA join (optional emit must not drop the sample)
.join(PARSEFIRSTMAPPING.out.<fasta_emit>, remainder: true)
// 2. rank/confirmation_status stay STRINGS — never .toInteger() here (NA crash guard)
def rank = new_meta.candidate_rank.toString()
// 3. the null-FASTA + confirmation_status guard
.filter { entry -> entry[0]['confirmation_status'] == 'pass' && entry[1] != null }
// 4. the id == sample assert
assert new_meta.id == new_meta.sample : "Metadata mismatch: id=${new_meta.id}, sample=${new_meta.sample}"
```

---

### `modules/local/parsefirstmapping/main.nf` (module, file-I/O emit) — D-01 / D-04 / D-03

**Analog:** itself. Current output emits (L19-29):
```groovy
tuple val(meta), path("*.parsefirstmapping.csv"), path("*major.fa"), emit: major_mapping, optional: true
tuple val(meta), path("*.parsefirstmapping.csv"), path("*minor.fa"), emit: minor_mapping, optional: true
tuple val(meta), path("*.parsefirstmapping.csv"),                    emit: csv,           optional: true
tuple val(meta), path("*.candidates.csv"),                           emit: candidates,    optional: true
```
Replace the two `major_mapping`/`minor_mapping` emits with an N-FASTA emit (e.g.
`emit: candidate_fasta` globbing `*.cand*.fa`, or a per-rank tuple). The R selection
script (`bin/summarize_mapping_to_all_references.R`) writes the per-rank FASTAs — confirm
it emits `*.cand{rank}.fa` (or update it) so the module glob matches.

**Stub block** (L57-90) — D-03 forbids hardcoded `.major.fa`/`.minor.fa`. Current (L79-80):
```bash
: > ${prefix}.major.fa
: > ${prefix}.minor.fa
```
Migrate to N cand-slot stub files (one per candidate row already emitted at L74-76).
The stub candidates.csv already has 2 rows (ranks 1 and 2) — emit `${prefix}.cand1.fa`
and `${prefix}.cand2.fa` to match. **Stub convention** (lowercase, deterministic, `: >`
or `printf`, no heredoc for data files) is established repo-wide.

---

### `modules/local/blastparse/main.nf` (module, file-I/O stub) — D-03

**Analog:** `modules/local/parsefirstmapping/main.nf` stub block (same convention).

**Output emits** (L22-23) and **stub** (L68-69) both hardcode `*major.fa`/`*minor.fa`:
```groovy
tuple val(meta), path("*major.fa")      , emit: major_fasta, optional: true
tuple val(meta), path("*minor.fa")      , emit: minor_fasta, optional: true
```
```bash
: > ${prefix}.major.fa
: > ${prefix}.minor.fa
```
Rename to cand-slot names consistent with the parsefirstmapping migration. **Caution:**
verify which downstream consumer reads BLASTPARSE's `major_fasta`/`minor_fasta` emits
before renaming the channel identifiers (D-04 says channel renames are discretionary/can
be a follow-up; the FILE names in the stub are the hard COMPAT-02 requirement).

---

### `bin/summarize.R` (utility R script, CSV join transform) — COMPAT-02 / D-02

**Analog:** the already-present `candidates_long` load (L570-580) — this is the join source.

**Reuse this loaded frame for the D-02 join** (L570-580, already in the script):
```r
candidates_long <- map_dfr(candidates_files, ~ read_csv(.x, col_types = cols(
  sample              = col_character(),
  candidate_rank      = col_integer(),
  candidate_ref       = col_character(),
  ...
))) %>%
  rename(sampleName = sample)
```

**Four slot-coupled sites to refactor** (all replace `first_major_minor` logic):

1. **Block 1 — withdup stats** (L246-306). Position-3 parse at **L261**:
   ```r
   tmp_df$first_major_minor[i] <- str_split(basename(stats_files[i]), "\\.")[[1]][3]
   ```
   drives `case_when(first_major_minor == "major"/"minor" ~ ...)` at L282-283, L288-289,
   L298-299, plus the `str_remove(..., "_major")/"_minor")` cleanup at L290-291.

2. **Block 2 — nodup stats** (L310-363). Position-3 parse at **L326**; same case_when
   pattern at L343-344 (`Major_genotype_mapping`/`Minor_genotype_mapping`), L347-349,
   L354-357.

3. **Block 3 — coverage** (L383-492). Position-3 parse at **L393**; case_when at L473-486.

4. **Block 4 — cv_by_ref** (L462-467) — **the easy-to-miss fourth site** (RESEARCH Pitfall 2).
   It strips the slot from a *reference token*, not position-3, so a grep for `[[1]][3]` misses it:
   ```r
   mutate(candidate_ref = str_remove(reference, "_(major|minor)$")) %>%
   ```

**Migration approach (D-02):** in each stats loop, parse `reference` (position 2) and
`sampleName` (position 1) as today, but instead of reading position-3 join each per-file
row to `candidates_long` on `(sampleName, candidate_ref)` and read `candidate_rank`, then
map `candidate_rank == 1 -> Major_*`, `candidate_rank == 2 -> Minor_*` (generalizable to N).
The cv_by_ref strip becomes `_cand[0-9]+$` (or is dropped if the slot leaves the reference field).

**Highest-risk constraint (RESEARCH Pitfall 1 + the join key at L366 / L1043):**
`Major_reference`/`Minor_reference` are a JOIN KEY downstream and MUST stay byte-identical
(cleaned ref name, no slot suffix):
```r
# L366
df_mapped_reads <- full_join(df_with_dups, df_nodups, join_by(sampleName, Major_reference, Minor_reference))
# L1043
left_join(df_coverage, join_by(sampleName, Major_reference, Minor_reference))
```
A mismatch here silently empties Major_* columns / Cartesians the join — no crash.

**COMPAT-03 (validate-only, D-05) — DO NOT modify these emit sites:**
- Legacy `Major_*`/`Minor_*` mapping-stat columns: produced across L282-488.
- Role columns `Major_role_*`/`Minor_role_*` + `overall_sample_call`: L695-737 (build),
  L1055-1057 + L1325-1329 (final Summary.csv select). Phase 9 only asserts both sets exist.

---

### `bin/tests/test_compat.R` (test — NEW) — TEST-01 / D-08

This is the only genuinely new file. It has TWO analogs (one per pattern); pick per-case
per D-08 discretion.

**Analog A — subprocess pattern** (`bin/tests/test_candidate_selection.R`). Use for
COMPAT-01/02/03 (system under test is the whole `summarize.R` script — monolithic, not
sourceable). Copy the self-location + fixture-in-tempdir + `system2` harness.

**Self-location boilerplate** (test_candidate_selection.R L46-51 — identical in both analogs):
```r
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
this_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
bin_dir  <- normalizePath(file.path(this_dir, ".."))
```

**fail/ok helper convention** (both analogs — no testthat; base-R asserts + exit code):
```r
fail <- function(msg) { cat("FAIL:", msg, "\n"); quit(status = 1) }
ok   <- function(msg) cat("PASS:", msg, "\n")
```

**Subprocess invocation** (test_candidate_selection.R L96-108):
```r
old <- getwd(); setwd(wd); on.exit(setwd(old), add = TRUE)
exit <- system2("Rscript",
  c(shQuote(script), shQuote(idx_path), shQuote(depth_path),
    shQuote(sampleName), shQuote(refs_path), minRead, minCov, n_candidates),
  stdout = FALSE, stderr = FALSE)
out <- if (file.exists(out_path)) read_csv(out_path, show_col_types = FALSE) else NULL
```
**Note for COMPAT-02:** the synthetic fixture must build `*.cand1.*.stats`/`*.cand2.*.stats`
files AND a matching `*.candidates.csv` (with `candidate_rank`) so the new join has its key,
then assert `summarize.R` yields populated `Major_*`/`Minor_*` rows (Pitfall 1 regression catch).

**Analog B — function-sourcing pattern** (`bin/tests/test_classify_roles.R`). Use for
COMPAT-04 if exercising the helper end-to-end. `source()` + in-memory tibble builder:
```r
source(file.path(bin_dir, "genotype_utils.R"))
source(file.path(bin_dir, "classify_roles.R"))
# builder + classify() wrapper — test_classify_roles.R L45-67
classify <- function(df, minRead = 500, minCov = 30) {
  classify_roles(score_candidates(df), minRead = minRead, minCov = minCov,
                 denovo_min_contig_length = 1000, denovo_min_kmer_cov = 2.0,
                 denovo_min_blast_identity = 90, match_level = "genotype")
}
```
**COMPAT-04 cases already exist at helper level** in test_classify_roles.R (Test 3 sim1
1a:1b allowed L122-132; Test 6 2k1b blocked L173-184). D-08 wants `test_compat.R` to exercise
the FULL path through `summarize.R` (integration), so the subprocess pattern is likely lower
friction here too — reuse Analog A and assert on `overall_sample_call`/role columns in the
emitted Summary.csv.

**COMPAT-01 fixture anchor** (D-06/D-07): assert ONLY core strain-call columns
(`Major_reference`, `Minor_reference`, `Major_genotype_mapping`, `Minor_genotype_mapping`,
`overall_sample_call`). `bin/tests/fixtures/flagoff_golden.csv` (2 rows: S1 co-infection
1a/2b, S2 monoinfection 3a) is the available anchor — note its columns are
`sampleName,Major_reference,Minor_reference,major_typable,minor_typable` (no genotype_mapping
or overall_sample_call), so the two D-07 golden cases must add those expected values.
**RESEARCH A2/Open-Q1: flag a `checkpoint:human-verify`** to confirm fixture values equal
true v1.0/v2.0 output before locking assertions.

**Auto-discovery (no CI change):** `run_all.sh` globs `"$here"/test_*.R` (L19) — `test_compat.R`
is picked up automatically; the `r-regression` CI job (`.github/workflows/ci.yml` L45-60) runs
`bash bin/tests/run_all.sh` in the pinned container. No YAML edit.

---

### `modules/local/{parsefirstmapping,blastparse}/tests/main.nf.test.snap` (regenerate) — D-03

**Analog:** tool-generated — never hand-edit md5s.

**Stale entries confirmed:**
- parsefirstmapping snap: `*_major.fa`/`*_minor.fa` and `Test_2.major.fa`/`Test_2.minor.fa`
  at L19, L61, L103, L112, L151, L160, L199, L208.
- blastparse snap: `toy.major.fa`/`toy.minor.fa` at L26, L34, L85, L93.

**Regenerate command** (nf-test 0.9.2 in conda env `NEXTFLOW`):
```bash
nf-test test modules/local/parsefirstmapping/tests/main.nf.test --profile test,docker --update-snapshot
nf-test test modules/local/blastparse/tests/main.nf.test --profile test,docker --update-snapshot
```
**Caution (RESEARCH Pitfall 4):** before committing, eyeball the snap diff — the ONLY
expected change is `.major.fa`/`.minor.fa` -> `.cand1.fa`/`.cand2.fa`. Any other md5 churn
(esp. `candidates.csv` content) is a real regression. RESEARCH A3: run the full nf-test
suite once and regen any OTHER snap that diffs only on the slot.

## Shared Patterns

### R test self-location + fail/ok harness
**Source:** `bin/tests/test_candidate_selection.R` L46-57, `bin/tests/test_classify_roles.R` L28-40
**Apply to:** `test_compat.R` (mandatory boilerplate so it runs from any cwd under `run_all.sh`)
```r
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
this_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
bin_dir  <- normalizePath(file.path(this_dir, ".."))
fail <- function(msg) { cat("FAIL:", msg, "\n"); quit(status = 1) }
ok   <- function(msg) cat("PASS:", msg, "\n")
cat("\nALL PASS\n")   # final line on success
```

### Slot-from-rank ternary (the migration hinge)
**Source:** `conf/modules_hcv.config` L131/228/240/248/257/263
**Apply to:** all 6 TARGETED_MAPPING ext.prefix closures
```groovy
// FROM: ...${rank==1?'major':(rank==2?'minor':"cand${rank}")}...
// TO:   ...cand${meta.candidate_rank.toInteger()}...   (meta1. for SAMTOOLS_DEPTH)
```

### String-typed candidate_rank in Groovy (NA crash guard)
**Source:** `workflows/hcvtyper.nf` L415-444 (Pitfall 3 comment)
**Apply to:** the N-FASTA fan-out generalization — keep `candidate_rank`/`confirmation_status`
as Strings in the workflow; only `.toInteger()` inside config closures where rank is guaranteed.

### Byte-identical Major_reference/Minor_reference join key
**Source:** `bin/summarize.R` L366, L1043
**Apply to:** all four D-02 refactor sites — the cleaned reference name must not gain a slot suffix.

### Stub-file convention (deterministic, no heredoc for data)
**Source:** `modules/local/parsefirstmapping/main.nf` L60-89
**Apply to:** parsefirstmapping + blastparse stub blocks — `: >` / `printf`, lowercase cand-slot
filenames, stable `versions.yml`.

### Container parity for the R suite (host has no tidyverse)
**Source:** `.github/workflows/ci.yml` L53-60
**Apply to:** every local run of `test_compat.R` / `run_all.sh`
```bash
docker run --rm -v "$PWD":/work -w /work \
  community.wave.seqera.io/library/r-seqinr_r-tidyverse:5358395134867368 \
  bash bin/tests/run_all.sh
```

## No Analog Found

| File | Role | Data Flow | Reason |
|------|------|-----------|--------|
| (none) | — | — | All Phase 9 files exist in-repo or have a direct sibling analog. The single new file (`test_compat.R`) is fully modeled by two existing test files. |

## Metadata

**Analog search scope:** `conf/`, `workflows/`, `subworkflows/local/targeted_mapping/`,
`modules/local/{parsefirstmapping,blastparse}/`, `bin/`, `bin/tests/`, `bin/tests/fixtures/`,
`.github/workflows/`
**Files scanned:** 13 (CONTEXT, RESEARCH, CLAUDE.md probe, 2 test analogs, modules_hcv.config,
summarize.R [4 ranges], run_all.sh, flagoff_golden.csv, hcvtyper.nf [3 ranges], 2 module main.nf,
2 snap files via grep, ci.yml via grep)
**Pattern extraction date:** 2026-06-14
**Note:** `./CLAUDE.md` does not exist; project conventions live in
`.planning/codebase/CONVENTIONS.md` (channel naming lowercase+underscores, stub block structure,
`task.ext.prefix` pattern, R naming) — planner should read it before implementing.
