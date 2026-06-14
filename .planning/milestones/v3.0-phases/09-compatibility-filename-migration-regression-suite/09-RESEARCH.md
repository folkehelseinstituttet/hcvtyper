# Phase 9: Compatibility, Filename Migration + Regression Suite - Research

**Researched:** 2026-06-14
**Domain:** Nextflow (nf-core) pipeline filename-slot migration + R (tidyverse) summary refactor + R/nf-test regression suite
**Confidence:** HIGH (codebase-verified; this is an internal-integration phase, no external-library research needed)

## Summary

Phase 9 is a pure-integration, non-breaking phase over an existing, well-instrumented codebase. Almost every fact needed to plan it is verifiable directly in the repository — there is no novel external library or framework to research. The four COMPAT requirements and TEST-01 decompose into four independent tracks already enumerated in `09-CONTEXT.md`: (1) the `.major.`/`.minor.` → `.cand1.`/`.cand2.`/`.cand{N}.` filename-slot migration plus the `summarize.R` parsing cutover, (2) validate-only legacy-column presence, (3) golden-baseline non-breaking proof, and (4) a new `bin/tests/test_compat.R` wired automatically into the existing `r-regression` CI job.

The single most important codebase discovery: **the `.major.`/`.minor.` filename slot is produced by `ext.prefix` closures in `conf/modules_hcv.config`, NOT inside the module `main.nf` files.** [VERIFIED: conf/modules_hcv.config L131,L228,L240,L248,L257,L263] Those six closures already contain a ternary that emits `cand${meta.candidate_rank}` for rank ≥ 3 — so the rank-1/rank-2 special-casing to `'major'`/`'minor'` is the only thing that needs to change to make every slot `.cand{rank}.`. The matching consumer is three `str_split(basename(...), "\\.")[[1]][3]` extractions in `bin/summarize.R` plus a fourth dependency: a `str_remove(reference, "_(major|minor)$")` in the `cv_by_ref` block and several `str_remove(..., "_major")/"_minor")` reference-cleanup steps. D-02 replaces the position-3 extractions with a join on `candidate_rank` from the already-loaded `candidates_long` frame.

The second discovery: the workflow fan-out in `workflows/hcvtyper.nf` (L406–444) iterates **all** candidate rows already, but picks the per-rank FASTA from the legacy `major_mapping`/`minor_mapping` emits with `rank == '1' ? major_fasta : minor_fasta` — a hard two-slot ceiling. The N-slot generalization (D-01) must replace this with an N-FASTA emit from `PARSEFIRSTMAPPING` and a rank-indexed lookup, then delete the `n_candidates > 2` guard at L85–87.

**Primary recommendation:** Plan the filename-slot migration (config closures + summarize.R candidate_rank join) and the summarize.R parse cutover as one lockstep wave, regenerate the two nf-test `.snap` files in the same wave, then add `bin/tests/test_compat.R` as a separate wave. Run the R suite via the pinned Seqera Docker container (tidyverse is not installed on the host); regenerate nf-test snapshots via the `nf-test` binary in the `NEXTFLOW` conda env.

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Output-filename slot naming (`.cand{N}.`) | Nextflow config (`conf/modules_hcv.config` ext.prefix) | Nextflow workflow (channel emit/fan-out) | The slot string is composed entirely in `ext.prefix` closures; the workflow only routes FASTAs per rank |
| N-ranked FASTA emission | Nextflow module (`parsefirstmapping/main.nf`) + R selection script | Nextflow workflow fan-out | The R script writes the per-rank FASTAs; the module declares the output glob; the workflow distributes them |
| Filename → candidate parsing | R summary (`bin/summarize.R`) | candidates CSV (data join key) | Parsing moves from filename position-3 string to a `candidate_rank` data join |
| Legacy `Major_*`/`Minor_*` column presence | R summary (`bin/summarize.R`) | — | Columns are emitted in summarize.R; Phase 9 only validates them, no structural change (D-05) |
| Golden-baseline reproduction | R regression test (`bin/tests/test_compat.R`) | R summary (system under test) | Test invokes/sources summarize.R logic and asserts core strain-call columns |
| COMPAT-04 exception preservation | R classifier (`bin/classify_roles.R` already wired) | R regression test (integration assertion) | Logic shipped in Phase 8; Phase 9 adds an end-to-end test through summarize.R |
| CI regression execution | GitHub Actions (`.github/workflows/ci.yml` r-regression job) | `run_all.sh` glob | New test auto-discovered by `test_*.R` glob; no CI YAML change needed |

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

**Filename slot migration (COMPAT-02)**
- **D-01 (Full N-slot generalization):** Phase 9 migrates `.major.`/`.minor.` → `.cand1.`/`.cand2.` AND generalizes parsefirstmapping to emit N ranked FASTAs (one per candidate), lifting the n_candidates > 2 guard. COMPAT-02 fully satisfied. The hcvtyper.nf:83 comment ("N>2 support arrives with Phase 9") is fulfilled.
- **D-02 (candidate_rank join in summarize.R):** The three stats-parsing blocks in `summarize.R` that extract `first_major_minor` from filename position 3 (`str_split(basename(...), "\\.")[[1]][3]`) are refactored to join on `candidate_rank` from the per-sample candidates CSV instead. This is the canonical, fragility-free approach — the candidates CSV already carries `candidate_rank` and is already used downstream. All `first_major_minor == "major"/"minor"` logic is replaced with `candidate_rank == 1` / `candidate_rank == 2` etc.
- **D-03 (module snapshot update):** nf-test `.snap` files for `parsefirstmapping`, `blastparse`, and any other affected modules are regenerated in Phase 9 to reflect the renamed outputs. Module tests must not hardcode `.major.fa`/`.minor.fa` after this phase.
- **D-04 (Nextflow channel names):** The `major_mapping`/`minor_mapping` channel emit names in parsefirstmapping are renamed to a per-rank or generic pattern (e.g., `candidate_fastas` or `cand_mapping_N`). The workflow `PARSEFIRSTMAPPING.out.major_mapping` / `.out.minor_mapping` join logic in `hcvtyper.nf` is generalized to collect all N ranked FASTAs.

**Legacy column aliasing (COMPAT-03)**
- **D-05 (validate-only, no structural change):** Phase 8 already retains `Major_*`/`Minor_*` mapping-stat columns alongside `Major_role_*`/`Minor_role_*`. Phase 9 adds a regression test asserting these legacy columns exist in Summary.csv output — no new column additions or renames in `summarize.R`.

**COMPAT-01 golden baseline test**
- **D-06 (R unit test, core columns only):** Add a `test_compat.R` case that runs `summarize.R` on an in-memory fixture and asserts that the core strain-call columns match expected v1.0/v2.0 values: `Major_reference`, `Minor_reference`, `Major_genotype_mapping`, `Minor_genotype_mapping`, `overall_sample_call`. New additive columns (role, dominance_score, etc.) are NOT part of the comparison.
- **D-07 (fixture scope):** Two golden cases are sufficient for COMPAT-01: (a) a simple monoinfection sample (1a dominant, no minor → `monoinfection`, `Major_reference` = 1a ref, no Minor_reference), and (b) a genuine co-infection sample (1a dominant + 1b co-infection confirmed by assembly support → `co-infection`, both references populated). These should reproduce the same strain-call output that the v1.0/v2.0 pipeline would have given.

**Regression suite structure (TEST-01)**
- **D-08 (one new test_compat.R):** A single `bin/tests/test_compat.R` covers all four COMPAT requirements:
  - COMPAT-01: golden baseline strain-call column assertion (D-06/D-07)
  - COMPAT-02: summarize.R stats-parsing correctly handles `cand1`/`cand2` filenames via candidate_rank join (smoke-test with synthetic stats files)
  - COMPAT-03: Summary.csv output contains legacy `Major_*`/`Minor_*` columns alongside `Major_role_*`/`Minor_role_*` columns
  - COMPAT-04: 1a/1b co-infection allowed, 2k/1b pair blocked — end-to-end integration case via `summarize.R` invocation
  Picked up automatically by `run_all.sh` glob; wired into the existing `r-regression` CI job (`.github/workflows/ci.yml:45`).

### Claude's Discretion
- Exact new channel emit names for the N-FASTA fan-out in parsefirstmapping (follow nf-core conventions: lowercase with underscores, e.g., `candidate_fasta`).
- Whether the N-FASTA emit is a tuple-per-candidate or a collected list — pick the pattern that minimizes change to the TARGETED_MAPPING fan-out logic in hcvtyper.nf.
- Exact `test_compat.R` fixture structure (build inline synthetic data or reuse the `flagoff_golden.csv` fixture as the v1.0 reference anchor).
- Whether `test_compat.R` invokes `summarize.R` via `system2("Rscript", ...)` (subprocess pattern from test_candidate_selection.R) or sources the helpers directly (function pattern from test_classify_roles.R) — pick the lower-friction approach given the candidate_rank refactor.

### Deferred Ideas (OUT OF SCOPE)
- **Dropping `Major_*`/`Minor_*` legacy columns:** deferred to the next milestone post-v3.0. Phase 9 only validates their presence.
- **Renaming `major_mapping`/`minor_mapping` channel emit identifiers to role-semantic names:** if the N-slot generalization makes these obsolete, the cleanup can be its own tidy-up commit rather than a Phase 9 hard requirement.
- **Cross-platform Singularity/Conda nf-test profiles:** not in Phase 9 scope; existing CI runs Docker only.
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| COMPAT-01 | Non-breaking — default flags + candidate count 2 reproduce v1.0/v2.0 reporting | `test_compat.R` golden cases (D-06/D-07); core strain-call columns `Major_reference`/`Minor_reference`/`Major_genotype_mapping`/`Minor_genotype_mapping`/`overall_sample_call` all assembled in summarize.R; `flagoff_golden.csv` available as anchor |
| COMPAT-02 | Output-filename slot `.major.`/`.minor.` → `.cand1.`/`.cand2.` with summarize.R parsing in lockstep | Slot produced by 6 `ext.prefix` closures in `conf/modules_hcv.config` (already have a `cand${rank}` arm); 3 position-3 parses + cv_by_ref `_(major|minor)$` strip in summarize.R; N-FASTA fan-out in `hcvtyper.nf` L406–444; guard at L85–87 to delete |
| COMPAT-03 | Legacy `Major_*`/`Minor_*` columns aliased one release alongside role columns | Phase 8 already emits both; validate-only test (D-05); `Major_role_*`/`Minor_role_*` at summarize.R L695–716, legacy `Major_*`/`Minor_*` at L282–488 |
| COMPAT-04 | 1a/1b co-infection + 2k/1b recombinant exceptions preserved | `is_valid_minor()` ported into `classify_roles.R` (Phase 8 D-12); already unit-tested in `test_classify_roles.R`; Phase 9 adds end-to-end integration case |
| TEST-01 | Regression suite extended; wired into CI | `run_all.sh` glob auto-discovers `test_*.R`; `r-regression` CI job runs `bash bin/tests/run_all.sh` in the pinned Seqera container; no CI YAML edit needed |
</phase_requirements>

## Standard Stack

This phase introduces **no new external dependencies**. It edits existing Nextflow config/modules/workflow and R scripts, and adds one R test file. The toolchain is already pinned and in use.

### Core
| Tool | Version | Purpose | Why Standard |
|------|---------|---------|--------------|
| Nextflow | 24.04.0 + latest-everything (CI matrix) | Pipeline DSL; `ext.prefix` closures own the filename slot | Already the pipeline runtime [VERIFIED: .github/workflows/ci.yml] |
| R + tidyverse | r-base 4.x, tidyverse (Seqera container tag `5358395134867368`) | `summarize.R` refactor + regression tests | Pinned Seqera container reused for both PARSEFIRSTMAPPING and the r-regression CI job [VERIFIED: conf modules + ci.yml] |
| nf-test | 0.9.2 | Module snapshot regeneration | Present in the `NEXTFLOW` conda env [VERIFIED: /home/jonbra/miniforge3/envs/NEXTFLOW/.../nf-test] |

### Supporting
| Tool | Version | Purpose | When to Use |
|------|---------|---------|-------------|
| Docker | 28.4.0 | Runs the pinned R container for the regression suite | Always (tidyverse not installed on host — see Environment Availability) |
| `run_all.sh` | repo script | Glob-runs every `bin/tests/test_*.R`, accumulates failures | Local dev + CI; new test auto-discovered |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| candidate_rank join (D-02) | Keep position-3 filename parse, just rename the slot string | Rejected by D-02 — filename-position parsing is the fragility this phase removes; a join on the already-loaded `candidates_long` is canonical |
| Subprocess test (`system2`) | Function-sourcing test | Discretion item; both patterns exist in the suite. summarize.R is a monolithic script (not sourceable helpers), so COMPAT-01/02/03 likely favor the subprocess pattern, COMPAT-04 the function pattern (mirrors test_classify_roles.R) |

**Installation:** None. No `npm`/`pip`/`cargo` install. All tooling present.

## Package Legitimacy Audit

Not applicable — this phase installs **no external packages**. All edits are to existing repo files; all tooling (Nextflow, the pinned R/tidyverse Seqera container, nf-test, Docker) is already provisioned. No legitimacy gate required.

## Architecture Patterns

### System Architecture Diagram

```
                        params.n_candidates (default 2)
                                   │
   reads ──► first-pass mapping ──► PARSEFIRSTMAPPING ──┬─► *.candidates.csv (rank 1..N, candidate_rank col)
                                  (summarize_mapping_     │
                                   to_all_references.R)   ├─► N ranked FASTAs  ◄── D-01: currently 2 slots
                                                          │   (major.fa/minor.fa → cand{rank}.fa)
                                                          │
                                                          ▼
        hcvtyper.nf fan-out (L406-444): splitCsv ALL rows, one element/candidate
            └── pick per-rank FASTA  ◄── D-01: rank==1?major:minor  →  rank-indexed N lookup
                                                          │
                                                          ▼
                                   TARGETED_MAPPING (per candidate)
                                                          │
        ext.prefix closures (conf/modules_hcv.config)  ◄── D-04/COMPAT-02: the SLOT is composed here
            "${id}.${reference}.${rank==1?major:rank==2?minor:cand{rank}}.{withdup|nodup}"
                                                          │
                          ┌───────────────────────────────┼────────────────────────────┐
                          ▼                                ▼                             ▼
                  *.<ref>.<slot>.withdup.stats   *.<ref>.<slot>.nodup.stats   *.<ref>.<slot>.nodup.tsv (depth)
                          │                                │                             │
                          └────────────────► bin/summarize.R ◄────────────────────────┘
                              parse blocks (L250-306, L314-363, L381-492):
                              str_split(basename, ".")[[1]][3] → first_major_minor
                                  ◄── D-02: REPLACE with join on candidate_rank
                                            (candidates_long already loaded L562-594)
                                                          │
                                                          ▼
                              Summary.csv: legacy Major_*/Minor_* (COMPAT-03 validate)
                                           + Major_role_*/Minor_role_* + overall_sample_call
```

File-to-implementation mapping is in the Component Responsibilities below; the diagram shows where the slot string is born (config) and where it is consumed (summarize.R).

### Recommended Project Structure (touchpoints, not new dirs)
```
conf/modules_hcv.config          # 6 ext.prefix closures: rank-1/2 → major/minor special-case is the slot source (COMPAT-02)
workflows/hcvtyper.nf            # L85-87 guard (delete), L406-444 N-FASTA fan-out (generalize)
modules/local/parsefirstmapping/main.nf   # emit N ranked FASTAs + rename channels (D-04); update stub
modules/local/blastparse/main.nf          # rename *.major.fa/*.minor.fa stubs (D-03)
bin/summarize.R                  # 3 position-3 parses + cv_by_ref strip → candidate_rank join (D-02)
bin/tests/test_compat.R          # NEW — all four COMPAT cases (D-08)
bin/tests/fixtures/              # optional new fixtures or reuse flagoff_golden.csv
modules/local/parsefirstmapping/tests/main.nf.test.snap   # regenerate (D-03)
modules/local/blastparse/tests/main.nf.test.snap          # regenerate (D-03)
```

### Pattern 1: Slot-from-rank ternary (the migration hinge)
**What:** The `.major.`/`.minor.`/`.cand{N}.` slot is the 3rd dot-field, composed in `ext.prefix` closures that already contain a full N-arm ternary.
**When to use:** All six TARGETED_MAPPING closures (SAMTOOLS_SORMADUP, SAMTOOLS_DEPTH, STATS_WITHDUP, STATS_MARKDUP, IVAR_CONSENSUS, CONSENSUS_DISTANCE).
**Current (two-slot special-case):**
```groovy
// Source: conf/modules_hcv.config L131 [VERIFIED: conf/modules_hcv.config]
ext.prefix = { "${meta.id}.${meta.reference}.${meta.candidate_rank.toInteger() == 1 ? 'major' : (meta.candidate_rank.toInteger() == 2 ? 'minor' : "cand${meta.candidate_rank}")}.nodup" }
```
**Migration target (uniform cand-slot):**
```groovy
// COMPAT-02: every rank → cand{rank}; drop the major/minor special-case
ext.prefix = { "${meta.id}.${meta.reference}.cand${meta.candidate_rank.toInteger()}.nodup" }
```
**Caution:** SAMTOOLS_DEPTH uses the `meta1.` namespace (joined input), not `meta.` — preserve that. [VERIFIED: conf/modules_hcv.config L224-228]

### Pattern 2: candidate_rank join replacing filename position-3 parse (D-02)
**What:** Each of the three stats loops sets `first_major_minor` from `str_split(basename, "\\.")[[1]][3]`, then drives every `Major_*`/`Minor_*` column off `first_major_minor == "major"/"minor"`.
**Current:**
```r
# Source: bin/summarize.R L261/L326/L393 [VERIFIED: bin/summarize.R]
tmp_df$first_major_minor[i] <- str_split(basename(stats_files[i]), "\\.")[[1]][3]
# ... mutate(Major_reference = case_when(first_major_minor == "major" ~ reference)) ...
```
**Migration target:** parse the slot to a rank integer (`cand1` → 1), or — preferred per D-02 — join each loop's per-file rows to `candidates_long` on `(sampleName, candidate_ref)` and read `candidate_rank`, then map `candidate_rank == 1 → Major_*`, `candidate_rank == 2 → Minor_*`. `candidates_long` is already read at L562-594 with `candidate_rank = col_integer()`. The reference token in the cov loop already matches `candidate_ref` after the `_(major|minor)$` strip (L464) — after migration the strip becomes `_cand[0-9]+$` (or is removed if the slot moves out of the reference field).
**Key constraint:** `Major_reference`/`Minor_reference` are a JOIN KEY downstream:
```r
# Source: bin/summarize.R L1043 [VERIFIED]
left_join(df_coverage, join_by(sampleName, Major_reference, Minor_reference))
```
The refactor must keep `Major_reference`/`Minor_reference` byte-identical in value (the cleaned reference name, no slot suffix) so this join still matches. This is the highest-risk silent-row-drop site.

### Pattern 3: Glob-discovered R regression test (D-08, TEST-01)
**What:** `run_all.sh` runs `for t in "$here"/test_*.R` and accumulates exit codes; CI runs that script in the pinned container. A new `test_compat.R` needs zero CI wiring.
**Test entry-point boilerplate (both existing patterns use it):**
```r
# Source: bin/tests/test_candidate_selection.R L46-51 [VERIFIED]
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
this_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
bin_dir  <- normalizePath(file.path(this_dir, ".."))
```

### Anti-Patterns to Avoid
- **Renaming the slot in config but not summarize.R (or vice-versa):** the lockstep failure mode — a mismatch silently drops candidate rows (no crash). This is the documented Phase 9 highest-risk concern [VERIFIED: STATE.md L124]. Land both in the same wave and prove with a full run + the COMPAT-02 test.
- **Coercing `candidate_rank`/`confirmation_status` to Integer in Groovy:** a single/no-candidate NA field crashes `NA.toInteger()`. The fan-out deliberately keeps these as Strings [VERIFIED: hcvtyper.nf L416-444 Pitfall 3 comment]. Preserve that when generalizing.
- **Forgetting the `cv_by_ref` `_(major|minor)$` strip** (summarize.R L464): it is a fourth slot-coupled site beyond the three position-3 parses. Miss it and per-candidate evenness silently fails to join → dominance scores shift → COMPAT-01 golden breaks.
- **Hand-editing `.snap` md5 values:** always regenerate via `nf-test --update-snapshot` against the renamed outputs (D-03).
- **Adding/removing legacy columns in COMPAT-03:** D-05 is validate-only; any structural change is out of scope.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Recovering candidate rank in summarize.R | A bespoke filename regex parser | Join on `candidate_rank` in the already-loaded `candidates_long` (D-02) | The CSV carries the integer rank authoritatively; filename parsing is the exact fragility being removed |
| nf-test snapshot md5s | Manual md5 computation/editing | `nf-test test <module> --update-snapshot` | Snapshots are tool-generated; hand-edits drift and break CI |
| CI wiring for the new test | New job/step in ci.yml | `run_all.sh` glob auto-discovery | The glob already runs all `test_*.R`; the r-regression job already invokes it |
| v1.0 reference values for COMPAT-01 | Re-derive expected numbers from scratch | Reuse `bin/tests/fixtures/flagoff_golden.csv` as the anchor (discretion) | Existing golden fixture already encodes the legacy column shape |
| Running R locally | Installing tidyverse on host | The pinned Seqera container via Docker | tidyverse is NOT on the host; CI uses the container — match it for parity |

**Key insight:** Phase 9 is deletion-and-join, not construction. The N-arm ternary, the candidates CSV with `candidate_rank`, the role classifier with COMPAT-04 exceptions, and the glob-discovered test runner already exist. The work is removing the two-slot special-cases and proving nothing regressed.

## Runtime State Inventory

This is a refactor/migration phase (filename-slot rename), so the inventory applies. The migration is **code/config only** — no persistent datastore, live-service, or OS-registered state embeds the `.major.`/`.minor.` slot.

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| Stored data | None — the `.major.`/`.minor.` slot appears only in transient Nextflow work-dir intermediate filenames consumed within a single run; no database, collection, or persisted key uses it. Verified by grepping config/bin/workflows for the slot and finding only `ext.prefix` closures + summarize.R parses. | None |
| Live service config | None — no external service (GLUE, Datadog, etc.) stores the candidate slot. NOTE: GLUE produces `GLUE_collected_report_major.tsv`/`_minor.tsv` [VERIFIED: bin/GLUE_json_parser.R L260,L263], but that is the GLUE **role-aggregation** naming, independent of the targeted-mapping intermediate slot, and is OUT OF SCOPE for COMPAT-02 (not a candidate filename). Confirm the plan does NOT touch it. | None (do not rename GLUE reports) |
| OS-registered state | None — no Task Scheduler / systemd / pm2 registration references the slot. | None |
| Secrets/env vars | None — no secret or env var name references major/minor candidate slots. | None |
| Build artifacts | nf-test `.snap` files for `parsefirstmapping` and `blastparse` contain hardcoded `*.major.fa`/`*.minor.fa` md5 entries [VERIFIED: modules/local/*/tests/*.snap] — these are committed test artifacts that go stale after the rename. | Regenerate both `.snap` files with `nf-test --update-snapshot` (D-03) |

**Canonical question — after every file is updated, what runtime systems still have the old string cached?** Answer: only the two committed `.snap` golden files (build artifacts), which are regenerated in-phase. No live or stored state survives the rename. The GLUE `_major.tsv`/`_minor.tsv` naming is a distinct concern and intentionally untouched.

## Common Pitfalls

### Pitfall 1: Lockstep mismatch silently drops rows
**What goes wrong:** Config emits `.cand1.` but summarize.R still parses position-3 expecting `major`, so `first_major_minor` becomes `"cand1"`, every `case_when(first_major_minor == "major" ~ ...)` returns NA, and the sample's Major_* columns silently empty. No error.
**Why it happens:** The writer (config) and the parser (summarize.R) are in different files and can be changed independently.
**How to avoid:** Land both changes in one wave; the COMPAT-02 test (D-08) must build synthetic `*.cand1.*.stats`/`*.cand2.*.stats` files and assert summarize.R yields populated Major_*/Minor_* rows.
**Warning signs:** Summary.csv rows present but Major_reference/Major_genotype_mapping all NA; the `join_by(..., Major_reference, Minor_reference)` at L1043 produces a Cartesian or empty join.

### Pitfall 2: Missed fourth slot-coupled site (cv_by_ref)
**What goes wrong:** The three position-3 parses get refactored but the `str_remove(reference, "_(major|minor)$")` at summarize.R L464 is overlooked, so per-candidate `cv_evenness` fails to join to `candidate_ref` → dominance scores change → COMPAT-01 golden mismatch on `overall_sample_call`.
**Why it happens:** This site strips the slot from a reference token rather than reading position-3, so a grep for `[[1]][3]` misses it.
**How to avoid:** Grep for both `[[1]][3]` AND `_(major|minor)` / `_major` / `_minor` across summarize.R; there are reference-cleanup `str_remove`s at L290-291, L356-357, L485-486 too.
**Warning signs:** `cv_by_ref` join produces all-NA `cv_evenness`; dominance scores differ from the Phase 8 `test_dominance_score.R` expectations.

### Pitfall 3: N-FASTA fan-out drops the dominant candidate
**What goes wrong:** Replacing the two-slot `rank == '1' ? major_fasta : minor_fasta` (hcvtyper.nf L431) with an N-lookup but losing the `remainder: true` joins (L407-408) → a single-candidate sample (no `_minor.fa`) gets dropped entirely, including its passing major.
**Why it happens:** The legacy emits are `optional: true`; a strict join silently drops samples missing an optional emit.
**How to avoid:** Preserve `remainder: true` semantics in whatever N-FASTA join replaces the two legacy joins; the per-rank null-FASTA guard (`entry[1] != null`) must remain.
**Warning signs:** Single-candidate (monoinfection) samples vanish from Summary.csv; the COMPAT-01 monoinfection golden case (D-07a) fails.

### Pitfall 4: nf-test snapshot regen masks a real regression
**What goes wrong:** Blindly running `--update-snapshot` accepts whatever the modules now emit, including an unintended change.
**Why it happens:** Snapshot update is a blanket "accept current output."
**How to avoid:** Before updating, eyeball the snap diff — the ONLY expected change is `*.major.fa`/`*.minor.fa` → `*.cand1.fa`/`*.cand2.fa` (and the candidates.csv md5 if its content changed). Any other md5 churn is a real regression to investigate.
**Warning signs:** Snap diff touches `candidates.csv` md5 unexpectedly, or stats/depth filenames change shape beyond the slot.

### Pitfall 5: SUMMARIZE positional ext.args mis-mapping
**What goes wrong:** Inserting any new arg into the SUMMARIZE `ext.args` string mid-list shifts every later positional index in summarize.R (args[9]/[10]/[11]... are read by index).
**Why it happens:** summarize.R reads args by fixed position [VERIFIED: bin/summarize.R L31-68], and the ext.args string at conf/modules_hcv.config L328 is positional.
**How to avoid:** Phase 9 should not need new SUMMARIZE args. If one is unavoidable, APPEND at the end (the documented Phase 8 T-08-04 rule).
**Warning signs:** `n_candidates` parses as a weight, role columns shift.

## Code Examples

### Loading candidate_rank in summarize.R (already present — reuse for the D-02 join)
```r
# Source: bin/summarize.R L570-580 [VERIFIED]
candidates_long <- map_dfr(candidates_files, ~ read_csv(.x, col_types = cols(
  sample              = col_character(),
  candidate_rank      = col_integer(),
  candidate_ref       = col_character(),
  candidate_subtype   = col_character(),
  candidate_genotype  = col_character(),
  candidate_reads     = col_double(),
  candidate_cov       = col_double(),
  confirmation_status = col_character()
))) %>%
  rename(sampleName = sample)
```

### Subprocess test pattern (model for COMPAT-01/02/03 — system under test is a whole script)
```r
# Source: bin/tests/test_candidate_selection.R L96-108 [VERIFIED]
old <- getwd(); setwd(wd); on.exit(setwd(old), add = TRUE)
exit <- system2("Rscript",
  c(shQuote(script), shQuote(idx_path), shQuote(depth_path),
    shQuote(sampleName), shQuote(refs_path), minRead, minCov, n_candidates),
  stdout = FALSE, stderr = FALSE)
out <- if (file.exists(out_path)) read_csv(out_path, show_col_types = FALSE) else NULL
```

### Function-sourcing test pattern (model for COMPAT-04 — exercises classify_roles end-to-end)
```r
# Source: bin/tests/test_classify_roles.R L63-67 [VERIFIED]
classify <- function(df, minRead = 500, minCov = 30) {
  classify_roles(score_candidates(df), minRead = minRead, minCov = minCov,
                 denovo_min_contig_length = 1000, denovo_min_kmer_cov = 2.0,
                 denovo_min_blast_identity = 90, match_level = "genotype")
}
```

### Running the R suite the way CI does (host has no tidyverse)
```bash
# Source: .github/workflows/ci.yml r-regression job [VERIFIED]
docker run --rm -v "$PWD":/work -w /work \
  community.wave.seqera.io/library/r-seqinr_r-tidyverse:5358395134867368 \
  bash bin/tests/run_all.sh
```

### Regenerating module snapshots (D-03)
```bash
# nf-test is in the NEXTFLOW conda env [VERIFIED: /home/jonbra/miniforge3/envs/NEXTFLOW/...]
nf-test test modules/local/parsefirstmapping/tests/main.nf.test --profile test,docker --update-snapshot
nf-test test modules/local/blastparse/tests/main.nf.test --profile test,docker --update-snapshot
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Filename position-3 string drives Major/Minor logic | Join on `candidate_rank` data column (D-02) | Phase 9 | Removes the fragile filename-coupling; enables N>2 |
| Fixed two FASTA slots (major.fa/minor.fa) | N ranked FASTAs, cand{rank}.fa (D-01) | Phase 9 | Lifts the `n_candidates > 2` guard (hcvtyper.nf L85-87) |
| `apply_denovo_layer()` confirmation path | `score_candidates()` + `classify_roles()` role model | Phase 8 (D-15) | Already shipped; Phase 9 only tests it |
| Asymmetric MAJOR/MINOR_MAPPING aliases | One per-candidate TARGETED_MAPPING fan-out | Phase 6 (06-03) | Already shipped; Phase 9 generalizes the FASTA pick to N |

**Deprecated/outdated:**
- The two-slot major/minor special-case in the 6 `ext.prefix` closures: replaced by uniform `cand{rank}` (COMPAT-02).
- `n_candidates > 2` guard at hcvtyper.nf L85-87: deleted in-phase once D-01 lands.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | The GLUE `_major.tsv`/`_minor.tsv` naming is out of scope for COMPAT-02 (it is role-aggregation output, not a candidate intermediate slot). | Runtime State Inventory | If a stakeholder considers GLUE report renaming part of COMPAT-02, scope expands; mitigated by 09-CONTEXT explicitly scoping COMPAT-02 to parsefirstmapping/blastparse/TARGETED_MAPPING slots only |
| A2 | Reusing `flagoff_golden.csv` values reproduces the true v1.0/v2.0 strain-call output for the COMPAT-01 cases. | COMPAT-01 / Don't Hand-Roll | If the golden fixture predates a behavior change, COMPAT-01 asserts the wrong baseline; planner should add a checkpoint to confirm the two D-07 fixture cases against a known v1.0 run or the handoff evidence table |
| A3 | No module other than parsefirstmapping and blastparse needs a `.snap` regen for the slot rename. | D-03 | A workflow-level (`tests/`) snapshot or another module could reference the slot; planner should run the full nf-test suite once and regen any other snap that diffs only on the slot |
| A4 | summarize.R remains a monolithic non-sourceable script (no extractable helper for the stats loops), so COMPAT-01/02/03 tests use the subprocess pattern. | Standard Stack / Alternatives | If the refactor extracts the parse into a sourceable helper, the function-pattern becomes viable and simpler; this is a discretion item |

## Open Questions

1. **Do the two COMPAT-01 golden cases need a real pipeline run to anchor expected values, or is the `flagoff_golden.csv` fixture authoritative?**
   - What we know: `flagoff_golden.csv` encodes Major_reference/Minor_reference/typable for two samples (S1 co-infection, S2 monoinfection) [VERIFIED: bin/tests/fixtures/flagoff_golden.csv].
   - What's unclear: whether those exact values equal what v1.0/v2.0 would emit for the D-07 1a/1b and 1a-mono cases, including `Major_genotype_mapping`/`overall_sample_call`.
   - Recommendation: planner adds a `checkpoint:human-verify` to confirm the fixture values against the handoff evidence table or a v1.0 reference run before locking the assertions (A2).

2. **Tuple-per-candidate vs collected-list for the N-FASTA emit (discretion D-04).**
   - What we know: the current fan-out joins two optional emits with `remainder: true` then picks by rank (hcvtyper.nf L407-431).
   - What's unclear: which emit shape minimizes change to the `flatMap` rank lookup.
   - Recommendation: a single `candidate_fasta` tuple emit `[meta(with candidate_rank), fasta]` joined to the candidates fan-out by full meta is the lowest-churn option — it lets the existing `flatMap` keep iterating candidate rows and attach the matching FASTA by rank, removing the major/minor branch entirely. Preserve the null-FASTA guard.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| Rscript (host) | Quick local R syntax checks | ✓ | 4.3.3 | — |
| R tidyverse (host) | Running the regression suite locally | ✗ | — | Run via the pinned Seqera Docker container (CI parity) |
| Docker | Running R suite + nf-test docker profile | ✓ | 28.4.0 | — |
| nf-test | Module snapshot regeneration (D-03) | ✓ (conda env `NEXTFLOW`) | 0.9.2 | — |
| nextflow (host bin) | Full pipeline run for COMPAT-01/04 end-to-end proof | ✗ (permission denied on `/usr/local/bin/nextflow`) | — | Use the `NEXTFLOW` conda env's nextflow, or run nf-test which bundles execution |

**Missing dependencies with no fallback:** None.
**Missing dependencies with fallback:**
- Host tidyverse → use the pinned `community.wave.seqera.io/library/r-seqinr_r-tidyverse:5358395134867368` container (this is exactly what CI does; use it for parity).
- Host `nextflow` binary has a permission error → use the conda `NEXTFLOW` env (which also provides nf-test 0.9.2) for any full-run or snapshot work.

## Validation Architecture

> `.planning/config.json` was not found at the standard path; nyquist_validation key absent → treated as ENABLED.

### Test Framework
| Property | Value |
|----------|-------|
| Framework (R) | base R asserts via `fail()`/`ok()` helpers + `run_all.sh` accumulator (custom, no testthat) |
| Framework (Nextflow) | nf-test 0.9.2 |
| Config file | `nf-test.config` (root); `bin/tests/run_all.sh` for R |
| Quick run command | `Rscript bin/tests/test_compat.R` (inside the R container) |
| Full suite command | `docker run --rm -v "$PWD":/work -w /work community.wave.seqera.io/library/r-seqinr_r-tidyverse:5358395134867368 bash bin/tests/run_all.sh` |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| COMPAT-01 | core strain-call columns reproduce v1.0/v2.0 | unit (subprocess) | `Rscript bin/tests/test_compat.R` | ❌ Wave 0 (new file) |
| COMPAT-02 | summarize.R parses `cand1`/`cand2` via candidate_rank join | unit (subprocess, synthetic stats) | `Rscript bin/tests/test_compat.R` | ❌ Wave 0 |
| COMPAT-03 | Summary.csv has legacy Major_*/Minor_* alongside role columns | unit (column-presence assert) | `Rscript bin/tests/test_compat.R` | ❌ Wave 0 |
| COMPAT-04 | 1a/1b allowed, 2k/1b blocked (end-to-end) | integration (function or subprocess) | `Rscript bin/tests/test_compat.R` | ❌ Wave 0 (helper-level exists in test_classify_roles.R) |
| COMPAT-02 | renamed module outputs (cand{N}.fa) | nf-test snapshot | `nf-test test modules/local/parsefirstmapping/tests/main.nf.test --profile test,docker` | ✅ (regenerate snap) |

### Sampling Rate
- **Per task commit:** `Rscript bin/tests/test_compat.R` (in container) for R changes; `nf-test test <changed module>` for module changes.
- **Per wave merge:** full R suite via `run_all.sh` in container + `nf-test test` for both touched modules.
- **Phase gate:** full R suite green + both module snapshots green + (ideally) one full `nextflow run -profile test,docker` confirming Summary.csv parses correctly with the new slots.

### Wave 0 Gaps
- [ ] `bin/tests/test_compat.R` — new file covering COMPAT-01..04 (D-08)
- [ ] (optional) `bin/tests/fixtures/` additions for the two D-07 golden cases, OR reuse `flagoff_golden.csv`
- [ ] Regenerate `modules/local/parsefirstmapping/tests/main.nf.test.snap`
- [ ] Regenerate `modules/local/blastparse/tests/main.nf.test.snap`
- [ ] Framework install: none — R container + nf-test already provisioned

## Security Domain

> `security_enforcement` config not located (absent = enabled). This phase processes no untrusted external input, handles no auth/secrets/PII, and adds no network surface — it edits an internal bioinformatics pipeline's filename conventions and test code.

### Applicable ASVS Categories

| ASVS Category | Applies | Standard Control |
|---------------|---------|-----------------|
| V2 Authentication | no | Phase touches no auth |
| V3 Session Management | no | No sessions |
| V4 Access Control | no | No access-control surface |
| V5 Input Validation | minimal | `candidate_rank.toInteger()` could throw on a malformed/NA field — already mitigated by keeping it a String in the fan-out and only `.toInteger()`-ing inside config closures where the value is known-present (Pitfall 2 from prior phases) |
| V6 Cryptography | no | No crypto; the only md5s are nf-test snapshot fixtures (integrity, not security) |

### Known Threat Patterns for {Nextflow + R pipeline}

| Pattern | STRIDE | Standard Mitigation |
|---------|--------|---------------------|
| Malformed `candidate_rank` (NA) crashing `.toInteger()` | Denial of Service (run abort) | Keep rank/confirmation_status as Strings in Groovy fan-out; only coerce where guaranteed non-null [VERIFIED: hcvtyper.nf L416-444] |
| Filename-injection via sample name | Tampering | Out of scope — sample names are operator-controlled samplesheet inputs, unchanged by this phase |

## Sources

### Primary (HIGH confidence)
- `conf/modules_hcv.config` (L115-264, L328) — the six `ext.prefix` slot closures + SUMMARIZE positional ext.args [VERIFIED]
- `workflows/hcvtyper.nf` (L85-87 guard, L380-448 fan-out) — N-FASTA distribution + guard to delete [VERIFIED]
- `bin/summarize.R` (L250-492 parse blocks, L562-594 candidates load, L687-738 role columns, L1031-1359 final/Summary.csv) [VERIFIED]
- `modules/local/parsefirstmapping/main.nf`, `modules/local/blastparse/main.nf` — output emits + stubs [VERIFIED]
- `subworkflows/local/targeted_mapping/main.nf` — reference enrichment, per-candidate mapping [VERIFIED]
- `bin/tests/test_candidate_selection.R`, `bin/tests/test_classify_roles.R`, `bin/tests/run_all.sh` — test patterns [VERIFIED]
- `.github/workflows/ci.yml` (r-regression job) — CI invocation + pinned container [VERIFIED]
- `modules/local/{parsefirstmapping,blastparse}/tests/main.nf.test.snap` — stale slot md5s [VERIFIED]
- `bin/tests/fixtures/flagoff_golden.csv` — golden anchor [VERIFIED]
- `.planning/phases/09-.../09-CONTEXT.md`, `.planning/REQUIREMENTS.md`, `.planning/STATE.md` — scope + decisions [VERIFIED]

### Secondary (MEDIUM confidence)
- Local environment probes (Rscript 4.3.3, no host tidyverse, Docker 28.4.0, nf-test 0.9.2 in conda env, nextflow host-bin permission error) [VERIFIED via shell]

### Tertiary (LOW confidence)
- None. No external/WebSearch sources were needed — this is an internal-integration phase.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — no new deps; all tooling verified present or container-pinned.
- Architecture: HIGH — every touchpoint and the slot-composition site verified in source.
- Pitfalls: HIGH — derived from in-code comments documenting prior-phase failure modes (Pitfalls 1-3,5) and direct code reading (Pitfall 2/4).
- COMPAT-01 baseline values: MEDIUM — fixture reuse is an assumption (A2); flagged for human-verify.

**Research date:** 2026-06-14
**Valid until:** 2026-07-14 (stable internal codebase; re-check if Phase 8 follow-ups alter summarize.R column layout or the candidates CSV schema)
