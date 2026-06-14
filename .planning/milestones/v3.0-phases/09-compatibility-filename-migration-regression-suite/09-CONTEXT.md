# Phase 9: Compatibility, Filename Migration + Regression Suite - Context

**Gathered:** 2026-06-14
**Status:** Ready for planning

<domain>
## Phase Boundary

Integrate and validate all v3.0 changes so the pipeline ships non-breaking. Four parallel tracks:

1. **COMPAT-02 — Filename slot migration + summarize.R parsing cutover:** rename the `.major.`/`.minor.` intermediate-file slots to `.cand1.`/`.cand2.`/`.cand{N}.` across parsefirstmapping, blastparse, and the TARGETED_MAPPING subworkflow. Generalize parsefirstmapping to emit N ranked FASTAs (lifting the n_candidates > 2 guard). Update `summarize.R`'s stats-parsing blocks to join on `candidate_rank` from the candidates CSV rather than extracting a slot name from filename position 3. Update affected nf-test module-level `.snap` files in the same commit.

2. **COMPAT-03 — Legacy column validation:** the old `Major_*`/`Minor_*` mapping-stat columns are already retained in Phase 8 alongside the new `Major_role_*`/`Minor_role_*` columns. Phase 9 validates this with a test — no new structural change to `summarize.R`.

3. **COMPAT-01/04 — Non-breaking verification:** prove that default flags + candidate_count=2 reproduce v1.0/v2.0 core strain-call output, and that the 1a/1b co-infection exception + 2k/1b recombinant suppression still hold (already wired in `classify_roles.R`).

4. **TEST-01 — Regression suite extension:** add `bin/tests/test_compat.R` covering all four COMPAT requirements. Automatically picked up by `run_all.sh` and wired into the existing `r-regression` CI job.

**Explicitly NOT in this phase:**
- Dropping the `Major_*`/`Minor_*` legacy columns — that's the "one release later" cleanup explicitly deferred in ROADMAP.md.
- Any changes to the classification logic (roles, scores, thresholds) — those shipped in Phase 8.
- Finishing the paused v2.0 HCVGLUE refactor (Phase 5).

</domain>

<decisions>
## Implementation Decisions

### Filename slot migration (COMPAT-02)
- **D-01 (Full N-slot generalization):** Phase 9 migrates `.major.`/`.minor.` → `.cand1.`/`.cand2.` AND generalizes parsefirstmapping to emit N ranked FASTAs (one per candidate), lifting the n_candidates > 2 guard. COMPAT-02 fully satisfied. The hcvtyper.nf:83 comment ("N>2 support arrives with Phase 9") is fulfilled.
- **D-02 (candidate_rank join in summarize.R):** The three stats-parsing blocks in `summarize.R` that extract `first_major_minor` from filename position 3 (`str_split(basename(...), "\\.")[[1]][3]`) are refactored to join on `candidate_rank` from the per-sample candidates CSV instead. This is the canonical, fragility-free approach — the candidates CSV already carries `candidate_rank` and is already used downstream. All `first_major_minor == "major"/"minor"` logic is replaced with `candidate_rank == 1` / `candidate_rank == 2` etc.
- **D-03 (module snapshot update):** nf-test `.snap` files for `parsefirstmapping`, `blastparse`, and any other affected modules are regenerated in Phase 9 to reflect the renamed outputs. Module tests must not hardcode `.major.fa`/`.minor.fa` after this phase.
- **D-04 (Nextflow channel names):** The `major_mapping`/`minor_mapping` channel emit names in parsefirstmapping are renamed to a per-rank or generic pattern (e.g., `candidate_fastas` or `cand_mapping_N`). The workflow `PARSEFIRSTMAPPING.out.major_mapping` / `.out.minor_mapping` join logic in `hcvtyper.nf` is generalized to collect all N ranked FASTAs.

### Legacy column aliasing (COMPAT-03)
- **D-05 (validate-only, no structural change):** Phase 8 already retains `Major_*`/`Minor_*` mapping-stat columns alongside `Major_role_*`/`Minor_role_*`. Phase 9 adds a regression test asserting these legacy columns exist in Summary.csv output — no new column additions or renames in `summarize.R`.

### COMPAT-01 golden baseline test
- **D-06 (R unit test, core columns only):** Add a `test_compat.R` case that runs `summarize.R` on an in-memory fixture and asserts that the core strain-call columns match expected v1.0/v2.0 values: `Major_reference`, `Minor_reference`, `Major_genotype_mapping`, `Minor_genotype_mapping`, `overall_sample_call`. New additive columns (role, dominance_score, etc.) are NOT part of the comparison — they're additions, not changes. The fixture should match the simple non-co-infection and genuine co-infection cases from the handoff evidence table.
- **D-07 (fixture scope):** Two golden cases are sufficient for COMPAT-01: (a) a simple monoinfection sample (1a dominant, no minor → `monoinfection`, `Major_reference` = 1a ref, no Minor_reference), and (b) a genuine co-infection sample (1a dominant + 1b co-infection confirmed by assembly support → `co-infection`, both references populated). These should reproduce the same strain-call output that the v1.0/v2.0 pipeline would have given.

### Regression suite structure (TEST-01)
- **D-08 (one new test_compat.R):** A single `bin/tests/test_compat.R` covers all four COMPAT requirements:
  - COMPAT-01: golden baseline strain-call column assertion (D-06/D-07)
  - COMPAT-02: summarize.R stats-parsing correctly handles `cand1`/`cand2` filenames via candidate_rank join (smoke-test with synthetic stats files)
  - COMPAT-03: Summary.csv output contains legacy `Major_*`/`Minor_*` columns alongside `Major_role_*`/`Minor_role_*` columns
  - COMPAT-04: 1a/1b co-infection allowed, 2k/1b pair blocked — this is already asserted in `test_classify_roles.R`, so `test_compat.R` includes an end-to-end integration case (via `summarize.R` invocation) that exercises the full path rather than just the helper function
  Picked up automatically by `run_all.sh` glob; wired into the existing `r-regression` CI job (`.github/workflows/ci.yml:45`).

### Claude's Discretion
- Exact new channel emit names for the N-FASTA fan-out in parsefirstmapping (follow nf-core conventions: lowercase with underscores, e.g., `candidate_fasta`).
- Whether the N-FASTA emit is a tuple-per-candidate or a collected list — pick the pattern that minimizes change to the TARGETED_MAPPING fan-out logic in hcvtyper.nf.
- Exact `test_compat.R` fixture structure (build inline synthetic data or reuse the `flagoff_golden.csv` fixture as the v1.0 reference anchor).
- Whether `test_compat.R` invokes `summarize.R` via `system2("Rscript", ...)` (subprocess pattern from test_candidate_selection.R) or sources the helpers directly (function pattern from test_classify_roles.R) — pick the lower-friction approach given the candidate_rank refactor.

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Phase scope & requirements
- `.planning/ROADMAP.md` §"Phase 9: Compatibility, Filename Migration + Regression Suite" — goal, 5 success criteria, dependency on Phase 8.
- `.planning/REQUIREMENTS.md` — COMPAT-01..04 and TEST-01 (all Phase 9 pending requirements).
- `.planning/PROJECT.md` — "Non-breaking" constraint, COMPAT requirements in Active section, Key Decisions table.

### Upstream phase decisions this phase depends on
- `.planning/phases/08-dominance-scoring-strain-role-classification/08-CONTEXT.md` — D-12 (`is_valid_minor()` ported as post-scoring special-cases in classify_roles.R — COMPAT-04 already wired), D-15 (legacy apply_denovo_layer retired), D-16 (wide Summary + enriched long candidates.csv), and the deferred list explicitly handing filename migration + legacy column aliasing to Phase 9.
- `.planning/phases/08-dominance-scoring-strain-role-classification/08-02-SUMMARY.md` — what Phase 8 Plan 02 actually shipped (key-files, key-decisions, patterns-established). Specifically: `Major_role_*`/`Minor_role_*` added additive alongside `Major_*`/`Minor_*`; `classify_roles.R` staged; affects list cites Phase 9 explicitly.
- `.planning/phases/06-neutral-candidate-selection/06-CONTEXT.md` — D-06 (legacy major/minor FASTA shim retained until Phase 9), Phase 6 comment in hcvtyper.nf about N>2 guard.

### Code touchpoints (primary change sites)
- `modules/local/parsefirstmapping/main.nf` — emit N ranked FASTAs instead of fixed `_major.fa`/`_minor.fa`; rename `major_mapping`/`minor_mapping` channel outputs.
- `modules/local/blastparse/main.nf` — rename `${prefix}.major.fa`/`${prefix}.minor.fa` stubs to cand-slot names.
- `workflows/hcvtyper.nf` — generalize the `PARSEFIRSTMAPPING.out.major_mapping`/`.minor_mapping` join and fan-out; lift the n_candidates > 2 guard at L85.
- `bin/summarize.R` — refactor the three stats-parsing blocks (~L250-260, ~L314-326, ~L381-393) that extract `first_major_minor` from filename position 3. Replace with a `candidate_rank`-based join against the per-sample candidates CSV. All `first_major_minor == "major"/"minor"` wrangling is replaced by `candidate_rank == 1` / `candidate_rank >= 2` etc.
- `bin/tests/test_compat.R` — new test file (D-08).
- `modules/local/parsefirstmapping/tests/main.nf.test.snap` — regenerate after filename rename.
- `modules/local/blastparse/tests/main.nf.test.snap` — regenerate after filename rename.

### Reference implementations / patterns
- `bin/tests/test_candidate_selection.R` — subprocess pattern (`system2("Rscript", ...)`) with synthetic fixture in tempdir; model for COMPAT-02 test case.
- `bin/tests/test_classify_roles.R` — function-sourcing pattern; model for COMPAT-04 test case.
- `bin/tests/fixtures/flagoff_golden.csv` — existing golden fixture; inspect as anchor for the COMPAT-01 v1.0/v2.0 reference values.
- `.github/workflows/ci.yml` §"r-regression" (L45-60) — existing CI job that runs `bash bin/tests/run_all.sh`; no changes needed, new test auto-discovered by glob.

### Conventions
- `.planning/codebase/CONVENTIONS.md` — channel naming (lowercase + underscores), stub block structure, `task.ext.prefix` pattern, R naming conventions.
- `.planning/codebase/TESTING.md` — nf-test framework, stub mode, snapshot update commands.

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **`candidate_rank` in candidates CSV:** already present as an R-emitted integer column; the refactored `summarize.R` stats-parsing joins on this directly, replacing the fragile filename-position extraction.
- **`is_valid_minor()` in classify_roles.R:** already encodes COMPAT-04 exceptions verbatim; `test_compat.R` exercises the full path end-to-end through `summarize.R` for integration coverage.
- **`flagoff_golden.csv`** (`bin/tests/fixtures/`): the v1.0 flag-off reference fixture; inspect its column structure as the template for the COMPAT-01 golden fixture.
- **`run_all.sh` glob pattern** (`"$here"/test_*.R`): picks up `test_compat.R` automatically — no CI wiring change needed.

### Established Patterns
- **Subprocess-style fixture test** (`test_candidate_selection.R`): build synthetic files in a `tempdir`, invoke the real R script via `system2("Rscript", ...)`, assert on emitted CSV. Use this for the COMPAT-02 summarize.R stats-parsing test.
- **Function-sourcing test** (`test_classify_roles.R`): source helpers + build in-memory tibbles. Use for COMPAT-04 role assertion.
- **nf-test stub pattern** (`parsefirstmapping/main.nf` stub block): create deterministic dummy outputs with `touch`/`> file`. The N-FASTA stub must emit N files (one per candidate slot).

### Integration Points
- `PARSEFIRSTMAPPING.out.major_mapping` / `.out.minor_mapping` join in `hcvtyper.nf` (~L407-408): this is where the per-rank FASTA fan-out logic lives. Generalize to N-slot by iterating the candidates CSV rows and looking up the per-rank FASTA by `candidate_rank`.
- `str_split(basename(stats_files[i]), "\\.")[[1]][3]` in `summarize.R` (~L261, L326, L393): the three fragile filename-position extraction sites that D-02 replaces with a `candidate_rank` join.

</code_context>

<specifics>
## Specific Ideas

- **n_candidates > 2 guard error message** (hcvtyper.nf:85) is removed entirely after D-01 is implemented — Phase 9 is explicitly what the message was pointing to.
- **nf-test snapshot regen** is cited in Phase 8's `08-02-SUMMARY.md` affects list as a known Phase 9 task — confirming this is expected work.
- **Golden baseline strain-call columns for COMPAT-01** (D-06): `Major_reference`, `Minor_reference`, `Major_genotype_mapping`, `Minor_genotype_mapping`, `overall_sample_call`. These are the minimum meaningful columns. Do NOT include new Phase 8 role columns in the COMPAT-01 comparison.

</specifics>

<deferred>
## Deferred Ideas

- **Dropping `Major_*`/`Minor_*` legacy columns:** explicitly deferred to the next milestone post-v3.0, per ROADMAP.md ("Drop the legacy `Major_*`/`Minor_*` column aliases after the one-release deprecation window (COMPAT-03)"). Phase 9 only validates their presence.
- **Renaming `major_mapping`/`minor_mapping` channel emit identifiers to role-semantic names:** if the N-slot generalization makes these obsolete, the cleanup can be its own tidy-up commit rather than a Phase 9 hard requirement.
- **Cross-platform Singularity/Conda nf-test profiles:** not in Phase 9 scope; existing CI runs Docker only.

</deferred>

---

*Phase: 09-compatibility-filename-migration-regression-suite*
*Context gathered: 2026-06-14*
