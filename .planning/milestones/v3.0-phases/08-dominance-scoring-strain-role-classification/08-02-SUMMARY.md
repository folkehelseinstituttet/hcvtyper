---
phase: 08-dominance-scoring-strain-role-classification
plan: 02
subsystem: pipeline-summarize
tags: [R, tidyverse, nextflow, dominance-score, strain-role, classify_roles, summarize, review_flag, candidates-csv]

# Dependency graph
requires:
  - phase: 08-01
    provides: "pure classify_roles.R helper (score_candidates + classify_roles + is_valid_minor), reconciled denovo floors 1000/2.0/90, score_weight_* params, cv_evenness + score_weight_* wiring (Task 1, beb9e37)"
  - phase: 07-per-genotype-assembly-support
    provides: "join_assembly_support() genotype-level join feeding candidate_support"
  - phase: 06-neutral-candidate-selection
    provides: "long-format *.candidates.csv (the classifier input)"
provides:
  - "summarize.R runs score_candidates() + classify_roles() over the Phase-7-joined candidate frame; each candidate gets role / dominance_score / role_reason and one overall_sample_call per sample"
  - "legacy apply_denovo_layer / minor_denovo_status / coinfection_flag consumption retired (D-15) — one confirmation system"
  - "Summary.csv carries overall_sample_call + Major_/Minor_role_* slots, no longer carries minor_denovo_status/coinfection_flag"
  - "enriched long candidates.csv written with every candidate incl. background + role/dominance_score/role_reason (CLASS-03)"
  - "review_flag rewired onto overall_sample_call / role_reason (any_refuted_denovo, any_uncorroborated)"
  - "classify_roles.R staged as a SUMMARIZE path() input + passed at the workflow call site; stub schema updated"
affects: [phase-09-compat-filename-migration, phase-09-Major-Minor-column-aliasing, phase-09-regression-suite, nf-test-snapshot-regen]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "R-emits-decision: the N-candidate role model is computed entirely in summarize.R over the long candidate frame; Nextflow only stages the helper"
    - "Per-sample role rollup (group_by/summarise) drives the rewired review_flag sentences instead of the retired per-row minor_denovo_status"
    - "Wide role slots (Major_/Minor_role_*) added ALONGSIDE legacy Major_*/Minor_* mapping slots — column aliasing deferred to Phase 9"

key-files:
  created: []
  modified:
    - bin/summarize.R
    - modules/local/summarize/main.nf
    - workflows/hcvtyper.nf
    - bin/tests/test_summarize_denovo.R

key-decisions:
  - "classify_roles() floor wired to min_targeted_read/min_targeted_cov (args[9]/[10] = configured minRead/minCov) — the b23e021 partial wiring called it with undefined minRead/minCov vars (Rule 1 bug fix)"
  - "review_flag rewired onto a per-sample role rollup (any_refuted_denovo / any_uncorroborated) + overall_sample_call + gate_flag + subtype-match, since refuted/uncorroborated reasons live on background candidates not surfaced in the wide slots"
  - "Wide role slots named Major_role_*/Minor_role_* (additive) rather than overwriting Major_*/Minor_* — the legacy column ALIASING is Phase 9 (COMPAT-03); Phase 8 only retires the legacy LOGIC (D-15)"

patterns-established:
  - "Pattern: enriched long candidates.csv is the CLASS-03 surface — write_csv(candidate_support) AFTER classification so backgrounds are never dropped"
  - "Pattern: D-15 retirement asserted in test_summarize_denovo.R Block-1 (no uncommented apply_denovo_layer call; sources+calls classify_roles) while denovo_layer.R remains sourced + unit-tested directly"

requirements-completed: [SCORE-01, SCORE-02, CLASS-01, CLASS-02, CLASS-03, CLASS-04]

# Metrics
duration: ~35min
completed: 2026-06-14
---

# Phase 8 Plan 02: Wire Dominance Score + Strain-Role Classifier into SUMMARIZE Summary

**summarize.R now runs the Phase-8 N-candidate role classifier over the Phase-7-joined candidate frame — emitting per-candidate role/dominance_score/role_reason + one overall_sample_call per sample, retiring the legacy apply_denovo_layer/minor_denovo_status/coinfection_flag path (D-15), rewiring review_flag onto roles, writing the enriched long candidates.csv, and staging classify_roles.R into the SUMMARIZE module.**

## Performance

- **Duration:** ~35 min (this session; Task 1 + the b23e021 partial Task 2 landed in a prior session)
- **Started:** 2026-06-14T08:48Z (this session)
- **Completed:** 2026-06-14T09:06Z
- **Tasks:** 3 (Task 1 pre-committed beb9e37; Task 2 completed; Task 3 executed)
- **Files modified:** 4

## Accomplishments
- Completed Task 2: classifier runs in summarize.R; legacy confirmation layer retired (D-15); roles mapped to wide Major_/Minor_role_* slots; overall_sample_call added; enriched long candidates.csv written; review_flag rewired onto roles.
- Fixed a latent Rule-1 bug left by the b23e021 partial wiring: classify_roles() was called with undefined `minRead`/`minCov` vars — rewired to `min_targeted_read`/`min_targeted_cov` (the configured minRead/minCov at args[9]/[10]).
- Executed Task 3: classify_roles.R staged as a path() input on the SUMMARIZE module and passed at the hcvtyper.nf call site; stub Summary.csv + summary_mqc.tsv schemas updated (minor_denovo_status out, overall_sample_call + 8 role columns in) with all four stub rows realigned to 92 fields.
- Full R unit suite green in the tidyverse container; `nextflow config` parses the edited module + workflow cleanly; a standalone runtime smoke test confirmed the false-4g→background/refuted_denovo and genuine→co-infection chain end-to-end.

## Task Commits

1. **Task 1: cv_evenness in cov loop + score_weight_* ext.args** - `beb9e37` (feat) — *committed in a prior session*
2. *(partial Task 2 prior wiring)* - `b23e021` (chore "updates") — source(classify_roles.R) + cv_by_ref join + score/classify calls
3. **Task 2: run classifier, retire legacy layer, rewire review_flag, write candidates.csv** - `2e9f095` (feat)
4. **Task 3: stage classify_roles.R + update stub schema** - `483a701` (feat)

**Plan metadata:** _(final docs commit — this SUMMARY + STATE + ROADMAP)_

## Files Created/Modified
- `bin/summarize.R` - classify_roles() floor fix (min_targeted_read/cov); enriched candidates.csv write; wide role-slot mapping (Major_/Minor_role_* + overall_sample_call); legacy apply_denovo_layer/coinfection_flag removed; review_flag rewired onto per-sample role rollup; reorder select() drops minor_denovo_status/coinfection_flag, adds overall_sample_call + role slots.
- `modules/local/summarize/main.nf` - `path(classify_roles)` input; stub Summary.csv heredoc + summary_mqc.tsv printf headers/data realigned to the new 92-field schema.
- `workflows/hcvtyper.nf` - `file("${projectDir}/bin/classify_roles.R")` passed at the SUMMARIZE call site.
- `bin/tests/test_summarize_denovo.R` - Block-1 sanity rewritten to assert the D-15 retirement contract (no uncommented apply_denovo_layer call; sources+calls classify_roles) instead of the old "summarize.R calls apply_denovo_layer()" assertion.

## Decisions Made
- **Floor source:** classify_roles() uses `min_targeted_read`/`min_targeted_cov` (args[9]/[10], which the SUMMARIZE ext.args populates from `${params.minRead}`/`${params.minCov}`). This is exactly the D-07/D-09 shared threshold set.
- **review_flag triggers:** derived from a per-sample rollup of the long classified frame (`any_refuted_denovo`, `any_uncorroborated`) plus `overall_sample_call`, `gate_flag`, and the existing subtype-match columns — because the `refuted_denovo`/`uncorroborated_kept` reasons attach to background candidates that the wide Major_/Minor_ slots do not surface.
- **Additive role columns:** introduced `Major_role_*`/`Minor_role_*` rather than overwriting the legacy `Major_*`/`Minor_*` mapping/coverage slots; the legacy column aliasing onto the role model is explicitly Phase 9 (COMPAT-03). Phase 8 retires only the legacy confirmation LOGIC (D-15).

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] classify_roles() called with undefined minRead/minCov**
- **Found during:** Task 2 (completing the b23e021 partial wiring)
- **Issue:** The prior partial commit `b23e021` called `classify_roles(candidate_support, minRead = minRead, minCov = minCov, ...)`, but `minRead`/`minCov` are not defined anywhere in summarize.R — the configured values arrive as `min_targeted_read`/`min_targeted_cov` (args[9]/[10], populated from `${params.minRead}`/`${params.minCov}`). At runtime this would error with "object 'minRead' not found".
- **Fix:** Rewired the call to `minRead = min_targeted_read, minCov = min_targeted_cov` with an explanatory comment tying it to the ext.args positional mapping and the D-07/D-09 shared-threshold decision.
- **Files modified:** bin/summarize.R
- **Verification:** `Rscript -e 'parse(...)'` PARSE_OK; standalone runtime smoke test exercised the floor with concrete reads/cov and produced the expected dominant/background/co-infection roles.
- **Committed in:** `2e9f095` (Task 2 commit)

**2. [Rule 3 - Blocking] test_summarize_denovo.R Block-1 asserted the retired contract**
- **Found during:** Task 2 (running the suite after removing the apply_denovo_layer call)
- **Issue:** `bin/tests/test_summarize_denovo.R:74` hard-asserted `summarize.R` must contain an `apply_denovo_layer(` call. D-15 retires that consumption, so the assertion would fail the suite even though the retirement is the intended behaviour. (Blocks 2-4 test `apply_denovo_layer()` as a sourced function directly — those stay valid because denovo_layer.R is still sourced.)
- **Fix:** Rewrote only the Block-1 sanity check to assert the D-15 contract: summarize.R still sources denovo_layer.R, must NOT have an uncommented apply_denovo_layer( call, and must source + call classify_roles(). Left the apply_denovo_layer() fixture tests (Blocks 2-4) untouched.
- **Files modified:** bin/tests/test_summarize_denovo.R
- **Verification:** `bash bin/tests/run_all.sh` in the tidyverse container — ALL R TESTS PASSED.
- **Committed in:** `2e9f095` (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (1 bug, 1 blocking)
**Impact on plan:** Both auto-fixes were necessary for correctness and to land the D-15 retirement the plan mandates. No scope creep — the test edit is the minimal change needed to encode the new contract; the legacy fixture coverage is preserved.

## Issues Encountered
- **Host R lacks tidyverse:** the system R (4.3.3) has no tidyverse and the stale R-4.1 user library fails to load under 4.3.3 (ICU `libicui18n.so.70` mismatch). Resolved by running the R suite inside the project's `jonbra/tidyverse_seqinr:2.0` docker image (`docker run --rm -v "$PWD":/work -w /work ... bash bin/tests/run_all.sh`), which matches the SUMMARIZE container's tidyverse.
- **Stub data-row field alignment:** an early `awk`-based rewrite of the stub `summary_mqc.tsv` data row expanded `\t` escapes into literal tabs and broke the printf line; recovered by rebuilding the data row programmatically from the verified 92-field Summary.csv header so all four stub rows (csv header/data, mqc header/data) are exactly 92 aligned fields. Verified with a Python field-count check.

## Known Stubs
None — the classifier is fully wired and consumed; `overall_sample_call` and the long `candidates.csv` are populated from real candidate data, not placeholders. The stub block in `modules/local/summarize/main.nf` is the Nextflow `stub:` directive (a fixed placeholder for `-stub-run`), not a data stub in the runtime path.

## Deferred Items (tracked, not blockers)
- **nf-test SUMMARIZE snapshot regen** for the new Summary.csv / summary_mqc.tsv schema (overall_sample_call + role columns in, minor_denovo_status out) — folded into the already-pending Phase-7 SUMMARIZE snapshot regen (docker/disk-constrained; run via `PATH=~/.nf-test:$PATH nf-test ... --profile docker`). Per CONTEXT scope this is NOT a Phase-8 gate.
- **Major_*/Minor_* legacy column ALIASING onto the role model** and the **.major./.minor. → .cand1./.cand2. filename migration** — Phase 9 (COMPAT-02/03).

## Next Phase Readiness
- The N-candidate role model is live end-to-end in SUMMARIZE: Phase-6 candidates → Phase-7 assembly-support join → cv_evenness → score_candidates() → classify_roles() → wide Summary.csv (overall_sample_call + role slots) + enriched long candidates.csv.
- Phase 9 can now layer the `Major_*`/`Minor_*` column aliasing + filename-slot migration + the golden-baseline regression suite on top of the retired-legacy / role-driven Summary.
- Outstanding for Phase 9: regenerate the deferred SUMMARIZE nf-test snapshot against the new schema.

## Self-Check: PASSED

- FOUND: bin/summarize.R, modules/local/summarize/main.nf, workflows/hcvtyper.nf, bin/tests/test_summarize_denovo.R
- FOUND commits: 2e9f095 (Task 2), 483a701 (Task 3)
- All R unit tests green (jonbra/tidyverse_seqinr:2.0); `nextflow config` parses; runtime smoke test of the classify→wide-mapping→review_flag chain succeeded.

---
*Phase: 08-dominance-scoring-strain-role-classification*
*Completed: 2026-06-14*
