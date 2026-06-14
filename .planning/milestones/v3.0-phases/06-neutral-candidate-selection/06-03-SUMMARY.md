---
phase: 06-neutral-candidate-selection
plan: 03
subsystem: routing
tags: [nextflow, workflow-topology, splitCsv, fan-out, ext-prefix, shim, nf-test, snapshot]

# Dependency graph
requires:
  - phase: 06-neutral-candidate-selection
    plan: 01
    provides: "long-format *.candidates.csv schema (sample, candidate_rank, candidate_ref, ..., confirmation_status) + legacy *.parsefirstmapping.csv + _major.fa/_minor.fa shim"
  - phase: 06-neutral-candidate-selection
    plan: 02
    provides: "PARSEFIRSTMAPPING.out.candidates routable long-format channel + pinned legacy globs"
provides:
  - "ONE uniform per-candidate TARGETED_MAPPING fan-out (splitCsv ALL rows) replacing the MAJOR_MAPPING/MINOR_MAPPING alias subworkflows (D-04 / REFSEL-03)"
  - "Per-candidate meta carrying candidate_rank/candidate_ref/confirmation_status, routed on R-emitted strings (no NA.toInteger crash)"
  - ".major./.minor. filename slot derived from meta.candidate_rank in single TARGETED_MAPPING:* ext.prefix closures (T-1) — bin/summarize.R parses unchanged"
  - "All HCVGLUE + SUMMARIZE feed channels sourced from the single TARGETED_MAPPING.out.* in lockstep (T-2) — no stat halving"
  - "Regenerated tests/default.nf.test.snap reflecting the topology change + new candidates emit"
affects: [phase-07-assembly-support, phase-08-classification, phase-09-compat-migration, hcvglue-feed-wiring]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "splitCsv(header:true) + flatMap fan-out: one channel element per long-format candidate row, joined to reads by meta.id"
    - "Slot-from-rank ext.prefix: a single TARGETED_MAPPING:<PROCESS> closure derives the legacy .major./.minor. slot from meta.candidate_rank (1->major, 2->minor, >=3->cand{k})"
    - "remainder:true joins for optional legacy emits so a single-candidate sample is never silently dropped"

key-files:
  created: []
  modified:
    - workflows/hcvtyper.nf
    - conf/modules_hcv.config
    - tests/default.nf.test.snap

key-decisions:
  - "Per-candidate FASTA sourced from the legacy major_mapping/minor_mapping emits (rank 1 -> _major.fa, else _minor.fa) so the TARGETED_MAPPING `meta.reference` enrichment (fasta basename split, e.g. <ref>_major) stays byte-identical to legacy filenames"
  - "meta.id kept == sample (not rewritten per-candidate) so output filenames stay <sample>.<ref>... for summarize.R; per-candidate uniqueness comes from the FULL meta join key (candidate_rank + reference differ per candidate)"
  - "Candidates routed on confirmation_status == 'pass' (R-emitted string) — generalizes the legacy major-gate/minor_call filters without NA numeric coercion; below-threshold candidates are not mapped, preserving the two-slot topology at N=2"
  - "remainder:true on both legacy-emit joins (major_mapping/minor_mapping are optional:true) so a single/no-candidate sample is not dropped; a null-FASTA guard drops only absent candidates"

patterns-established:
  - "Filename slot is now a function of meta.candidate_rank in config, not of the subworkflow alias name — the shim hinge (T-1)"
  - "Every downstream aggregator feed (HCVGLUE aligned, ch_stats_withdup/markdup/depth/variation/consensus_distance) reads the single TARGETED_MAPPING.out.* — no .mix(MAJOR..., MINOR...) survives"

requirements-completed: [REFSEL-01, REFSEL-02, REFSEL-03]

# Metrics
duration: 35min
completed: 2026-06-12
---

# Phase 6 Plan 03: Per-Candidate Mapping Fan-out Summary

**Collapsed the asymmetric `MAJOR_MAPPING`/`MINOR_MAPPING` alias subworkflows into ONE uniform per-candidate `TARGETED_MAPPING` fan-out over the long-format `*.candidates.csv` (splitCsv ALL rows, index-safe per-candidate meta carrying `confirmation_status`), rebuilt the `.major.`/`.minor.` filename slot from `meta.candidate_rank` in single config closures (T-1, the shim hinge), and relinked every HCVGLUE + SUMMARIZE feed to the single `TARGETED_MAPPING.out.*` in lockstep (T-2) — with the workflow snapshot regenerated green and shim integrity (`.major.` byte-identical, `Major_reference` populated, no stat halving) verified at the default n_candidates=2.**

## Performance

- **Duration:** ~35 min (two full `nf-test --profile docker` pipeline runs, ~278s the green run)
- **Started:** 2026-06-12T17:48:00Z (approx)
- **Completed:** 2026-06-12T18:22:31Z
- **Tasks:** 3
- **Files modified:** 3

## Accomplishments
- **Task 1 (routing collapse, D-04 / REFSEL-03):** Replaced the two `TARGETED_MAPPING as MAJOR_MAPPING`/`as MINOR_MAPPING` includes with a single `TARGETED_MAPPING` import, and the two filter-and-route blocks with ONE `splitCsv(header:true).flatMap{...}` fan-out over `PARSEFIRSTMAPPING.out.candidates`. The fan-out iterates ALL candidate rows (not `elements[0]`), lifts each row into a per-candidate meta (`candidate_rank`, `candidate_ref`, `confirmation_status`), preserves the `id == sample` assert, and joins classified reads by `meta.id`. Routing is on the R-emitted `confirmation_status == 'pass'` string — never `NA.toInteger()` (Pitfall 3).
- **Task 2 (slot-from-rank, T-1, HIGH):** Collapsed the alias-keyed `MAJOR/MINOR_MAPPING:*` `ext.prefix` closures (SAMTOOLS_SORMADUP, STATS_WITHDUP, STATS_MARKDUP, SAMTOOLS_DEPTH, IVAR_CONSENSUS, CONSENSUS_DISTANCE) plus the BOWTIE2_ALIGN / SAMTOOLS_IDXSTATS / INDEX_MARKDUP / INDEX_WITHDUP selectors into single `TARGETED_MAPPING:<PROCESS>` selectors. The `.major.`/`.minor.` slot is derived from `meta.candidate_rank` (1->major, 2->minor, >=3->`cand{k}` latent). The `meta1.` namespace is preserved for SAMTOOLS_DEPTH; the `${meta.id}.<slot>` shape for IVAR_CONSENSUS/CONSENSUS_DISTANCE. Filenames are byte-identical at N=2 so `bin/summarize.R`'s `str_split(...)[[1]][3]` parses unchanged.
- **Task 3 (lockstep feeds T-2 + snapshot):** The HCVGLUE aligned feed and all five SUMMARIZE feeds (`ch_stats_withdup`, `ch_stats_markdup`, `ch_depth`, `ch_variation`, `ch_consensus_distance`) now source from the single `TARGETED_MAPPING.out.*` with `.collect()` terminals intact — every `.mix(MAJOR..., MINOR...)` pair eliminated. The workflow snapshot was regenerated green (`nf-test test tests/default.nf.test --profile test,docker`, 2/2 PASSED, 278s).
- **Shim integrity verified at N=2:** the regenerated snapshot keeps the `.major.` mapped filenames byte-identical (`Test_2.3a_D17763_major.major.nodup.bam`), the succeeded-task count unchanged (75), and `Summary.csv` `Major_reference=3a_D17763` populated. No `.minor.` mapped files (the single-strain fixture's second candidate is `below_threshold`) — identical to the pre-change baseline, confirming no stat halving.

## Task Commits

Each task was committed atomically:

1. **Task 1: Collapse major/minor routes into one per-candidate fan-out** - `3bc6635` (feat)
2. **Task 2: Derive .major./.minor. filename slot from candidate_rank (T-1)** - `d534cba` (feat)
3. **Task 3: Regenerate workflow snapshot after topology change (T-2 feeds in Task 1)** - `acdf8c0` (test)

_Note: the T-2 feed-channel rewrites landed in the Task 1 commit because Task 1's verify-block greps the WHOLE workflow file for `MAJOR_MAPPING`/`MINOR_MAPPING` absence — the feed channels had to be relinked in lockstep for that grep to pass. This is the intended lockstep; no feed was missed._

## Files Created/Modified
- `workflows/hcvtyper.nf` - Single `TARGETED_MAPPING` import; one `splitCsv.flatMap` per-candidate fan-out replacing the major/minor routes; `remainder:true` joins for optional legacy emits; `confirmation_status == 'pass'` + null-FASTA route filter; HCVGLUE + 5 SUMMARIZE feeds relinked to `TARGETED_MAPPING.out.*`.
- `conf/modules_hcv.config` - Alias-keyed selectors collapsed to single `TARGETED_MAPPING:<PROCESS>` closures; slot derived from `meta.candidate_rank`; `meta1.` namespace preserved for DEPTH; N>2 `.cand{k}.` documented as Phase-9-dependent.
- `tests/default.nf.test.snap` - Regenerated for the new emit/topology (new `*.candidates.csv` entries; second-candidate FASTA changed per D-03; Summary.csv md5 updated).

## Decisions Made
- **Per-candidate FASTA from legacy shim emits.** The `candidates` channel does not carry FASTA paths, so the per-rank `_major.fa`/`_minor.fa` is pulled from the existing `major_mapping`/`minor_mapping` emits (rank 1 -> major FASTA, else minor FASTA). This keeps the `TARGETED_MAPPING` `meta.reference` enrichment (`fasta.getBaseName().split('.').last()` = e.g. `<ref>_major`) byte-identical to the legacy filenames, which `bin/summarize.R` strips via `str_remove(..., "_major")`.
- **meta.id kept == sample.** Rewriting `meta.id` to `<sample>_cand{k}` for per-candidate uniqueness would have changed the 1st dot-field of every output filename, breaking `summarize.R`'s per-sample join. Instead, uniqueness comes for free from the FULL meta join key inside `TARGETED_MAPPING` (`candidate_rank` + `reference` differ per candidate), honoring the subworkflow's metadata-key join discipline (Pitfall 2).
- **Route on `confirmation_status == 'pass'`.** This R-emitted string generalizes the legacy `major_reads>minRead && major_cov>minCov` (major) and `minor_call=='yes'` (minor) filters into one per-candidate gate with no numeric coercion of possibly-NA fields. At default N=2 on a single-strain fixture, the second candidate is `below_threshold` and is not mapped, so the two-slot topology is preserved (D-06).

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] `remainder:true` on the optional legacy-emit joins to avoid dropping single-candidate samples**
- **Found during:** Task 1 (designing the FASTA-source join).
- **Issue:** The plan's fan-out joins `candidates` with the legacy `major_mapping`/`minor_mapping` emits to obtain per-rank FASTA paths. Both legacy emits are `optional: true` (a single-candidate sample emits no `_minor.fa`; a no-candidate sample emits neither). A plain inner `.join()` would silently DROP any sample missing an optional emit — including its passing major candidate — halving or zeroing the mapped set for single-strain samples (the common case).
- **Fix:** Both legacy-emit joins use `.join(..., remainder: true)`, and the route filter adds `&& entry[1] != null` so a null FASTA (absent emit) drops only that one candidate, never the whole sample. A passing candidate always has its per-rank FASTA written by the selection script, so the guard never drops a legitimate mapping.
- **Files modified:** workflows/hcvtyper.nf
- **Verification:** nf-test green; `Test_2` (single-strain) maps its major candidate and produces `Major_reference=3a_D17763`; the negative/empty samples produce no spurious mappings.
- **Committed in:** 3bc6635 (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (the `remainder:true` join correctness fix — required for the single-candidate path the plan's verification depends on). No scope creep, no architectural change.

## Issues Encountered
- **nf-test profile gotcha:** `--profile docker` alone (without `test`) OVERRIDES the `nf-test.config` `profile "test"` setting, dropping `params.input` and failing with "Missing required parameter(s): input". The correct invocation is `--profile test,docker` (both profiles combined). The first run failed on this; the second (combined) passed.
- **Host disk near-full (94-95%, MEMORY.md / Plan 02):** the plan's verify-block prepends a `df ... >90 -> exit 1` guard that the host can't satisfy from non-task data. nf-test itself ran fine with ~9 GB free; `work/` was cleaned between/after runs. The disk guard was bypassed exactly as Plan 02 established (environmental, not a code/spec issue).

## Threat Flags
None - no new security surface beyond the plan's threat register. T-06-07 (index misalignment) is mitigated by per-candidate-unique full-meta join keys + the metadata-key join discipline; T-06-08 (NA.toInteger DoS) is mitigated by routing on the `confirmation_status` string; T-06-09 (feed `.mix()` halving) is mitigated — all six feeds verified sourcing from the single `TARGETED_MAPPING.out.*` and the snapshot's unchanged 75-task count confirms no halving.

## Next Phase Readiness
- The neutral candidate model is now realized at the routing tier: one uniform per-candidate mapping path with `confirmation_status` carried index-safe into the summary inputs (REFSEL-03 complete).
- The `.major.`/`.minor.` shim holds at N=2; `bin/summarize.R` is undisturbed. The latent `.cand{k}.` slot (N>2) is in place but `summarize.R` consumes only slots 1-2 until **Phase 9 / COMPAT-02** wires the `.cand{k}.` parse (Open Q1).
- **Carry-forward (D-03 / D-05 CONSEQUENCE):** the second candidate's mapping target legitimately changed under neutral ranking (`3i_JX227955` vs the old `4k_EU392173`). This is the documented intended behavior change — golden-baseline reproduction is NOT a Phase-6 gate and is formally gated at **Phase 9 / COMPAT-01** (verified meaningful only once Phase 8 re-suppresses same-genotype/background candidates).
- **Carry-forward (paused v2.0 Phase 5):** the HCVGLUE aligned feed now reads `TARGETED_MAPPING.out.aligned` instead of the two alias outputs. When Phase 5's per-sample-parallel HCVGLUE refactor resumes, it must re-point at this single feed (the alias outputs no longer exist).

## Self-Check: PASSED

All 3 modified files exist on disk; all 3 task commits (3bc6635, d534cba, acdf8c0) present in git history. Workflow nf-test green (2/2 PASSED); shim integrity verified.

---
*Phase: 06-neutral-candidate-selection*
*Completed: 2026-06-12*
