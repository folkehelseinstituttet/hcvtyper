# HCVTyper — Living Retrospective

## Milestone: v1.0 — De novo Minor Confirmation

**Shipped:** 2026-06-08
**Phases:** 4 | **Plans:** 12 | **Tasks:** 19
**Git range:** `refactor(01-01)` → `feat(04-03)` (2026-06-06 → 2026-06-07, ~2 days)

### What Was Built

De novo/BLAST evidence and a major-gate now drive minor-strain reporting end-to-end, behind `--denovo_confirm_minor` (default ON, flag-OFF reproduces legacy). Major-gate (Change 1) blocks minors behind a failed major; the confirm/refute/fall-back layer (Change 2) cross-checks the candidate minor against calibrated substantial-contig evidence at genotype level; `minor_denovo_status` surfaces the basis for each call; a 4-file R regression suite (`bin/tests/run_all.sh`) wired into CI locks the correctness contract.

### What Worked

- **Risk-phased decomposition.** Hardening + dead-code removal (Phase 1) → zero-behaviour-change plumbing (Phase 2) → decision logic (Phase 3) → regression guard (Phase 4). Each phase was independently verifiable; the additive-only plumbing in Phase 2 de-risked the load-bearing change in Phase 3.
- **R-emits-decision, Nextflow-routes.** Moving the gate decision into the R script as `minor_call`/`gate_flag` columns killed the `NA.toInteger()` Groovy crash class instead of patching around it.
- **Pure sourced helpers.** `classify_minor_denovo()` / `apply_denovo_layer()` as side-effect-free sourced functions let the tests call the *real* code (no drift), which made retroactive verification (UAT via the live suite) genuinely meaningful.
- **Evidence-table-derived fixtures.** Tests assert real outcomes for named samples (ERR1810453, sim1_1a_single, etc.), so the suite doubles as executable documentation of the benchmarked failure modes.

### What Was Inefficient

- **Verification ran out of order.** Phases 03/04 were left unverified (`human_needed` / no VERIFICATION.md) at milestone-close time; UAT + security + audit all had to be done retroactively in one session before the close could proceed. Running verify/secure as part of each execute-phase would have avoided the scramble.
- **Skip-assembly path never exercised until Phase 2 plan 3.** Two latent plumbing bugs (undefined `BLASTPARSE.out`; vanishing denovo columns) only surfaced when the `--skip_assembly` nf-test case was finally added.
- **Dangling fixture symlinks.** Phase 3's nf-test depended on symlinks to a deleted external work dir; had to be restored from in-repo `minimal_test/`.

### Patterns Established

- Gate/decision logic lives in the R layer and is surfaced as explicit status columns; Nextflow filters on those strings.
- New params are declared in all four locations (nextflow.config, schema, modules_hcv.config, R arg passthrough) — the `minDenovoLength` latent-gap anti-pattern is the cautionary tale.
- Flag-gated behaviour ships with a committed golden-baseline differential proving the OFF path reproduces legacy output.
- Threat models authored at plan time → `/gsd-secure-phase` verifies mitigations against code (State A/B), short-circuiting when `threats_open: 0` and the register was plan-authored.

### Key Lessons

- Verify and secure each phase *during* execution, not at milestone close — retroactive verification works but compresses risk into one session.
- For an analysis pipeline with no network/auth surface, the real security surface is the CI job (image pinning, no PR-input interpolation) — worth an explicit threat in the register.
- A passing R unit/regression suite that calls the real functions is acceptable verification evidence when a full docker pipeline run is impractical (no FASTAs on disk, host at 95%).

### Cost Observations

- Sessions: milestone close + retroactive verification done in 1 session.
- Model mix: Opus (orchestration, audits, security verification, doc evolution).
- Notable: the Phase-4 R suite (~30s, no docker) substituted for a multi-GB docker pipeline run, keeping the close feasible on a near-full disk.

---


## Milestone: v3.0 — Strain Model Redraw

**Shipped:** 2026-06-14
**Phases:** 4 (6-9) | **Plans:** 11 | **Tasks:** ~31

### What Was Built
- Neutral read-recruitment candidate ranking (`cand_1..cand_n`); `is_valid_minor()` validity filtering deleted from selection
- Per-genotype assembly support roll-up (`*.assembly_support.csv`) + genotype-level join helper (`assembly_support_join.R`)
- Breadth-evenness-weighted dominance score (evenness 3× reads) + strain-role classification (dominant/co-infection/background) with `overall_sample_call`
- `.cand{rank}.` filename-slot migration (6 closures + FASTA loop + `summarize.R` join parse), lockstepped in one phase
- Legacy `Major_*`/`Minor_*` column aliases + regression suite (`test_compat.R`) covering all 4 COMPAT requirements

### What Worked
- Staging pure sourced R helpers (`assembly_support_join.R`, `classify_roles.R`) as `path()` inputs to SUMMARIZE — clean separation, easy to unit-test
- TDD with subprocess-contract tests (system2 on inline fixtures) — caught real bugs before integration
- Locking the highest-risk change (filename-slot migration + parse cutover) as the final phase with integration test coverage

### What Was Inefficient
- Two post-verification fix commits needed for Phase 9 (CR-01/CR-02 missed during execution) — code review after each plan would have caught these earlier
- Phase 9 VERIFICATION.md not updated after gap closure; leaves stale status artifact

### Patterns Established
- Score-then-classify pipeline: `score_candidates()` → `classify_roles()` over the joined candidate frame (Phase 7 output)
- Integration-gate pattern: every helper staged into the SUMMARIZE module as a `path()` input, verified at module arity
- Calibration-as-test: evidence-table fixtures drive the score-weight defaults in `test_dominance_score.R`

### Key Lessons
- The highest-risk integration (consumer/producer rename lockstep) is best scheduled as the final phase and verified by a subprocess-contract regression test, not just nf-test snapshots
- Static analysis can verify 13/13 truths for R logic that requires Docker to run — a good stopping point when the host lacks tidyverse
- Deferred snapshot regen (ENOSPC) is a recurring bottleneck for this pipeline; consider a CI job or dedicated disk quota

### Cost Observations
- 3 days, 45 commits, 30 files, +3152/-529 LOC
- All phases executed inline (no worktrees); gsd-verifier and integration-checker ran as subagents

## Cross-Milestone Trends

| Metric | v1.0 |
|--------|------|
| Phases | 4 |
| Plans | 12 |
| Tasks | 19 |
| Audit status | tech_debt (20/20 reqs) |
| Duration | ~2 days |

_First milestone — trends accumulate from v2 onward._
