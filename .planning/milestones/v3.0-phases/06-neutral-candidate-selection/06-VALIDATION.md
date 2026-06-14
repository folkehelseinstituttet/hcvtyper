---
phase: 6
slug: neutral-candidate-selection
status: ready
nyquist_compliant: true
wave_0_complete: false
created: 2026-06-12
---

# Phase 6 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.
> Mirrored from 06-RESEARCH.md §"Validation Architecture".

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | nf-test 0.9.3 (module/workflow snapshots) + Rscript subprocess-contract tests |
| **Config file** | `nf-test.config` (repo root); R harness `bin/tests/run_all.sh` |
| **Quick run command** | `bash bin/tests/run_all.sh` (R logic, seconds) |
| **Module run command** | `PATH=~/.nf-test:$PATH nf-test test modules/local/parsefirstmapping/tests/main.nf.test --profile docker` |
| **Full suite command** | `PATH=~/.nf-test:$PATH nf-test test --profile docker` |
| **Estimated runtime** | R tests ~seconds; per-module nf-test ~1–3 min; full workflow nf-test ~5–15 min (docker) |

> Host disk is near-full (~93%+, MEMORY.md). Prefer the fast R subprocess tests during iteration; gate full nf-test runs behind a disk-headroom check and `sudo rm -rf work/` + `docker volume prune -f` on ENOSPC.

---

## Sampling Rate

- **After every task commit:** `bash bin/tests/run_all.sh` (+ the single `parsefirstmapping` nf-test when the module/script changed)
- **After every plan wave:** full module + `tests/default.nf.test` workflow snapshot under `--profile docker`
- **Before `/gsd-verify-work`:** full suite green
- **Phase gate:** full suite green — **explicitly NOT** golden-baseline reproduction (D-05; that gate is Phase 9 / COMPAT-01)
- **Max feedback latency:** R tests < ~30s; full workflow gate minutes (run at wave boundaries only)

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|----------|-----------|-------------------|-------------|--------|
| 06-01-01 | 01 | 1 | REFSEL-01/02/03 | Neutral top-N-by-reads ranking, distinct-subtype dedup, long-format CSV + legacy shim columns, `n_candidates` param | source + R unit | `bash bin/tests/run_all.sh` | ⚠️ script exists, rewrite | ⬜ pending |
| 06-01-02 | 01 | 1 | REFSEL-01 | R subprocess contract: ranking + dedup + shim reconstruction | unit (R) | `Rscript bin/tests/test_candidate_selection.R` | ❌ W0 | ⬜ pending |
| 06-01-03 | 01 | 1 | REFSEL-02 | `params.n_candidates` declared (schema + config), default 2, integer ≥ 1 | schema/source | `nextflow run . -profile test,docker --validate_params` (smoke) + grep schema | ⚠️ param not yet declared | ⬜ pending |
| 06-02-01 | 02 | 2 | REFSEL-01/03 | PARSEFIRSTMAPPING surfaces `candidates` emit; passes `n_candidates`; legacy emits intact | module (source) | grep `emit: candidates` / `emit: major_mapping` / `emit: minor_mapping` in main.nf | ⚠️ exists, extend | ⬜ pending |
| 06-02-02 | 02 | 2 | REFSEL-01 | nf-test asserts candidates emit + long-format header; snapshot regen (T-3 churn expected) | module (nf-test) | `nf-test test modules/local/parsefirstmapping/tests/main.nf.test --profile docker` | ⚠️ exists, new assertions + snap | ⬜ pending |
| 06-03-01 | 03 | 3 | REFSEL-03 | Single uniform per-candidate fan-out replacing MAJOR/MINOR_MAPPING; `confirmation_status` carried via meta-key joins | workflow (source) | grep `! MAJOR_MAPPING\|MINOR_MAPPING` + `TARGETED_MAPPING.out` count in hcvtyper.nf | ⚠️ exists, rewrite | ⬜ pending |
| 06-03-02 | 03 | 3 | REFSEL-03 / shim (T-1) | `.major.`/`.minor.` slot derived from `candidate_rank` in `ext.prefix` closures | config (source) | grep `candidate_rank` in `conf/modules_hcv.config` ext.prefix | ⚠️ alias-keyed today | ⬜ pending |
| 06-03-03 | 03 | 3 | REFSEL-03 / shim (T-2) | Lockstep feed-channel `.mix()` collapse for HCVGLUE + SUMMARIZE; workflow snapshot regen | workflow (nf-test) | `nf-test test tests/default.nf.test --profile docker` | ⚠️ exists, snap regen | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `bin/tests/test_candidate_selection.R` — covers REFSEL-01 (ranking + distinct-subtype dedup + shim column reconstruction). Auto-globbed by `run_all.sh`. *(written in Plan 06-01 Task 2)*
- [ ] New nf-test assertions in `modules/local/parsefirstmapping/tests/main.nf.test` for the `candidates` emit + long-format CSV; regenerate `.snap`. *(Plan 06-02 Task 2)*
- [ ] Update `parsefirstmapping` GATE-03/04 (2k1b / no-minor) assertions to reflect no-validity-filtering selection (D-05). *(Plan 06-02 Task 2)*
- [ ] Regenerate `tests/default.nf.test.snap` after the topology change (expected snapshot churn). *(Plan 06-03 Task 3)*
- [ ] Optional fixture: a 2-distinct-subtype sample where the 2nd candidate differs from today's minor, to lock the NEW D-03/D-05 behavior.

> `wave_0_complete: false` until these test artifacts exist and pass; they are authored inside the plan tasks above (not a separate Wave 0 plan).

---

## What Phase 6 verification SHOULD assert (per D-05, instead of ROADMAP criterion #5)

1. At default `n_candidates=2`, the run completes without crash on the regression fixtures (topology intact — two mapping slots).
2. Mapping outputs still carry `.major.`/`.minor.` filename slots and `summarize.R` still produces populated `Major_reference`/`Minor_reference` columns (shim integrity).
3. The legacy `major_*`/`minor_*`/`minor_call`/`gate_flag` columns are present in the parsefirstmapping CSV.
4. A new per-candidate `confirmation_status` field exists in the long-format table.
5. For an unambiguous single-subtype-dominant fixture, the selected `cand_1` == today's major. Cases where `is_valid_minor()` previously suppressed a minor are EXPECTED to differ — **do not assert golden reproduction** (Phase 9 / COMPAT-01 gate).

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| `params.n_candidates` validation surfaces a clear error on invalid input (e.g. 0 or non-integer) | REFSEL-02 | nf-schema param-validation error path is awkward to snapshot | `nextflow run . -profile test,docker --n_candidates 0` → expect validation failure with a bounds message |

---

## Validation Sign-Off

- [x] All tasks have `<automated>` verify or Wave 0 dependencies
- [x] Sampling continuity: no 3 consecutive tasks without automated verify
- [x] Wave 0 covers all MISSING references (authored within plan tasks)
- [x] No watch-mode flags
- [x] Feedback latency acceptable (R tests fast; full nf-test gated at wave boundaries)
- [x] `nyquist_compliant: true` set in frontmatter

**Approval:** approved 2026-06-12
