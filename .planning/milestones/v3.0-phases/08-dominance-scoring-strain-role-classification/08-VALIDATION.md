---
phase: 8
slug: dominance-scoring-strain-role-classification
status: complete
nyquist_compliant: true
wave_0_complete: true
created: 2026-06-13
audited: 2026-06-14
---

# Phase 8 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Rscript unit tests in `bin/tests/` (testthat-style `ok()`/`fail()` assertions, exit 1 on fail) |
| **Test runner** | `bash bin/tests/run_all.sh` (auto-globs `test_*.R`) |
| **Quick run command** | `docker run --rm -v "$PWD":/work -w /work jonbra/tidyverse_seqinr:2.0 bash bin/tests/run_all.sh` |
| **Full suite command** | same — all tests are R unit tests in `bin/tests/` |
| **Estimated runtime** | ~30–60 seconds |
| **Note** | Host R (4.3.3) lacks tidyverse; suite must run in `jonbra/tidyverse_seqinr:2.0` container (or `community.wave.seqera.io/library/r-seqinr_r-tidyverse`) |

---

## Sampling Rate

- **After every task commit:** Run the quick run command for the touched module/helper
- **After every plan wave:** Run the full suite command
- **Before `/gsd-verify-work`:** Full suite must be green
- **Max feedback latency:** 180 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 8-01-T1 | 01 | 1 | SCORE-01 / SCORE-02 (config baseline) | — | N/A | grep verify | `grep -E 'denovo_min_contig_length\s*=\s*1000\|denovo_min_kmer_cov\s*=\s*2\.0' nextflow.config` | ✅ | ✅ green |
| 8-01-T2 | 01 | 1 | SCORE-01 / SCORE-02 / CLASS-01..04 (helper) | T-08-01 | Typed zero-row guard — returns typed frame, never stop() | unit | `docker run --rm -v "$PWD":/work -w /work jonbra/tidyverse_seqinr:2.0 Rscript bin/tests/test_dominance_score.R` | ✅ `bin/classify_roles.R` | ✅ green |
| 8-01-T3 | 01 | 1 | SCORE-01 / SCORE-02 / CLASS-01..04 (unit tests) | T-08-01 | Zero-row / NULL → typed frame (Test 9) | unit | `docker run --rm -v "$PWD":/work -w /work jonbra/tidyverse_seqinr:2.0 bash bin/tests/run_all.sh` | ✅ `bin/tests/test_dominance_score.R`, `bin/tests/test_classify_roles.R` | ✅ green |
| 8-02-T1 | 02 | 2 | SCORE-01 / SCORE-02 (ext.args wiring) | T-08-04 | score_weight_* appended AFTER n_candidates — no positional re-map | grep verify | `grep -E 'n_candidates.*score_weight_evenness' conf/modules_hcv.config` | ✅ `conf/modules_hcv.config` | ✅ green |
| 8-02-T2 | 02 | 2 | CLASS-01..04 / D-15 retirement | T-08-05 | Legacy apply_denovo_layer retired; classify_roles() runs; enriched candidates.csv written | unit + grep | `docker run --rm -v "$PWD":/work -w /work jonbra/tidyverse_seqinr:2.0 Rscript bin/tests/test_summarize_denovo.R` | ✅ `bin/tests/test_summarize_denovo.R` | ✅ green |
| 8-02-T3 | 02 | 2 | CLASS-01..04 (module staging) | — | classify_roles.R staged; stub schema updated (overall_sample_call in, minor_denovo_status out) | grep verify | `grep -q 'classify_roles' modules/local/summarize/main.nf && grep -q 'overall_sample_call' modules/local/summarize/main.nf && ! grep -q 'minor_denovo_status' modules/local/summarize/main.nf` | ✅ `modules/local/summarize/main.nf` | ✅ green |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [x] Function test(s) for the dominance-score helper (SCORE-01, SCORE-02) — breadth-evenness weighting vs raw read count → `bin/tests/test_dominance_score.R` (Tests 1–5, incl. 4g HEADLINE)
- [x] Function test(s) for the strain-role classifier (CLASS-01..CLASS-04) — including the benchmarked false-4g-minor (→ background) and genuine partial co-infection (→ preserved) cases → `bin/tests/test_classify_roles.R` (Tests 1–9)
- [ ] SUMMARIZE module snapshot regen for the Summary.csv schema change (overall_sample_call + role columns in, minor_denovo_status out) — **deferred**: folded into the already-pending Phase-7 SUMMARIZE nf-test snapshot regen (docker/disk constrained); NOT a Phase-8 gate per CONTEXT scope

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| End-to-end overall sample call on real benchmark FASTQs | CLASS-04 | Full pipeline run is heavy (docker/disk constrained) | Run pipeline on the false-4g-minor and genuine-coinfection samples; confirm roles + overall call in Summary output |
| cv_evenness = 0 when reference has zero mean depth (T-08-06) | T-08-06 | Guard is inside the cov loop over pipeline depth files; requires real SAMTOOLS_DEPTH output | Run SUMMARIZE step on a no-mapping sample; verify cv_evenness column is 0.0, not NaN/Inf, in the candidate frame |

---

## Validation Sign-Off

- [x] All tasks have `<automated>` verify or Wave 0 dependencies
- [x] Sampling continuity: no 3 consecutive tasks without automated verify
- [x] Wave 0 automated requirements covered
- [x] No watch-mode flags
- [x] Feedback latency < 180s
- [x] `nyquist_compliant: true` set in frontmatter

**Approval:** complete

---

## Validation Audit 2026-06-14

| Metric | Count |
|--------|-------|
| Tasks audited | 6 |
| Gaps found | 0 |
| Requirements COVERED | 6 (SCORE-01, SCORE-02, CLASS-01, CLASS-02, CLASS-03, CLASS-04) |
| Requirements PARTIAL | 0 |
| Requirements MISSING | 0 |
| Manual-only added | 1 (T-08-06 cv_evenness zero-mean guard) |
| Suite result | ALL PASS (`bash bin/tests/run_all.sh` in `jonbra/tidyverse_seqinr:2.0`) |

**Test files confirmed green:**
- `bin/tests/test_dominance_score.R` — SCORE-01/02 (5 tests incl. 4g HEADLINE)
- `bin/tests/test_classify_roles.R` — CLASS-01..04 + D-11/D-12 (9 tests)
- `bin/tests/test_summarize_denovo.R` — D-15 retirement + integration contract
- All pre-existing tests (test_assembly_support_join, test_assembly_support, test_candidate_selection, test_coinfection, test_compat, test_denovo_confirm, test_major_gate) — ALL PASS

**Note:** VALIDATION.md was in draft state (pre-execution) when audited. Updated to reflect completed execution: test files exist, suite is green, all Wave-0 requirements fulfilled. One deferred item (nf-test SUMMARIZE snapshot regen) remains tracked from Phase-7 and is not a Phase-8 gate.
