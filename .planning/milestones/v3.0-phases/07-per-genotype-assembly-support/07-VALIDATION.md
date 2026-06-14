---
phase: "07"
slug: per-genotype-assembly-support
status: complete
nyquist_compliant: true
wave_0_complete: true
created: 2026-06-14
---

# Phase 07 — Validation Strategy

> Per-phase validation contract reconstructed from Plan + Summary artifacts (State B).
> All tests verified green via Docker on 2026-06-14.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | R subprocess-contract + function-level tests (base R + tidyverse) |
| **Container** | `community.wave.seqera.io/library/r-seqinr_r-tidyverse:5358395134867368` |
| **Config file** | `.github/workflows/ci.yml` (Run R regression guard step) |
| **Quick run command** | `docker run --rm -v "$PWD":/work -w /work community.wave.seqera.io/library/r-seqinr_r-tidyverse:5358395134867368 Rscript bin/tests/test_assembly_support.R` |
| **Full suite command** | `docker run --rm -v "$PWD":/work -w /work community.wave.seqera.io/library/r-seqinr_r-tidyverse:5358395134867368 bash bin/tests/run_all.sh` |
| **Estimated runtime** | ~2–3 minutes |

---

## Sampling Rate

- **After every task commit:** Run quick run command (single test file for the task)
- **After every plan wave:** Run full suite command (`bash bin/tests/run_all.sh`)
- **Before `/gsd-verify-work`:** Full suite must be green
- **Max feedback latency:** ~180 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 07-01-01 | 01 | 1 | ASUP-01 | T-07-01, T-07-02 | Empty/no-hit input writes typed header-only CSV, exits 0; NA contig names do not abort | subprocess | `Rscript bin/tests/test_assembly_support.R` | ✅ | ✅ green |
| 07-01-02 | 01 | 1 | ASUP-01 | — | BLASTPARSE emits non-optional `*.assembly_support.csv`; stub writes real header | grep | `grep -q 'emit: support' modules/local/blastparse/main.nf && grep -c 'assembly_support.csv' modules/local/blastparse/main.nf` | ✅ | ✅ green |
| 07-01-03 | 01 | 1 | ASUP-01 | T-07-01 | Subprocess test pins Case A (single-best-by-length), Case B (multi-hit de-dup), Case C (empty no-abort) | subprocess | `Rscript bin/tests/test_assembly_support.R` | ✅ | ✅ green |
| 07-02-01 | 02 | 2 | ASUP-02 | T-07-03, T-07-05 | No candidate row dropped (left-join); unmatched gets `assembly_support="none"` + NA; zero-row inputs return typed frames | unit | `Rscript bin/tests/test_assembly_support_join.R` | ✅ | ✅ green |
| 07-02-02 | 02 | 2 | ASUP-02 | T-07-03, T-07-04 | summarize.R reads both CSVs with typed-empty fallbacks; legacy `apply_denovo_layer`/`review_flag` byte-unchanged (D-04) | parse+grep | `Rscript -e 'invisible(parse("bin/summarize.R"))' && grep -q 'join_assembly_support' bin/summarize.R && grep -q 'PARSEFIRSTMAPPING.out.candidates' workflows/hcvtyper.nf && echo WIRED_OK` | ✅ | ✅ green |
| 07-02-03 | 02 | 2 | ASUP-02 | T-07-05 | Criterion #2 (3b←3a genotype match), criterion #3 (no row loss + explicit none), criterion #4 (cand_2 equals today's denovo_minor_contig_length) | unit | `Rscript bin/tests/test_assembly_support_join.R` | ✅ | ✅ green |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

Existing infrastructure covers all phase requirements. No new test framework or fixture scaffolding was needed — the subprocess-contract pattern (`test_candidate_selection.R`) and function-level pattern (`test_summarize_denovo.R`) from prior phases were reused directly.

---

## Manual-Only Verifications

All phase behaviors have automated verification.

---

## Validation Sign-Off

- [x] All tasks have `<automated>` verify or Wave 0 dependencies
- [x] Sampling continuity: no 3 consecutive tasks without automated verify
- [x] Wave 0 covers all MISSING references (none — all COVERED)
- [x] No watch-mode flags
- [x] Feedback latency < 180s
- [x] `nyquist_compliant: true` set in frontmatter

**Approval:** approved 2026-06-14

---

## Validation Audit 2026-06-14

| Metric | Count |
|--------|-------|
| Gaps found | 0 |
| Resolved | 0 |
| Escalated | 0 |
| Status | NYQUIST-COMPLIANT (all tasks COVERED) |

Full suite result: `ALL R TESTS PASSED` (10 test files, 0 failures) verified via Docker on 2026-06-14.
