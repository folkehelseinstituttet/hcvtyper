---
phase: 9
slug: compatibility-filename-migration-regression-suite
status: complete
nyquist_compliant: true
wave_0_complete: true
created: 2026-06-14
audited: 2026-06-14
---

# Phase 9 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework (R)** | base R asserts via `fail()`/`ok()` helpers + `run_all.sh` accumulator (custom, no testthat) |
| **Framework (Nextflow)** | nf-test 0.9.2 |
| **Config file** | `nf-test.config` (root); `bin/tests/run_all.sh` for R |
| **Quick run command** | `docker run --rm -v "$PWD":/work -w /work community.wave.seqera.io/library/r-seqinr_r-tidyverse:5358395134867368 Rscript bin/tests/test_compat.R` |
| **Full suite command** | `docker run --rm -v "$PWD":/work -w /work community.wave.seqera.io/library/r-seqinr_r-tidyverse:5358395134867368 bash bin/tests/run_all.sh` |
| **Estimated runtime** | ~60 seconds (R suite); ~5–10 min (full nf-test with Docker) |

---

## Sampling Rate

- **After every task commit:** Run quick R command above for R changes; `conda run -n NEXTFLOW nf-test test <changed module> --profile test,docker` for module changes
- **After every plan wave:** Run full suite command above + both module nf-test runs
- **Before `/gsd-verify-work`:** Full suite must be green + both `.snap` files regenerated
- **Max feedback latency:** ~60 seconds (R only) / ~10 minutes (nf-test modules)

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 09-01-01 | 01 | 1 | COMPAT-02 | — | N/A | config edit | `conda run -n NEXTFLOW nf-test test modules/local/parsefirstmapping/tests/main.nf.test --profile test,docker` | ✅ snap has cand1/cand2 | ✅ green |
| 09-01-02 | 01 | 1 | COMPAT-02 | — | N/A | module edit | `conda run -n NEXTFLOW nf-test test modules/local/parsefirstmapping/tests/main.nf.test --profile test,docker` | ✅ snap has cand1/cand2 | ✅ green |
| 09-01-03 | 01 | 1 | COMPAT-02 | — | N/A | workflow edit | visual diff of hcvtyper.nf fan-out + nf-test | ✅ (manual-only) | ✅ manual |
| 09-01-04 | 01 | 1 | COMPAT-02 | — | N/A | R refactor | `docker run --rm -v "$PWD":/work -w /work community.wave.seqera.io/library/r-seqinr_r-tidyverse:5358395134867368 Rscript bin/tests/test_compat.R` | ✅ bin/tests/test_compat.R | ✅ green |
| 09-01-05 | 01 | 1 | COMPAT-02 | — | N/A | snap regen | `conda run -n NEXTFLOW nf-test test modules/local/parsefirstmapping/tests/main.nf.test --profile test,docker --update-snapshot` | ✅ snap regenerated | ✅ green |
| 09-02-01 | 02 | 2 | COMPAT-01 | — | N/A | unit (subprocess) | `docker run --rm -v "$PWD":/work -w /work community.wave.seqera.io/library/r-seqinr_r-tidyverse:5358395134867368 Rscript bin/tests/test_compat.R` | ✅ bin/tests/test_compat.R | ✅ green |
| 09-02-02 | 02 | 2 | COMPAT-03 | — | N/A | unit (column-presence) | `docker run --rm -v "$PWD":/work -w /work community.wave.seqera.io/library/r-seqinr_r-tidyverse:5358395134867368 Rscript bin/tests/test_compat.R` | ✅ bin/tests/test_compat.R | ✅ green |
| 09-02-03 | 02 | 2 | COMPAT-04 | — | N/A | integration (function) | `docker run --rm -v "$PWD":/work -w /work community.wave.seqera.io/library/r-seqinr_r-tidyverse:5358395134867368 Rscript bin/tests/test_compat.R` | ✅ bin/tests/test_compat.R | ✅ green |
| 09-02-04 | 02 | 2 | TEST-01 | — | N/A | full R suite | `docker run --rm -v "$PWD":/work -w /work community.wave.seqera.io/library/r-seqinr_r-tidyverse:5358395134867368 bash bin/tests/run_all.sh` | ✅ bin/tests/run_all.sh | ✅ green |
| 09-03-01 | 03 | 3 | COMPAT-02 | — | N/A | R integration (candidate_rank join in stats loops) | `docker run --rm -v "$PWD":/work -w /work community.wave.seqera.io/library/r-seqinr_r-tidyverse:5358395134867368 Rscript bin/tests/test_compat.R` | ✅ COMPAT-02 smoke in test_compat.R | ✅ green |
| 09-03-02 | 03 | 3 | COMPAT-02 | — | N/A | R integration (cv_by_ref cand-slot strip) | `docker run --rm -v "$PWD":/work -w /work community.wave.seqera.io/library/r-seqinr_r-tidyverse:5358395134867368 Rscript bin/tests/test_compat.R` | ✅ COMPAT-02 smoke asserts suffix-free refs | ✅ green |
| 09-04-01 | 04 | 4 | COMPAT-01 | — | N/A | checkpoint (human-verify golden fixture) | Manual — human confirms compat_golden.csv values before assertions locked | ✅ (manual-only) | ✅ manual |
| 09-04-02 | 04 | 4 | COMPAT-01/02/03 | — | N/A | subprocess (summarize.R end-to-end) | `docker run --rm -v "$PWD":/work -w /work community.wave.seqera.io/library/r-seqinr_r-tidyverse:5358395134867368 Rscript bin/tests/test_compat.R` | ✅ bin/tests/test_compat.R | ✅ green |
| 09-04-03 | 04 | 4 | COMPAT-04 | — | N/A | TDD integration (1a/1b allowed + 2k1b suppressed, subprocess + helper) | `docker run --rm -v "$PWD":/work -w /work community.wave.seqera.io/library/r-seqinr_r-tidyverse:5358395134867368 Rscript bin/tests/test_compat.R` | ✅ bin/tests/test_compat.R | ✅ green |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky · ✅ manual (manual-only verification)*

---

## Wave 0 Requirements

- [x] `bin/tests/test_compat.R` — new file covering COMPAT-01..04 (D-08) ✅ committed feat(09-04)
- [x] `bin/tests/fixtures/compat_golden.csv` — two D-07 golden cases (MONO monoinfection + COINF co-infection) ✅ committed test(09-04)
- [x] Regenerate `modules/local/parsefirstmapping/tests/main.nf.test.snap` via `--update-snapshot` ✅ snap has cand1/cand2 slots
- [x] Regenerate `modules/local/blastparse/tests/main.nf.test.snap` via `--update-snapshot` ✅ snap exists
- [x] Framework install: none — R container + nf-test already provisioned ✅

*All Wave 0 items complete as of 2026-06-14.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| COMPAT-01 golden fixture values match true v1.0/v2.0 run | COMPAT-01 | Fixture accuracy depends on which run anchored `flagoff_golden.csv` (A2) — automated test asserts fixture values but cannot verify the fixture itself is correct | Compare D-07 expected values (`Major_reference`, `Minor_reference`, `overall_sample_call`) against the v1.0/v2.0 handoff evidence table or a reference run before locking test assertions (`checkpoint:human-verify`) |
| nf-test snap diff covers only `.major.fa`/`.minor.fa` → `.cand1.fa`/`.cand2.fa` changes | COMPAT-02 | `--update-snapshot` accepts ALL current output blindly; unexpected md5 churn = real regression | Eyeball snap diff before committing; ONLY expected change is the slot string in filename md5s; any other md5 change = investigate |
| Full `nextflow run -profile test,docker` confirms Summary.csv parses correctly | COMPAT-01/02 | Full pipeline end-to-end not captured by R unit tests alone | Run `conda run -n NEXTFLOW nextflow run . -profile test,docker` on the test dataset; inspect Summary.csv for populated `Major_reference`/`Minor_reference` with new slot names |

---

## Validation Sign-Off

- [x] All tasks have `<automated>` verify or Wave 0 dependencies
- [x] Sampling continuity: no 3 consecutive tasks without automated verify
- [x] Wave 0 covers all MISSING references
- [x] No watch-mode flags
- [x] Feedback latency < 60s (R), < 10min (nf-test)
- [x] `nyquist_compliant: true` set in frontmatter

**Approval:** 2026-06-14

---

## Validation Audit 2026-06-14

| Metric | Count |
|--------|-------|
| Gaps found | 0 |
| Resolved | 0 |
| Escalated (manual-only) | 2 |
| Tasks added (plans 03/04 missing from map) | 5 |

All requirements covered. Plans 03/04 tasks (09-03-01, 09-03-02, 09-04-01, 09-04-02, 09-04-03) added to Per-Task Map — all COVERED via `test_compat.R` or designated Manual-Only (human-verify checkpoint, visual snap diff).
