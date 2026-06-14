---
phase: 06-neutral-candidate-selection
verified: 2026-06-12T00:00:00Z
status: passed
score: 3/3 must-haves verified
overrides_applied: 0
---

# Phase 6: Neutral Candidate Selection — Verification Report

**Phase Goal:** Reference selection during the run becomes dominance-neutral — it ranks N candidate references by read recruitment and maps each one independently — so no major/minor semantics are baked in before the evidence is in. This is the structural foundation the later classification depends on.
**Verified:** 2026-06-12
**Status:** PASSED
**Re-verification:** No — initial verification

---

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | After first-pass mapping, candidate references are emitted as a neutrally-ranked ordered set (cand_1..cand_n) ordered by read recruitment, with no "major"/"minor" label assigned during the run | VERIFIED | `bin/summarize_mapping_to_all_references.R` ranks by `arrange(desc(reads)) %>% head(n = n_candidates) %>% mutate(candidate_rank = row_number())`; `is_valid_minor()` deleted; `*.candidates.csv` long-format emitted; all R tests pass (`bash bin/tests/run_all.sh` → ALL PASS) |
| 2 | The number of candidates is controlled by a single parameter (default 2); running with the default produces the same two candidate slots the pipeline selects today | VERIFIED | `params.n_candidates` declared as `integer, default 2, minimum 1` in all three required files (`nextflow_schema.json`, `nextflow.config`, `conf/modules_hcv.config`); R script accepts it as positional arg 7 with `as.integer` coercion and NA-safe default; module passes `${params.n_candidates}` positionally; at N=2 the snapshot produces `Test_2.3a_D17763_major.major.nodup.bam` and a populated `Summary.csv` |
| 3 | Per-candidate `confirmation_status` replaces the old `gate_flag`/`minor_call` (REFSEL-03), and each candidate is mapped independently via a uniform TARGETED_MAPPING fan-out | VERIFIED | `MAJOR_MAPPING`/`MINOR_MAPPING` aliases absent from `workflows/hcvtyper.nf` and `conf/modules_hcv.config`; single `TARGETED_MAPPING` import and `splitCsv(header:true).flatMap` fan-out over `PARSEFIRSTMAPPING.out.candidates`; `confirmation_status` carried in `new_meta` per candidate; all six downstream feeds (HCVGLUE aligned + 5 SUMMARIZE channels) source from `TARGETED_MAPPING.out.*`; workflow snapshot regenerated green (2/2 PASSED, 75 tasks, unchanged) |

**Score:** 3/3 truths verified

---

## Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `bin/summarize_mapping_to_all_references.R` | Neutral top-N candidate ranking + long-format table + legacy shim reconstruction | VERIFIED | Exists, substantive (neutral ranking logic at lines 137-162, `confirmation_status` at lines 145-160, legacy shim at lines 167-201, FASTA write at lines 210-221), no `is_valid_minor`, sourced from `genotype_utils.R` via relative path |
| `bin/tests/test_candidate_selection.R` | R subprocess-contract test for ranking, distinct-subtype dedup, shim column reconstruction, NA-fill | VERIFIED | Exists, parses as valid R, asserts `candidate_rank` and `confirmation_status`, covers ranking/dedup/single-candidate/empty-input cases, auto-discovered by `run_all.sh` glob |
| `nextflow_schema.json` | Typed `n_candidates` param declaration | VERIFIED | Contains `"n_candidates": {"type": "integer", "default": 2, "minimum": 1, "description": "..."}` |
| `nextflow.config` | `n_candidates` default 2 | VERIFIED | Line 53: `n_candidates = 2 // Number of neutrally-ranked candidate references...` |
| `conf/modules_hcv.config` | `n_candidates` default in params{} block + slot-from-rank ext.prefix closures | VERIFIED | Line 21: `n_candidates = 2`; lines 130-263: single `TARGETED_MAPPING:<PROCESS>` closures deriving slot from `meta.candidate_rank.toInteger()` |
| `modules/local/parsefirstmapping/main.nf` | `candidates` emit + `n_candidates` positional arg + updated stub | VERIFIED | Line 28: `tuple val(meta), path("*.candidates.csv"), emit: candidates, optional: true`; line 46: `${params.n_candidates}`; legacy globs pinned to `*.parsefirstmapping.csv`; stub writes 2-row candidates CSV with correct header |
| `modules/local/parsefirstmapping/tests/main.nf.test` | Updated nf-test assertions for the candidates emit | VERIFIED | All 5 test cases assert `process.out.candidates` and check long-format header contains `candidate_rank`/`confirmation_status` |
| `modules/local/parsefirstmapping/tests/main.nf.test.snap` | Regenerated snapshot for new emit contract | VERIFIED | Exists (8058 bytes, 2026-06-12), contains `candidates` entries with `*.candidates.csv:md5,...` for all 5 fixture samples |
| `workflows/hcvtyper.nf` | Uniform per-candidate fan-out replacing the two alias subworkflow calls | VERIFIED | No `MAJOR_MAPPING`/`MINOR_MAPPING` references; single `TARGETED_MAPPING` import; `splitCsv.flatMap` fan-out at lines 381-436; all feeds at lines 446, 474-496 |
| `tests/default.nf.test.snap` | Regenerated workflow snapshot | VERIFIED | Exists, contains `.major.` filename entries (`Test_2.3a_D17763_major.major.nodup.bam`), `*.candidates.csv` entries, `Summary.csv` entry — 75 tasks, 2/2 PASSED |

---

## Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `bin/summarize_mapping_to_all_references.R` | `bin/genotype_utils.R` | `source("genotype_utils.R")` relative staged path | WIRED | Line 10: `source("genotype_utils.R")`; `genotype_from_subtype()` called at line 69 |
| `bin/tests/run_all.sh` | `bin/tests/test_candidate_selection.R` | `test_*.R` glob | WIRED | `run_all.sh` line 19: `for t in "$here"/test_*.R`; file auto-discovered and confirmed green |
| `modules/local/parsefirstmapping/main.nf` | `bin/summarize_mapping_to_all_references.R` | `${params.n_candidates}` positional arg | WIRED | Line 46: `${params.n_candidates}` appended after `${params.minCov}`; R script reads it as `args[7]` |
| `modules/local/parsefirstmapping/main.nf` | `workflows/hcvtyper.nf` fan-out | `emit: candidates` consumed by `splitCsv` | WIRED | Module declares `emit: candidates`; workflow at line 396: `PARSEFIRSTMAPPING.out.candidates` then `.join(...).flatMap{... .splitCsv(...)}` |
| `workflows/hcvtyper.nf` | `PARSEFIRSTMAPPING.out.candidates` | `splitCsv(header:true)` per-candidate fan-out | WIRED | Line 402: `candidates_csv.splitCsv(header: true, sep:',')` iterates ALL rows |
| `conf/modules_hcv.config` ext.prefix | `bin/summarize.R` filename 3rd-dot-field parse | `meta.candidate_rank == 1 ? 'major' : 'minor'` slot literal | WIRED | Config closures at lines 131, 228, 240, 248, 257, 263 all derive slot via `.toInteger() == 1 ? 'major' : 'minor'`; snapshot confirms `Test_2.3a_D17763_major.major.nodup.bam` is byte-identical to pre-change baseline |
| `workflows/hcvtyper.nf TARGETED_MAPPING.out.*` | HCVGLUE + SUMMARIZE inputs | Single subworkflow output replacing `.mix(MAJOR..., MINOR...)` | WIRED | `TARGETED_MAPPING.out` appears 8 times in `hcvtyper.nf`; lines 446, 474, 475, 476, 495, 496 cover all six required feeds |

---

## Data-Flow Trace (Level 4)

| Artifact | Data Variable | Source | Produces Real Data | Status |
|----------|---------------|--------|--------------------|--------|
| `bin/summarize_mapping_to_all_references.R` | `candidates_long` | `read_tsv(idxstats)` + `group_by(Subtype)` ranking | Yes — reads from process-computed idxstats, produces per-candidate rows | FLOWING |
| `modules/local/parsefirstmapping/main.nf` | `*.candidates.csv` | R script invocation via `Rscript` in `script:` block | Yes — R script writes the file; `emit: candidates` captures it | FLOWING |
| `workflows/hcvtyper.nf` fan-out | `new_meta` per-candidate channel element | `PARSEFIRSTMAPPING.out.candidates` → `splitCsv.flatMap` | Yes — each CSV row becomes a live channel element with `confirmation_status` | FLOWING |

---

## Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| R selection script parses as valid R | `Rscript -e 'invisible(parse("bin/summarize_mapping_to_all_references.R"))'` | exit 0 | PASS |
| `is_valid_minor` deleted from R script | `grep -q 'is_valid_minor' bin/summarize_mapping_to_all_references.R` | no match (exit 1 for grep = deleted) | PASS |
| R unit test suite passes | `bash bin/tests/run_all.sh` | "ALL R TESTS PASSED" | PASS |
| `n_candidates` in all three param files | `grep -q '"n_candidates"' nextflow_schema.json && grep -q 'n_candidates' nextflow.config && grep -q 'n_candidates' conf/modules_hcv.config` | All present | PASS |
| `emit: candidates` in PARSEFIRSTMAPPING | `grep -q 'emit: candidates' modules/local/parsefirstmapping/main.nf` | Found at line 28 | PASS |
| Legacy emits retained | `grep -q 'emit: major_mapping' && grep -q 'emit: minor_mapping'` in `main.nf` | Both present, globs pinned to `*.parsefirstmapping.csv` | PASS |
| No `MAJOR_MAPPING`/`MINOR_MAPPING` in workflow | `grep -q 'MAJOR_MAPPING\|MINOR_MAPPING' workflows/hcvtyper.nf` | No match | PASS |
| `splitCsv` and `candidate_rank` in workflow | Grep confirms | Both present | PASS |
| No `MAJOR_MAPPING`/`MINOR_MAPPING` in config | `grep -q 'MAJOR_MAPPING\|MINOR_MAPPING' conf/modules_hcv.config` | No match | PASS |
| Slot-from-rank in config | `grep -q 'candidate_rank' conf/modules_hcv.config` | 8 closures deriving slot | PASS |
| `TARGETED_MAPPING.out` feeds (>=5) | `grep -c 'TARGETED_MAPPING.out' workflows/hcvtyper.nf` | 8 | PASS |
| All 9 task commits exist in git | `git log --oneline -12` | 52fda8d, ea71811, 43904de, f2c2560, d8751a7, 3bc6635, d534cba, acdf8c0 all present (plus a post-commit fix 6c0baf9) | PASS |
| Workflow snapshot reflects `.major.` filenames | Grep `default.nf.test.snap` for `major` | `Test_2.3a_D17763_major.major.nodup.bam` etc. present | PASS |

---

## Probe Execution

Step 7c: No `scripts/*/tests/probe-*.sh` files exist; phase is a Nextflow/R refactor, not a CLI/migration phase. nf-test (module + workflow) was run by the executor and the regenerated snapshots are committed and verified above.

---

## Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|------------|-------------|--------|----------|
| REFSEL-01 | 06-01, 06-02, 06-03 | Reference selection ranks candidates neutrally (cand_1..cand_n) by read recruitment, no major/minor dominance semantics during the run | SATISFIED | Neutral read-recruitment ranking in R script; `is_valid_minor` deleted; `candidates` channel emitted; workflow routes on `confirmation_status == 'pass'` string |
| REFSEL-02 | 06-01, 06-02 | Candidate count is a parameter, default 2 (reproduces today's two-slot behaviour) | SATISFIED | `params.n_candidates` declared as `integer, default 2, minimum 1` across typed-param triple; module passes `${params.n_candidates}` positionally |
| REFSEL-03 | 06-03 | Each candidate is independently targeted-mapped, replacing asymmetric MAJOR/MINOR_MAPPING aliases; `gate_flag`/`minor_call` plumbing becomes per-candidate `confirmation_status` | SATISFIED | Single `TARGETED_MAPPING` fan-out via `splitCsv.flatMap`; `confirmation_status` in per-candidate meta; all six downstream feeds relinked in lockstep |

All three phase requirements are satisfied. No orphaned requirements for Phase 6 were found in REQUIREMENTS.md.

---

## Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| No TBD/FIXME/XXX markers found in phase-modified files | — | — | — | — |

No blocking debt markers detected. The code-review file (`06-REVIEW.md`) documents CR-01 (n_candidates >= 3 wrong-FASTA selection) and WR-01 through WR-05. These are carry-forward findings:

- **CR-01** (n >= 3 wrong FASTA): Confirmed in codebase at `workflows/hcvtyper.nf` line 418-423 — the binary `rank == '1' ? major_fasta : minor_fasta` FASTA selection produces wrong-reference mapping for rank-3 candidates. Per the verified scoping note, this is the Phase-9 boundary (COMPAT-02). The default is `n_candidates = 2` and there is no user-facing documentation advertising n >= 3 as supported. This is NOTED as a carry-forward design boundary, not a Phase-6 success-criterion violation.
- **WR-01** (tie non-determinism), **WR-02** (major-pass comment overstates guarantee), **WR-03** (n >= 3 GLUE/summarize drop): All relate to the n >= 3 unsupported path or future classification work; none affect the n=2 success criteria.
- **WR-04 / WR-05** (`1:length()` in `summarize.R`): Pre-existing; not introduced by this phase.

---

## Human Verification Required

None. All phase success criteria are verifiable from the codebase alone:

- Parameter declarations: checked via grep and JSON parse.
- R script logic: confirmed by parsing and unit-test results.
- Module emit contract: confirmed by grep and nf-test snapshot.
- Workflow routing topology: confirmed by grep and workflow snapshot file content.
- Shim integrity at N=2: confirmed by snapshot entries showing `.major.` filenames and `Summary.csv`.

---

## Gaps Summary

No gaps. All three success criteria (REFSEL-01, REFSEL-02, REFSEL-03) are fully implemented and verified:

1. Neutral read-recruitment ranking with no validity filtering is in the R script and confirmed by the passing R unit test suite.
2. `params.n_candidates` (integer, default 2, minimum 1) is declared in all three required locations and wired into the module.
3. The uniform per-candidate `TARGETED_MAPPING` fan-out with `confirmation_status` per candidate replaces the asymmetric aliases; the `.major.`/`.minor.` shim holds at N=2; all downstream feed channels are relinked.

The n >= 3 FASTA-selection defect (CR-01 in 06-REVIEW.md) is correctly scoped as a deferred Phase-9 concern and does not violate any Phase-6 success criterion.

---

_Verified: 2026-06-12_
_Verifier: Claude (gsd-verifier)_
