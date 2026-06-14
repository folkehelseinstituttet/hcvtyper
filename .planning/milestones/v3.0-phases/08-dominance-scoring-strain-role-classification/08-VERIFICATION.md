---
phase: 08-dominance-scoring-strain-role-classification
verified: 2026-06-14T12:00:00Z
status: human_needed
score: 13/13
overrides_applied: 0
human_verification:
  - test: "Run the full test suite in the Docker container"
    expected: "bash bin/tests/run_all.sh exits 0 with PASS lines for test_dominance_score.R and test_classify_roles.R including the 4g-headline and all evidence-table cases"
    why_human: "tidyverse is not installed in the host shell — tests require the pipeline Docker container (community.wave.seqera.io/library/r-gridextra_r-png_r-seqinr_r-tidyverse:3536dd50a17de0ab). The R parse check passes (PARSE_OK) and all logic was verified by static analysis; this is a runtime confirmation, not a code gap."
---

# Phase 08: Dominance Scoring + Strain-Role Classification — Verification Report

**Phase Goal:** Replace the 2-strain major/minor confirmation model with an N-candidate dominance scoring + strain-role classification model. Every candidate (dominant / co-infection / background) gets a combined dominance score (breadth-evenness weighted) and a role. Background/artefact candidates are surfaced explicitly, never dropped. One overall_sample_call is reported per sample.

**Verified:** 2026-06-14T12:00:00Z
**Status:** human_needed
**Re-verification:** No — initial verification

---

## Goal Achievement

### Observable Truths

Truths are derived from the ROADMAP.md Success Criteria (5 items) plus the must-haves from both PLAN frontmatter files. All are verified against the actual codebase.

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | Each candidate gets a combined dominance score over mapped read count, k-mer coverage, and mapping coverage breadth, with breadth-evenness weighted strongly above raw read count | VERIFIED | `classify_roles.R:164–168`: `score = wr*reads_term + we*breadth_frac + we*even_fac + wk*kmer_term` with defaults `we=3.0 > wr=1.0`; cv_evenness computed in `summarize.R:413–414` cov loop from `cov$X3` before it is discarded |
| 2 | Each candidate is classified into exactly one strain role — dominant / co-infection / background | VERIFIED | `classify_roles.R:262–312`: every code path (dominant, below_floor, corroborated, refuted_denovo, uncorroborated_kept) assigns exactly one of the three roles; `test_classify_roles.R:199–205` asserts no NA roles and only the three valid values (CLASS-01) |
| 3 | A non-dominant candidate clearing the floor with assembly support classifies co-infection/corroborated | VERIFIED | `classify_roles.R:277–279`: `own_substantial TRUE → role="co-infection", reason="corroborated"`; test_classify_roles.R Test 2 exercises ERR1810447 full 9207bp 2b and ERR1810453 partial 2b at kmer~5 (passes reconciled 2.0 floor) |
| 4 | A non-dominant candidate clearing the floor with NO support, when the dominant DID assemble substantially, classifies background/refuted_denovo | VERIFIED | `classify_roles.R:280–283`: `dom_substantial TRUE → role="background", reason="refuted_denovo"`; test_classify_roles.R Test 1 (false 4g case) |
| 5 | A non-dominant candidate clearing the floor, no own support, dominant ALSO did not assemble substantially → co-infection/uncorroborated_kept | VERIFIED | `classify_roles.R:284–289`: `!dom_substantial → role="co-infection", reason="uncorroborated_kept"`; test_classify_roles.R Test 4 (IVT extreme-ratio case) |
| 6 | Same-genotype (non-1a/1b) candidate or 2k1b pair classifies background per ported is_valid_minor() | VERIFIED | `classify_roles.R:295–308`: D-12 exceptions applied after provisional corroboration; `is_valid_minor()` (L69–88) preserves all three original rules; test_classify_roles.R Tests 5 and 6 cover both demotion paths |
| 7 | An overall sample call (monoinfection / co-infection / indeterminate) is derived from roles and reported once per sample | VERIFIED | `classify_roles.R:315–321`: D-14 overall call; `summarize.R:1279,1325`: `overall_sample_call` in `pmap_chr` list and `select()` column order; module stub header L74/L79 confirms field presence |
| 8 | Background/artefact candidates are surfaced in the output with explicit role_reason, never silently dropped | VERIFIED | `summarize.R:671–678`: `write_csv(candidate_support, file="candidates.csv")` writes ALL candidates including background; `modules/local/summarize/main.nf:38`: `emit: candidates`; `conf/modules_hcv.config:343–347`: `publishDir` for `candidates.csv` |
| 9 | summarize.R sources classify_roles.R and runs score_candidates() + classify_roles() over the Phase-7-joined candidate frame | VERIFIED | `summarize.R:25`: `source("classify_roles.R")`; `summarize.R:645–669`: `score_candidates()` then `classify_roles()` called on `candidate_support` (the Phase-7-joined frame) |
| 10 | cv_evenness factor is computed inside the cov loop from per-position depth (cov$X3) before X3 is discarded | VERIFIED | `summarize.R:380–381`: `tmp_df` widened to ncol=8 with `cv_evenness` in colnames; `summarize.R:409–414`: `cv_raw = sd(cov$X3)/mean(cov$X3)` computed while `cov` is in scope, with zero-mean guard → `tmp_df$cv_evenness[i]` |
| 11 | score_weight_* params reach summarize.R via ext.args appended AFTER n_candidates (no positional re-map) | VERIFIED | `conf/modules_hcv.config:328`: ext.args ends `...${params.n_candidates} ${params.score_weight_evenness} ${params.score_weight_reads} ${params.score_weight_kmercov} ${params.score_evenness_k}`; `summarize.R:54,63–68`: n_candidates at args[11], weights at args[12–15] |
| 12 | The legacy apply_denovo_layer / classify_minor_denovo / coinfection_flag path is retired | VERIFIED | No uncommented `apply_denovo_layer(` call exists in summarize.R (grep confirms zero matches); `minor_denovo_status` and `coinfection_flag` appear only in comments at summarize.R:639,1123,1226,1254,1324 — no live code paths; `select()` at summarize.R:1311+ does not list either column |
| 13 | classify_roles.R is staged as a path() input to SUMMARIZE and passed at the workflow call site | VERIFIED | `modules/local/summarize/main.nf:34`: `path(classify_roles)`; `workflows/hcvtyper.nf:532`: `file("${projectDir}/bin/classify_roles.R")` |

**Score:** 13/13 truths verified

---

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `bin/classify_roles.R` | Pure sourced helper: score_candidates() + classify_roles() + is_valid_minor() | VERIFIED | All three functions present (L69, L114, L186); no commandArgs, no top-level file I/O; defensive `if (!exists("group_by"))` guard (L54); `stopifnot(match_level %in% c("genotype","subtype"))` (L193); `Rscript -e 'parse("bin/classify_roles.R")'` exits 0 |
| `bin/tests/test_dominance_score.R` | SCORE-01/02 unit coverage incl. the 4g-loses-to-even-minor headline assertion | VERIFIED | 335 lines; 5 test assertions including Test 4 headline (false 4g at 53279 reads vs genuine minor at 4199 reads) and Test 4b (no-kmer variant); sources real classify_roles.R (L29) |
| `bin/tests/test_classify_roles.R` | CLASS-01..04 unit coverage incl. false-4g→background and genuine-coinfection→preserved | VERIFIED | 220 lines; 9 tests covering all must-have scenarios; sources real classify_roles.R (L34); tests exercises null input (Test 9) |

---

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `bin/tests/test_classify_roles.R` | `bin/classify_roles.R` | `source()` | WIRED | L33–34: `source(file.path(bin_dir, "genotype_utils.R"))` then `source(file.path(bin_dir, "classify_roles.R"))` |
| `bin/tests/test_dominance_score.R` | `bin/classify_roles.R` | `source()` | WIRED | L28–29: same self-locating source pattern |
| `bin/classify_roles.R` | `bin/genotype_utils.R` | `genotype_from_subtype()` | WIRED | `classify_roles.R:70–73` alias vars; `test_classify_roles.R:52`: `candidate_genotype = genotype_from_subtype(subtype)` uses it in test fixture builder |
| `bin/summarize.R` | `bin/classify_roles.R` | `source("classify_roles.R")` then calls | WIRED | L25: source; L645: `score_candidates()`; L655: `classify_roles()` |
| `conf/modules_hcv.config` | `bin/summarize.R` | positional ext.args | WIRED | L328: score_weight_* appended after n_candidates; summarize.R L63–68 parses at args[12..15] |
| `modules/local/summarize/main.nf` | `bin/classify_roles.R` | `path()` staging | WIRED | L34: `path(classify_roles)` input declaration |
| `workflows/hcvtyper.nf` | `bin/classify_roles.R` | call site file() | WIRED | L532: `file("${projectDir}/bin/classify_roles.R")` |

---

### Data-Flow Trace (Level 4)

| Artifact | Data Variable | Source | Produces Real Data | Status |
|----------|---------------|--------|--------------------|--------|
| `summarize.R` (cv_evenness) | `tmp_df$cv_evenness` | `cov$X3` (per-position depth from SAMTOOLS_DEPTH `-aa`) | YES — computed from real depth data `sd(cov$X3)/mean(cov$X3)` with zero guard | FLOWING |
| `summarize.R` (score_candidates) | `candidate_support$dominance_score` | `score_candidates()` over joined `candidate_support` frame (Phase-7 output) | YES — real Phase-7 candidate+support data flows through `join_assembly_support()` → `score_candidates()` | FLOWING |
| `summarize.R` (classify_roles) | `candidate_support$role`, `$role_reason`, `$overall_sample_call` | `classify_roles()` over scored frame | YES — real dominance scores and assembly-support metrics drive classification | FLOWING |
| `summarize.R` (candidates.csv) | `write_csv(candidate_support)` | entire scored+classified long frame | YES — all candidates including background written | FLOWING |

---

### Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| classify_roles.R parses without error | `Rscript -e 'invisible(parse("bin/classify_roles.R")); cat("PARSE_OK\n")'` | `PARSE_OK` | PASS |
| summarize.R parses without error | `Rscript -e 'invisible(parse("bin/summarize.R")); cat("PARSE_OK\n")'` | `PARSE_OK` | PASS |
| is_valid_minor() preserves 3 rules (static) | Inspect L76–88 | Rule 1: allow 1a+1b cross-subtype; Rule 2: block 2k1b vs {1,2,2k1b}; Rule 3: different genotype required | PASS |
| score_weight_evenness > score_weight_reads | `grep score_weight_evenness nextflow.config` | `score_weight_evenness = 3.0`, `score_weight_reads = 1.0` — evenness dominates | PASS |
| n_candidates not re-mapped by new args | Inspect ext.args line in modules_hcv.config | n_candidates at position 8 in ext.args → args[11] in R; score_weight_evenness at position 9 → args[12] | PASS |
| Run unit suite | `bash bin/tests/run_all.sh` | SKIP — tidyverse not installed in host shell; tests require Docker container | SKIP (human needed) |

---

### Probe Execution

No probe scripts declared in PLAN or SUMMARY. Step 7c: SKIPPED (no probes configured for this phase).

---

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|------------|-------------|--------|----------|
| SCORE-01 | 08-01, 08-02 | Dominance score is a combined score over mapped read count, k-mer coverage, and mapping coverage breadth | SATISFIED | `classify_roles.R:114–169` implements `score_candidates()` combining log10(reads), breadth fraction, cv_evenness, and log10(1+kmer_cov); wired in `summarize.R:645–653` |
| SCORE-02 | 08-01, 08-02 | Breadth evenness/uniformity is weighted strongly (read counts alone are contamination-prone) | SATISFIED | `evenness=3.0 > reads=1.0` in both nextflow.config defaults and classify_roles.R defaults; score formula: both `breadth_frac` and `even_fac` carry the `we` weight; test_dominance_score.R Test 4 asserts the headline: 53279-read spiky 4g (cv_evenness~0.20) scores below 4199-read even minor (cv_evenness~0.80) |
| CLASS-01 | 08-01, 08-02 | Each candidate is classified into a strain role: dominant / co-infection / background | SATISFIED | `classify_roles.R:262–312`: every path assigns exactly one of three roles; test_classify_roles.R Test 8 asserts no NA roles and only valid values |
| CLASS-02 | 08-01, 08-02 | A non-dominant candidate is reported as co-infection only if it clears an abundance floor AND has genotype-level assembly support; otherwise background | SATISFIED | `classify_roles.R:277–289`: corroboration check at D-10 ANDed floors; asymmetric refute at D-11; false-4g → background/refuted_denovo; partial ERR1810453 at kmer~5 → co-infection/corroborated (passes reconciled 2.0 floor, would fail old 10.0 floor as test_classify_roles.R Test 2b proves) |
| CLASS-03 | 08-01, 08-02 | Background/artefact candidates are surfaced explicitly with their reason, never silently dropped | SATISFIED | `summarize.R:671–678`: `write_csv(candidate_support)` writes ALL candidates including background; module `emit: candidates` at main.nf:38; publishDir for candidates.csv at modules_hcv.config:343–347; zero-row/NULL guard in classify_roles.R:197–207 returns typed frame never stops |
| CLASS-04 | 08-01, 08-02 | An overall sample call is derived from roles (monoinfection / co-infection / indeterminate) | SATISFIED | `classify_roles.R:315–321`: D-14 logic; `summarize.R:1279,1325`: `overall_sample_call` in pmap_chr trigger list and final select() column; stub Summary.csv header at main.nf:74 includes `overall_sample_call`; `minor_denovo_status` absent from all live code |

---

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| `modules/local/summarize/main.nf` | 82, 85 | Comment words "placeholder" | INFO | These appear inside the `stub:` block only, which is the nf-test stub, not the real `script:` block. The real process calls `summarize.R`. Acceptable. |
| `modules/local/summarize/main.nf` | 86 | Stub candidates.csv header spells `domainance_score` (transposed 'a' and 'i') | WARNING | Typo in stub header only; real `write_csv(candidate_support)` output uses the actual column name `dominance_score` from the R data frame. The stub is for nf-test runs. Does not affect real pipeline output. |

No TBD, FIXME, or XXX markers found in any file modified by this phase.

---

### Human Verification Required

#### 1. Full R Unit Test Suite in Docker

**Test:** Pull the pipeline container and run `bash bin/tests/run_all.sh` from the repository root
```bash
docker run --rm -v $(pwd):/work -w /work \
  community.wave.seqera.io/library/r-gridextra_r-png_r-seqinr_r-tidyverse:3536dd50a17de0ab \
  bash bin/tests/run_all.sh
```
**Expected:** Exit 0 with PASS lines for all tests including:
- `PASS: Test4 (SCORE-02 HEADLINE): spiky high-read 4g scores below a genuine even minor with fewer reads`
- `PASS: Test1 (D-11): false 4g -> background/refuted_denovo; sample monoinfection`
- `PASS: Test2 (CLASS-02): full + partial genuine 2b both -> co-infection/corroborated`
- `PASS: Test4 (D-11): de-novo-failed-for-both genuine minor -> co-infection/uncorroborated_kept`
- `ALL PASS` from both test_dominance_score.R and test_classify_roles.R

**Why human:** tidyverse is not installed in the host environment. All code logic has been verified by static analysis and parse check; this is a runtime execution confirmation that the functions produce correct numeric outputs when called with the evidence-table fixture data. The test code is substantive (not stubs), sources the real helpers, and contains no inline re-implementations — but it cannot be executed without the container.

---

### Gaps Summary

No gaps found. All 13 truths are verified by codebase evidence. All 6 requirements satisfied. All artifacts are substantive and wired. The only item requiring human action is a runtime test suite execution in Docker, which is a normal operational constraint for this containerized pipeline, not a code gap.

**Notable observations (not blockers):**
- The stub `candidates.csv` header in `modules/local/summarize/main.nf:86` has a typo (`domainance_score`). The real output column name in the R frame is `dominance_score` — this only affects nf-test stub runs, not the live pipeline. Recommend fixing in Phase 9 or a patch commit.
- The nf-test snapshot regen (deferred per PLAN.md) is correctly tracked as a deferred item, not a Phase 8 gate.

---

_Verified: 2026-06-14T12:00:00Z_
_Verifier: Claude (gsd-verifier)_
