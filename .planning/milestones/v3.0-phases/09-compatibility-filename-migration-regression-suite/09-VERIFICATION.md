---
phase: 09-compatibility-filename-migration-regression-suite
verified: 2026-06-14T20:00:00Z
status: gaps_found
score: 4/5 roadmap success criteria verified (SC-1 partial due to open CR-01)
gaps:
  - truth: "A full run produces correctly-parsed candidate rows AND variation-plot PNGs with the new slot names (SC-1 / COMPAT-02)"
    status: partial
    reason: "CR-01 (from 09-REVIEW.md): summarize.R lines 941-942 still split variation-plot PNGs using grepl(\"major\", ...) / grepl(\"minor\", ...). After the Phase-9 slot rename the PNG filenames carry _nodup (position 3 of the BAM basename) instead of major/minor, so both major_plots and minor_plots are always empty. Variation_plot_major.png and Variation_plot_minor.png are never written. Upstream root cause is plot_bam_variation.R line 24 (position-3 slot parse). This is a silent regression with no error emitted."
    artifacts:
      - path: "bin/summarize.R"
        issue: "Lines 941-942: grepl(\"major\", variation_plot_files) / grepl(\"minor\", ...) never match after the cand-slot rename"
      - path: "bin/plot_bam_variation.R"
        issue: "Line 24: major_minor <- unlist(str_split(basename(bam_file), pattern = \"\\\\.\"))[[1]][3] extracts position 3 (\"nodup\"), not the cand slot"
    missing:
      - "Update summarize.R lines 941-942 to match on cand: major_plots <- variation_plot_files[grepl(\"_cand1\\\\.\", ...)] and minor_plots <- variation_plot_files[grepl(\"_cand2\\\\.\", ...)]"
      - "Update plot_bam_variation.R line 24 to extract the cand-slot token from the correct dot-field position"
  - truth: "The -stub-run of the full workflow exercises the candidate_fasta fan-out path (TEST-01 / COMPAT-02 stub coverage)"
    status: failed
    reason: "CR-02 (from 09-REVIEW.md): PARSEFIRSTMAPPING stub writes ${prefix}.cand1.fa / ${prefix}.cand2.fa (period before cand) but the declared emit glob is *_cand*.fa (underscore before cand). No stub FASTA matches the glob; candidate_fasta is perpetually empty in any -stub-run. Confirmed by parsefirstmapping snapshot line 175: 'candidate_fasta': []. The workflow fan-out join at hcvtyper.nf:400 has zero stub coverage."
    artifacts:
      - path: "modules/local/parsefirstmapping/main.nf"
        issue: "Lines 82-83 stub writes ${prefix}.cand1.fa/${prefix}.cand2.fa; line 27 emit glob *_cand*.fa requires underscore before cand"
      - path: "modules/local/parsefirstmapping/tests/main.nf.test.snap"
        issue: "Stub snapshot at line 175 records 'candidate_fasta': [] confirming the glob mismatch"
    missing:
      - "Rename stub FASTA outputs to ${prefix}.stubref_cand1.fa and ${prefix}.stubref_cand2.fa (or any name matching *_cand*.fa)"
      - "Regenerate the parsefirstmapping stub snapshot after the rename"
---

# Phase 9: Compatibility, Filename Migration + Regression Suite — Verification Report

**Phase Goal:** The redraw ships non-breaking. The output-filename slot migration and its consumer-side parsing change land together, legacy columns are aliased for one release, the preserved HCV exceptions still hold, and the regression suite proves all of it in CI.
**Verified:** 2026-06-14T20:00:00Z
**Status:** gaps_found
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths (Roadmap Success Criteria)

| # | Truth (Roadmap SC) | Status | Evidence |
|---|---|---|---|
| SC-1 | Filename slot migrates .major./.minor. -> .cand1./.cand2. in lockstep with summarize.R parsing | PARTIAL | Producer (6 closures, R loop, module emits) and consumer (candidate_rank join) verified. BLOCKER: variation-plot grid silently suppressed post-rename (CR-01 from 09-REVIEW.md). |
| SC-2 | Legacy Major_*/Minor_* columns aliased alongside new role-based columns | VERIFIED | test_compat.R COMPAT-03 block asserts both column sets co-present; summarize.R emits both. |
| SC-3 | Pipeline reproduces v1.0/v2.0 reporting on golden baseline (COMPAT-01) | VERIFIED | compat_golden.csv exists, human-approved; test_compat.R COMPAT-01 block drives real summarize.R and compares five core strain-call columns for MONO and COINF cases. |
| SC-4 | 1a/1b co-infection exception and 2k/1b recombinant suppression preserved (COMPAT-04) | VERIFIED | test_compat.R COMPAT-04 block (Part A subprocess + Part B helper-level) asserts 1b role == co-infection and 2k1b role == background/recombinant_2k1b. |
| SC-5 | CI regression suite extended and green (TEST-01) | PARTIAL | test_compat.R exists, auto-discovered by run_all.sh test_*.R glob, prints ALL PASS in pinned container per SUMMARY. Caveat: -stub-run coverage of fan-out path broken by CR-02. |

**Score:** 4/5 truths verified (SC-1 and SC-5 partial due to CR-01 and CR-02)

### Required Artifacts

| Artifact | Expected | Status | Details |
|---|---|---|---|
| `conf/modules_hcv.config` | 6 TARGETED_MAPPING ext.prefix closures with uniform cand{rank} slot | VERIFIED | 5 meta. + 1 meta1. closures confirmed at L131, L228, L240, L248, L257, L263; 0 non-comment 'major'/'minor' literals |
| `bin/summarize_mapping_to_all_references.R` | N-candidate FASTA write loop with _cand{k}.fa naming | VERIFIED | Loop present (seq_along), writes paste0(sampleName, ".", ref, "_cand", k, ".fa"); 0 non-comment _major.fa/_minor.fa |
| `modules/local/parsefirstmapping/main.nf` | N-FASTA candidate_fasta emit + cand-slot stub | VERIFIED (emit) / STUB MISMATCH (stub) | emit: candidate_fasta at L27 with glob *_cand*.fa confirmed; stub at L82-83 writes .cand1.fa/.cand2.fa (period) which does NOT match the glob — CR-02 BLOCKER |
| `modules/local/blastparse/main.nf` | Collapsed candidate_fasta emit; stub writes cand-slot FASTAs | VERIFIED | emit: candidate_fasta at L22 with *cand*.fa; stub writes ${prefix}.cand1.fa / ${prefix}.cand2.fa; no legacy major.fa/minor.fa |
| `workflows/hcvtyper.nf` | n_candidates>2 guard deleted; rank-indexed fan-out using candidate_fasta | VERIFIED | Guard absent (grep returns 0); join at L400 uses PARSEFIRSTMAPPING.out.candidate_fasta with remainder: true; FASTA pick by _cand${rank}. at L429; all 4 anti-regression invariants present |
| `modules/local/parsefirstmapping/tests/main.nf.test` | No legacy emit names; asserts candidate_fasta and _cand1.fa/_cand2.fa | VERIFIED | 0 non-comment major_mapping/minor_mapping references; candidate_fasta in all 5 snapshot maps; L117 asserts 2k1b_AB031663_cand1.fa; L121-122 assert _cand1.fa present and _cand2.fa absent for monoinfection |
| `modules/local/parsefirstmapping/tests/main.nf.test.snap` | Regenerated with cand-slot filenames; 0 legacy slot strings | VERIFIED | 0 major.fa/minor.fa/major_fasta/minor_fasta/major_mapping/minor_mapping; real-run entries carry _cand1.fa/_cand2.fa; stub entry records candidate_fasta:[] (CR-02 artifact) |
| `modules/local/blastparse/tests/main.nf.test.snap` | Regenerated with cand-slot filenames; 0 legacy slot strings | VERIFIED | 0 legacy strings; toy.cand1.fa/toy.cand2.fa confirmed; candidate_fasta keys present |
| `bin/summarize.R` | candidate_rank-join parse; 0 first_major_minor == "major"/"minor"; _cand[0-9]+ strip | VERIFIED (join) / BLOCKER (variation plot) | candidate_rank_lookup join at L287-288 and used in all three stats loops + cv_by_ref + consensus-distance block; 0 non-comment first_major_minor == "major"/"minor"; _cand[0-9]+$ strip at L326, L399, L518, L537; BUT lines 941-942 still grepl("major"/"minor") on variation-plot filenames — CR-01 BLOCKER |
| `bin/tests/test_compat.R` | Regression suite covering COMPAT-01/02/03/04; auto-discovered; prints ALL PASS | VERIFIED | File exists (379 lines); entrypoint boilerplate confirmed; COMPAT-01/02/03/04 assertion blocks present; ends with cat("\nALL PASS\n"); run_all.sh uses test_*.R glob |
| `bin/tests/fixtures/compat_golden.csv` | Human-verified two-row golden anchor | VERIFIED | File exists; two rows: MONO (1a_M62321, NA, 1a, NA, monoinfection) and COINF (1a_M62321, 1b_D90208, 1a, 1b, co-infection); human-approved per 09-04-SUMMARY.md |

### Key Link Verification

| From | To | Via | Status | Details |
|---|---|---|---|---|
| conf/modules_hcv.config closures | meta.candidate_rank | cand${meta.candidate_rank.toInteger()} | WIRED | 6 closures confirmed; 5 using meta., 1 using meta1. |
| bin/summarize_mapping_to_all_references.R | modules/local/parsefirstmapping/main.nf | *_cand*.fa FASTA glob matches _cand{k}.fa | WIRED | R writes <sample>.<ref>_cand{k}.fa; module glob *_cand*.fa matches (underscore present in real output) |
| workflows/hcvtyper.nf fan-out | PARSEFIRSTMAPPING.out.candidate_fasta | .join(candidate_fasta, remainder: true) + _cand${rank}. pick | WIRED | L400 join confirmed; L429 FASTA pick by basename contains "_cand${rank}." |
| bin/summarize.R stats loops | candidate_rank_lookup (candidates_long) | left_join on (sampleName, candidate_ref) after str_remove _cand[0-9]+$ | WIRED | Confirmed at L326-327, L399-400, L518-519, L537 |
| bin/tests/test_compat.R | bin/summarize.R | system2("Rscript", ...) subprocess on staged tempdir | WIRED | run_summarize() harness confirmed |
| bin/tests/run_all.sh | bin/tests/test_compat.R | test_*.R glob auto-discovery | WIRED | run_all.sh uses for t in "$here"/test_*.R glob |
| bin/summarize.R lines 941-942 | variation-plot PNG files | grepl("major"/"minor", ...) | BROKEN (CR-01) | After cand-slot rename, no PNG matches grepl("major") or grepl("minor"); both grids always empty |
| parsefirstmapping stub | candidate_fasta emit glob *_cand*.fa | : > ${prefix}.cand1.fa / .cand2.fa | BROKEN (CR-02) | Period before cand does not match underscore-requiring glob; candidate_fasta always empty in stub runs |

### Data-Flow Trace (Level 4)

| Artifact | Data Variable | Source | Produces Real Data | Status |
|---|---|---|---|---|
| bin/summarize.R Major_reference | candidate_rank (from join) | candidate_rank_lookup left_join on candidate_ref | Yes — DB-backed (candidates CSV from pipeline) | FLOWING |
| bin/summarize.R variation_plot_files | major_plots / minor_plots | grepl("major"/"minor") filter on file list | No — filter always returns empty list after rename | STATIC (CR-01) |
| parsefirstmapping stub candidate_fasta | candidate_fasta channel | *.cand1.fa/*.cand2.fa stub files | No — glob never matches stub filenames | DISCONNECTED (CR-02) |

### Behavioral Spot-Checks

Step 7b: SKIPPED for most checks (test suite requires Docker; no server entry point). Existence checks run instead.

| Behavior | Command | Result | Status |
|---|---|---|---|
| 6 cand closures in config | grep -c 'cand\${meta' conf/modules_hcv.config | 6 | PASS |
| 0 legacy major/minor literals in config | grep -v comment + grep -c 'major'/'minor' | 0 | PASS |
| candidate_fasta emit in parsefirstmapping | grep -n 'emit: candidate_fasta' | L27 confirmed | PASS |
| n_candidates>2 guard absent | grep -c 'params.n_candidates > 2' workflows/hcvtyper.nf | 0 | PASS |
| All 4 anti-regression invariants in hcvtyper.nf | grep remainder:true, candidate_rank.toString(), confirmation_status==pass&&entry[1]!=null, assert new_meta.id==new_meta.sample | All 4 found at L400, L425, L442, L417 | PASS |
| candidate_rank join in summarize.R | grep -c 'first_major_minor == "major"' (non-comment) | 0 | PASS |
| _cand[0-9]+ strip in summarize.R | grep -n '_cand\[0-9\]+\$' | L326, L399, L518, L537 | PASS |
| grepl("major"/"minor") for variation plots | grep -n 'grepl.*major.*grepl.*minor' summarize.R:941-942 | Lines 941-942 found — UNFIXED | FAIL (CR-01) |
| stub FASTA glob match | stub writes .cand1.fa, glob requires _cand | Mismatch confirmed via snapshot line 175: candidate_fasta:[] | FAIL (CR-02) |
| test_compat.R exists and substantive | wc -l | 379 lines | PASS |
| compat_golden.csv exists with two rows | cat | MONO + COINF rows confirmed | PASS |
| ALL PASS tail in test_compat.R | tail -3 | cat("\nALL PASS\n") confirmed | PASS |

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|---|---|---|---|---|
| COMPAT-01 | 09-04 | v1.0/v2.0 golden reproduction | SATISFIED | compat_golden.csv human-verified; test_compat.R COMPAT-01 block drives real summarize.R; human checkpoint approved |
| COMPAT-02 | 09-01/02/03 | Lockstep .major./.minor. -> .cand{rank}. rename with summarize.R parse | PARTIAL | Producer/consumer join verified; variation-plot grid broken (CR-01 open); stub fan-out path unexercised (CR-02 open) |
| COMPAT-03 | 09-04 | Legacy Major_*/Minor_* columns aliased alongside new role columns | SATISFIED | test_compat.R COMPAT-03 asserts both column sets present; summarize.R confirmed to emit both |
| COMPAT-04 | 09-04 | 1a/1b exception preserved; 2k/1b suppression preserved | SATISFIED | test_compat.R COMPAT-04 blocks (subprocess + helper-level) confirmed; classify_roles.R wired via source() |
| TEST-01 | 09-04 | Regression suite in CI; auto-discovered; no CI YAML edit | PARTIAL | test_compat.R auto-discovered by run_all.sh test_*.R glob; no CI YAML edit needed; stub-run fan-out coverage absent (CR-02) |

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|---|---|---|---|---|
| bin/summarize.R | 941-942 | grepl("major", ...) / grepl("minor", ...) on variation_plot_files — unfixed legacy slot filter | BLOCKER | Variation-plot PNGs (Variation_plot_major.png, Variation_plot_minor.png) never written; silent regression with no error |
| bin/plot_bam_variation.R | 24 | major_minor <- unlist(str_split(basename(bam_file), "\\."))[3] — position-3 parse still reads "nodup" instead of cand-slot | BLOCKER | Upstream root cause of CR-01; PNG filename contains _nodup not _cand{rank} |
| modules/local/parsefirstmapping/main.nf | 82-83 | Stub writes ${prefix}.cand1.fa / ${prefix}.cand2.fa (period) while emit glob is *_cand*.fa (underscore) | BLOCKER | candidate_fasta always empty in -stub-run; workflow fan-out has no stub coverage |
| bin/tests/test_compat.R | 143-144 | gate_flag = "pass" (real script emits "ok"); minor_call = "co-infection"/"monoinfection" (real emits "yes"/"no") | WARNING | Spurious review_flag trigger on every test case; schema mismatch on minor_call — future assertions on these columns will see wrong values |

Note: WR-02 (1:length(id_files) loop crash when id/ is empty) and WR-03 (doubled _cand{rank} slot in intermediate filenames) are pre-existing or benign-in-production warnings identified in 09-REVIEW.md. They do not block the phase goal and are not introduced by Phase 9.

### Human Verification Required

None — all assertions are programmatically verifiable or deferred to the gap closure items above.

### Gaps Summary

Two blockers were identified in the post-execution code review (09-REVIEW.md, findings CR-01 and CR-02) and confirmed in the codebase:

**CR-01 (BLOCKER): Variation-plot grid silently suppressed.** `bin/summarize.R` lines 941-942 use `grepl("major", variation_plot_files)` and `grepl("minor", ...)` to split variation-plot PNGs into major/minor grids. After the Phase-9 cand-slot rename, the PNG filenames no longer contain "major" or "minor" (they contain "_nodup" at position 3 of the BAM basename, as `plot_bam_variation.R` line 24 has not been updated). Both `major_plots` and `minor_plots` lists are always empty, and neither grid PNG is written. No error is emitted. The fix requires two coordinated changes: update `plot_bam_variation.R` line 24 to extract the cand-slot token, and update `summarize.R` lines 941-942 to match on `_cand1\\.` / `_cand2\\.`.

**CR-02 (BLOCKER): Parsefirstmapping stub FASTA glob mismatch.** The `candidate_fasta` output is declared with glob `*_cand*.fa` (underscore before "cand") at line 27, but the stub block at lines 82-83 creates `${prefix}.cand1.fa` and `${prefix}.cand2.fa` (period before "cand"). No stub file matches the glob. The `candidate_fasta` channel is always empty in `-stub-run` mode, confirmed by the snapshot at line 175 (`"candidate_fasta": []`). This means the workflow fan-out at `hcvtyper.nf:400` has no stub-run coverage — a full-workflow `-stub-run` cannot exercise the fan-out path that was the primary structural change of Phase 9. The fix is to rename the stub files to match the glob (e.g. `${prefix}.stubref_cand1.fa` / `${prefix}.stubref_cand2.fa`) and regenerate the snapshot.

These two gaps both affect the COMPAT-02 / TEST-01 requirements. The variation-plot suppression (CR-01) is the more operationally serious gap as it affects every production run output silently. The stub mismatch (CR-02) affects CI/stub-run coverage and review confidence in the fan-out path.

**One warning is also noted (WR-01):** `test_compat.R` line 144 passes `gate_flag = "pass"` (should be `"ok"`) and line 143 passes `minor_call = "co-infection"/"monoinfection"` (should be `"yes"/"no"`) in its synthetic `parsefirstmapping.csv` fixture. This causes every test case to spuriously trigger the review_flag path in summarize.R. Current assertions do not check review_flag so tests still pass, but the fixture schema mismatch will mislead any future assertion on this field.

---

_Verified: 2026-06-14T20:00:00Z_
_Verifier: Claude (gsd-verifier)_
