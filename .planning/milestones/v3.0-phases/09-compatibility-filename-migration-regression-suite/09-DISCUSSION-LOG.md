# Phase 9: Compatibility, Filename Migration + Regression Suite - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-06-14
**Phase:** 09-compatibility-filename-migration-regression-suite
**Areas discussed:** Filename slot migration scope, COMPAT-03 aliasing semantics, COMPAT-01 golden baseline test format, Regression suite structure for TEST-01

---

## Filename Slot Migration Scope

### N>2 support

| Option | Description | Selected |
|--------|-------------|----------|
| Full N-slot generalization | Rename .major./.minor. → .cand1./.cand2., generalize parsefirstmapping to emit N ranked FASTAs, lift the n>2 guard. Matches hcvtyper.nf:83 comment. | ✓ |
| 2-slot rename only | Just rename .major.→.cand1. and .minor.→.cand2. in stubs and summarize.R parser. Keep guard, smaller scope. | |

**User's choice:** Full N-slot generalization
**Notes:** hcvtyper.nf comment at L83 explicitly says "N>2 support arrives with Phase 9 filename-slot migration."

### summarize.R parser approach

| Option | Description | Selected |
|--------|-------------|----------|
| Refactor to candidate_rank join | Use candidate_rank from candidates CSV as join key; replaces fragile filename-position extraction. Cleaner, handles N slots naturally. | ✓ |
| Keep filename-position parsing | Keep extracting slot name from position 3, update strings from major/minor to cand1/cand2/etc. Minimal but fragile. | |

**User's choice:** Refactor to candidate_rank join

### Module-level nf-test snapshots

| Option | Description | Selected |
|--------|-------------|----------|
| Update module snapshots in Phase 9 | Regenerate affected .snap files for parsefirstmapping, blastparse stubs. Keeps CI consistent. | ✓ |
| Defer snapshot updates | Leave module-level .snap files unchanged; only update pipeline-level test + R tests. | |

**User's choice:** Update module snapshots in Phase 9

---

## COMPAT-03 Aliasing Semantics

| Option | Description | Selected |
|--------|-------------|----------|
| Validate only — no new columns | Old Major_*/Minor_* columns already present (Phase 8 kept them). Phase 9 adds a regression test confirming their presence. No structural change to summarize.R. | ✓ |
| Re-point Major_*/Minor_* to role-based values | Rename mapping-stat Major_reference → Major_legacy_reference and add new Major_reference = dominant candidate ref from role classification. Bigger change. | |

**User's choice:** Validate only — no new columns

---

## COMPAT-01 Golden Baseline Test Format

| Option | Description | Selected |
|--------|-------------|----------|
| R unit test with CSV fixture | test_compat.R runs summarize.R on a small in-memory fixture; asserts core strain-call columns match expected v1.0/v2.0 values. Fast, no Nextflow needed. | ✓ |
| nf-test pipeline snapshot | Regenerate tests/default.nf.test snapshot after Phase 9; treat snapshot pass as COMPAT-01 proof. Slower, requires full pipeline run. | |

**User's choice:** R unit test with CSV fixture

### Core columns scope

| Option | Description | Selected |
|--------|-------------|----------|
| Core strain-call columns only | Assert: Major_reference, Minor_reference, Major_genotype_mapping, Minor_genotype_mapping, overall_sample_call. New columns not compared. | ✓ |
| All non-new Summary.csv columns | Assert every column from v1.0/v2.0. Stricter but brittle. | |

**User's choice:** Core strain-call columns only

---

## Regression Suite Structure (TEST-01)

| Option | Description | Selected |
|--------|-------------|----------|
| One new test_compat.R file | Single file covering COMPAT-01 (golden baseline), COMPAT-02 (summarize.R cand-slot parsing), COMPAT-03 (legacy columns present), COMPAT-04 (1a/1b + 2k1b exceptions end-to-end). Auto-picked up by run_all.sh. | ✓ |
| Add assertions to existing files | Extend test_classify_roles.R, test_summarize_denovo.R, test_candidate_selection.R with compat assertions. No new file, harder to isolate. | |

**User's choice:** One new test_compat.R file

---

## Claude's Discretion

- Exact new channel emit names for N-FASTA fan-out in parsefirstmapping (nf-core conventions).
- Whether N-FASTA emit is tuple-per-candidate or collected list — minimize change to TARGETED_MAPPING fan-out.
- Exact test_compat.R fixture structure (inline synthetic data or flagoff_golden.csv as v1.0 reference anchor).
- Whether test_compat.R uses subprocess pattern (system2) or function-sourcing — pick lower-friction given candidate_rank refactor.

## Deferred Ideas

- Dropping legacy Major_*/Minor_* columns — explicitly deferred to next milestone post-v3.0.
- Renaming major_mapping/minor_mapping channel identifiers to role-semantic names — optional tidy-up, not Phase 9 hard requirement.
- Cross-platform Singularity/Conda nf-test profiles — existing CI Docker only, out of scope.
