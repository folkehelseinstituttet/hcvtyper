---
phase: "09"
slug: compatibility-filename-migration-regression-suite
status: verified
threats_open: 0
asvs_level: 1
created: 2026-06-14
---

# Phase 09 — Security

> Per-phase security contract: threat register, accepted risks, and audit trail.

---

## Trust Boundaries

| Boundary | Description | Data Crossing |
|----------|-------------|---------------|
| samplesheet → pipeline | Operator-controlled sample names flow into filenames; unchanged by this phase | sample metadata (non-sensitive) |
| PARSEFIRSTMAPPING emit → workflow fan-out | The candidate_fasta tuple (meta + per-rank FASTA list) crosses into the Nextflow channel routing; produced upstream in the same pipeline | candidate FASTA paths |
| committed `.snap` golden files → CI | Tool-generated md5 fixtures consumed by nf-test in CI; stale slot strings break the build | md5 hashes / test fixtures |
| candidates CSV → summarize.R | R-emitted candidate_rank/candidate_ref consumed as the join key; produced upstream in the same pipeline | CSV with non-sensitive pipeline metadata |
| candidates CSV / stats files → summarize.R test | Synthetic candidate_rank/candidate_ref + cand-slot stats files staged in tests; same boundary the production pipeline crosses | synthetic fixtures (no PII) |

---

## Threat Register

| Threat ID | Category | Component | Disposition | Mitigation | Status |
|-----------|----------|-----------|-------------|------------|--------|
| T-09-01 | Denial of Service | `candidate_rank.toInteger()` in ext.prefix closures | accept | Rank is guaranteed non-null for a mapped candidate; never coerced in the workflow fan-out (kept String) — existing Pitfall-3 guard preserved. | closed |
| T-09-SC | Tampering | npm/pip/cargo installs (supply chain) | accept | No package installs in this phase; all sub-plans run in the already-pinned `NEXTFLOW` conda env or pinned Seqera container. | closed |
| T-09-02 | Tampering — silent sample loss | `workflows/hcvtyper.nf` fan-out (`remainder: true` + rank String + null-FASTA guard) | mitigate | All four anti-regression invariants preserved verbatim: (1) `remainder: true` on `candidate_fasta` join; (2) rank kept as String, never `.toInteger()`; (3) `filter { entry -> entry[0]['confirmation_status'] == 'pass' && entry[1] != null }`; (4) `assert new_meta.id == new_meta.sample`. Verified in 09-02-SUMMARY.md. | closed |
| T-09-02-SNAP | Tampering — false-green CI | Regenerated `.snap` md5 fixtures | mitigate | Diff-guarded per Pitfall-4: only slot-rename content changed (`*_major.fa`/`*_minor.fa` → `*_cand1.fa`/`*_cand2.fa`); every non-filename md5 byte-identical. Both module nf-tests re-run green without `--update-snapshot`. Verified in 09-02-SUMMARY.md. | closed |
| T-09-03 | Tampering — silent data loss | `Major_reference`/`Minor_reference` join keys (L366, L1043 of summarize.R) | mitigate | Cleaned-ref values kept byte-identical via `_cand[0-9]+$` strip; `cv_by_ref` strip updated in lockstep. Verified in 09-03-SUMMARY.md. | closed |
| T-09-04 | Tampering — silent data loss | `summarize.R` candidate_rank join under the renamed cand-slot | mitigate | COMPAT-02 smoke assertion in `test_compat.R` fails loudly if `.cand1.`/`.cand2.` filenames silently empty `Major_*/Minor_*`. This test is the executable guard for the highest-risk integration site (RESEARCH Pitfall 1). Verified in 09-04-SUMMARY.md. | closed |
| T-09-04-FP | Repudiation — false-green baseline | COMPAT-01 golden values in `compat_golden.csv` | mitigate | Blocking human-verify checkpoint completed before assertions locked: reference accessions `1a_M62321`/`1b_D90208` confirmed real genotype-1 entries; `Major_/Minor_genotype_mapping` tokens confirmed against actual `summarize.R` `separate()` output. Human-approved ("Yes to both") in Task 1 of Plan 04. Verified in 09-04-SUMMARY.md. | closed |

*Status: open · closed*
*Disposition: mitigate (implementation required) · accept (documented risk) · transfer (third-party)*

---

## Accepted Risks Log

| Risk ID | Threat Ref | Rationale | Accepted By | Date |
|---------|------------|-----------|-------------|------|
| AR-09-01 | T-09-01 | `candidate_rank.toInteger()` is only called inside `ext.prefix` closures where rank is guaranteed non-null by the upstream `candidate_fasta` emit contract; the Pitfall-3 guard prevents the DoS scenario at the workflow fan-out layer. Risk is structurally excluded by design. | Plan threat model | 2026-06-14 |
| AR-09-SC | T-09-SC | Phase introduces no new dependencies. All execution environments (NEXTFLOW conda env, Seqera r-seqinr_r-tidyverse container) are pre-pinned in the repo. Supply-chain surface is unchanged from the project baseline. | Plan threat model | 2026-06-14 |

---

## Security Audit Trail

| Audit Date | Threats Total | Closed | Open | Run By |
|------------|---------------|--------|------|--------|
| 2026-06-14 | 7 | 7 | 0 | gsd-secure-phase (short-circuit: register_authored_at_plan_time=true, threats_open=0) |

---

## Sign-Off

- [x] All threats have a disposition (mitigate / accept / transfer)
- [x] Accepted risks documented in Accepted Risks Log
- [x] `threats_open: 0` confirmed
- [x] `status: verified` set in frontmatter

**Approval:** verified 2026-06-14
