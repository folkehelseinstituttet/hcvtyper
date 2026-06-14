---
phase: 09-compatibility-filename-migration-regression-suite
plan: 01
subsystem: targeted-mapping-filename-migration
tags: [COMPAT-02, D-01, D-04, filename-slot, nextflow, r]
requires:
  - "Phase 6 neutral candidate selection (candidate_rank in meta, *.candidates.csv)"
provides:
  - "Uniform .cand{rank}. filename slot across all 6 TARGETED_MAPPING ext.prefix closures"
  - "<sample>.<ref>_cand{k}.fa per-rank FASTAs written for every selected candidate (N, not 2)"
  - "PARSEFIRSTMAPPING N-FASTA candidate_fasta emit (replaces major_mapping/minor_mapping)"
  - "Cand-slot stubs in parsefirstmapping + blastparse (no .major.fa/.minor.fa)"
affects:
  - "workflows/hcvtyper.nf (Plan 02 must rewire the fan-out to candidate_fasta + lift N>2 guard)"
  - "bin/summarize.R (Plan 03 must parse the .cand{rank}. slot in lockstep)"
tech-stack:
  added: []
  patterns:
    - "Slot-from-rank: cand${meta.candidate_rank.toInteger()} (meta1. for SAMTOOLS_DEPTH)"
    - "N-candidate FASTA write loop (seq_along) replacing fixed two-slot writes"
    - "Stub convention: : > cand-slot FASTAs matching the 2-row stub candidates.csv"
key-files:
  created: []
  modified:
    - conf/modules_hcv.config
    - bin/summarize_mapping_to_all_references.R
    - modules/local/parsefirstmapping/main.nf
    - modules/local/blastparse/main.nf
decisions:
  - "blastparse major_fasta/minor_fasta collapsed to single candidate_fasta emit (no downstream consumers confirmed)"
  - "workflows/hcvtyper.nf NOT modified — workflow fan-out rewire + N>2 guard lift is Plan 02 (same-wave lockstep)"
metrics:
  duration: ~3min
  completed: 2026-06-14
  tasks: 3
  files: 4
---

# Phase 9 Plan 01: Targeted-Mapping Filename Migration (producer side) Summary

Migrated the producer side of the `.major.`/`.minor.` → `.cand{rank}.` filename-slot
rename: 6 `ext.prefix` closures, the FASTA-write loop in the selection R script, and the
`parsefirstmapping`/`blastparse` module emits + stubs now all use the uniform `cand{rank}`
slot, and PARSEFIRSTMAPPING exposes an N-FASTA `candidate_fasta` channel for the Plan 02 fan-out.

## What Was Built

### Task 1 — 6 ext.prefix closures → uniform cand{rank} (`conf/modules_hcv.config`)
Dropped the `rank==1?'major':(rank==2?'minor':"cand${rank}")` ternary in all six
TARGETED_MAPPING closures, replacing it with `cand${...candidate_rank.toInteger()}`:
- SAMTOOLS_SORMADUP (`meta.`, `.nodup`)
- SAMTOOLS_DEPTH (`meta1.` namespace preserved, `.nodup`)
- STATS_WITHDUP (`meta.`, `.withdup`)
- STATS_MARKDUP (`meta.`, `.nodup`)
- IVAR_CONSENSUS (`meta.`, no `${meta.reference}` field)
- CONSENSUS_DISTANCE (`meta.`, no `${meta.reference}` field)

`.toInteger()` kept inside the closures (rank guaranteed present for a mapped candidate; V5 note).

### Task 2 — N-candidate FASTA write (`bin/summarize_mapping_to_all_references.R`)
Replaced the fixed two-slot `_major.fa`/`_minor.fa` write with a `for (k in seq_along(selected_refs))`
loop emitting `<sample>.<ref>_cand{k}.fa` for every selected candidate. Preserved the
`if (length(selected_refs) > 0)` guard and the single `read.fasta(file = references)` read so
the no-mapping branch writes nothing and never crashes (D-06 / Pitfall 3). Updated the comment block.

### Task 3 — N-FASTA emit + module stubs (`modules/local/{parsefirstmapping,blastparse}/main.nf`)
- **parsefirstmapping**: collapsed `major_mapping`/`minor_mapping` into a single
  `tuple val(meta), path("*.parsefirstmapping.csv"), path("*_cand*.fa"), emit: candidate_fasta, optional: true`.
  `csv`/`candidates`/`versions` emits unchanged. Stub now writes `${prefix}.cand1.fa` / `${prefix}.cand2.fa`.
- **blastparse**: confirmed no `BLASTPARSE.out.major_fasta`/`.minor_fasta` consumers in `workflows/`
  or `subworkflows/`, so collapsed both fasta emits into one `emit: candidate_fasta` globbing `*cand*.fa`.
  Stub now writes `${prefix}.cand1.fa` / `${prefix}.cand2.fa`; all other stub outputs unchanged.

## Verification

| Check | Result |
|-------|--------|
| 6 `cand${...candidate_rank.toInteger()}` closures (5 meta. + 1 meta1.) | PASS |
| No non-comment `'major'`/`'minor'` slot literal in config | PASS |
| R script writes `_cand{k}.fa` in a loop; no non-comment `_major.fa`/`_minor.fa` | PASS |
| no-mapping guard + single `read.fasta` preserved | PASS |
| parsefirstmapping `emit: candidate_fasta`; no non-comment `major.fa`/`minor.fa` | PASS |
| blastparse no non-comment `major.fa`/`minor.fa`; stub writes cand1/cand2 | PASS |
| No BLASTPARSE fasta-emit consumers (rename safe) | PASS |

## Deviations from Plan

None — plan executed exactly as written.

Note (expected handoff, not a deviation): `workflows/hcvtyper.nf` L407-408 still references the
removed `PARSEFIRSTMAPPING.out.major_mapping`/`minor_mapping` emits. Rewiring the fan-out to the
new `candidate_fasta` emit and lifting the `n_candidates > 2` guard is Plan 02 scope (the file is
not in this plan's `files_modified`). This producer/consumer pair is designed to land in the same
wave; the workflow is transiently inconsistent until Plan 02 commits.

## Known Stubs

None introduced. The module `stub:` blocks are existing -stub-run scaffolding (deterministic dummy
outputs), migrated to cand-slot filenames; not data stubs in the production path.

## Self-Check: PASSED

- FOUND: conf/modules_hcv.config
- FOUND: bin/summarize_mapping_to_all_references.R
- FOUND: modules/local/parsefirstmapping/main.nf
- FOUND: modules/local/blastparse/main.nf
- FOUND commit 84a043f (Task 1)
- FOUND commit 7599299 (Task 2)
- FOUND commit 517a9d5 (Task 3)
