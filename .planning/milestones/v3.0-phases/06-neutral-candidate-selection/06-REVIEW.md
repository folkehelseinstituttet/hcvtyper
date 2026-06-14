---
phase: 06-neutral-candidate-selection
reviewed: 2026-06-12T00:00:00Z
depth: standard
files_reviewed: 7
files_reviewed_list:
  - workflows/hcvtyper.nf
  - conf/modules_hcv.config
  - modules/local/parsefirstmapping/main.nf
  - bin/summarize_mapping_to_all_references.R
  - bin/summarize.R
  - bin/tests/test_candidate_selection.R
  - modules/local/parsefirstmapping/tests/main.nf.test
findings:
  critical: 1
  warning: 5
  info: 3
  total: 9
status: issues_found
---

# Phase 6: Code Review Report

**Reviewed:** 2026-06-12
**Depth:** standard
**Files Reviewed:** 7
**Status:** issues_found

## Summary

This phase collapses the asymmetric MAJOR_MAPPING / MINOR_MAPPING alias routes into a single
per-candidate `TARGETED_MAPPING` fan-out over a long-format candidates CSV. At the **default
`n_candidates = 2`** the redesign is internally consistent: the splitCsv fan-out, the
whole-meta join inside `TARGETED_MAPPING`, the `candidate_rank → major/minor` slot derivation in
`conf/modules_hcv.config`, and the `summarize.R` field-3 parse all line up, and the per-rank FASTA
selection (`rank 1 → _major.fa`, `else → _minor.fa`) matches the references the R script actually
writes. The channel-arity lockstep on the HCVGLUE / SUMMARIZE feeds is correct — each former
`.mix(MAJOR…, MINOR…)` pair was correctly collapsed to a single `TARGETED_MAPPING.out.*` collect,
so stats are not halved.

The serious defect is that the design **breaks as soon as the user sets `n_candidates ≥ 3`** — a
configuration the pipeline explicitly exposes as a typed parameter (`conf/modules_hcv.config:21`).
A rank-3 candidate that passes its thresholds is mapped against the **rank-2 reference FASTA**,
silently producing a biologically wrong strain call. The remaining findings are tie-break
non-determinism in candidate ranking, the documented-but-unguarded `n_candidates ≥ 3` slot/GLUE
gaps, and pre-existing latent bugs in `summarize.R` that this phase did not introduce but which sit
in reviewed files.

## Critical Issues

### CR-01: `n_candidates ≥ 3` maps rank-3+ candidate reads against the WRONG reference FASTA

> **RESOLVED 2026-06-12 (commit `141cdc7`).** Guarded rather than fully fixed: the
> N-slot FASTA/filename migration remains Phase 9 (COMPAT-02) scope, but `n_candidates > 2`
> now fails loudly instead of emitting wrong-reference data. Two layers: `nextflow_schema.json`
> sets `maximum: 2` (nf-schema rejects `--n_candidates 3` at validation — verified via
> `-preview`), and `workflows/hcvtyper.nf` has a workflow-entry `error` guard with a
> Phase-9-pointer message as backstop. The silent-wrong-results footgun is closed.

**File:** `workflows/hcvtyper.nf:418-423`, `bin/summarize_mapping_to_all_references.R:210-220`

**Issue:** `n_candidates` is a user-settable typed parameter (default 2, `conf/modules_hcv.config:21`),
and the selection script emits one candidate row per distinct subtype up to `head(n = n_candidates)`.
But the per-rank FASTA writer in `summarize_mapping_to_all_references.R` only ever writes two files:
`<sample>.<ref>_major.fa` for `selected_refs[1]` and `<sample>.<ref>_minor.fa` for `selected_refs[2]`.
There is no `_cand3.fa` / rank-3 FASTA.

In the workflow fan-out the FASTA is chosen with a binary rule:

```groovy
def fasta = (rank == '1') ? major_fasta : minor_fasta
```

So for a rank-3 candidate whose `confirmation_status == 'pass'` (which passes the
`filter { … 'pass' && entry[1] != null }` guard because `minor_fasta` is non-null), the workflow
hands `minor_fasta` — the **rank-2 reference** — to `TARGETED_MAPPING` while `meta.candidate_ref`
holds the rank-3 reference. The rank-3 candidate's reads are then aligned to the rank-2 genome.
`TARGETED_MAPPING` derives `meta.reference` from the rank-2 FASTA basename, so the emitted stats,
consensus, idxstats and GLUE inputs are all labelled/aligned against the wrong reference. This is a
silent data-correctness failure (wrong genotype evidence), not a crash. It directly contradicts the
project core value that a reported strain must be backed by its own orthogonal evidence.

**Fix:** Make the FASTA selection rank-general and fail loudly on a missing rank FASTA. Either (a)
have the R script write a per-rank FASTA for every selected candidate (e.g.
`<sample>.<ref>_cand<rank>.fa`) and select it by rank in the workflow, or (b) explicitly clamp the
fan-out to the two ranks the shim supports and reject `n_candidates > 2` until the cand-k path is
wired. Minimal guard in the flatMap:

```groovy
def rank = new_meta.candidate_rank.toString()
def fasta = (rank == '1') ? major_fasta : (rank == '2' ? minor_fasta : null)
assert fasta != null : "No per-rank FASTA for candidate_rank=${rank} (n_candidates>2 not yet wired); refusing to cross-map ${new_meta.candidate_ref}"
```

(The existing `entry[1] != null` filter would then drop the rank-3 candidate instead of cross-mapping
it — still incomplete, but no longer silently wrong.)

## Warnings

### WR-01: Candidate ranking is non-deterministic on tied read counts

**File:** `bin/summarize_mapping_to_all_references.R:139-143`

**Issue:** `ranked <- subtype_reads %>% … arrange(desc(reads)) %>% head(n = n_candidates) %>%
mutate(candidate_rank = row_number())`. When two subtypes have equal total `reads`, `arrange()`
leaves their relative order to the input frame (a stable but data-incidental order). The candidate
that becomes `candidate_rank == 1` (the "major") versus `2` (the "minor") is then decided by an
arbitrary tie order. Because `candidate_rank` drives the `.major.` / `.minor.` filename slot, the
shim's major/minor assignment — and which strain is gated as the dominant one — can flip between
otherwise-identical inputs. For a co-infection benchmark this is a reproducibility hazard.

**Fix:** Add a deterministic secondary (and tertiary) sort key, e.g.
`arrange(desc(reads), desc(percent_gt_4_int), candidate_ref)` before `head()`/`row_number()`, so a
tie resolves the same way every run.

### WR-02: Per-candidate `confirmation_status` does not encode the major-pass dependency

**File:** `bin/summarize_mapping_to_all_references.R:148-160`, `workflows/hcvtyper.nf:426-434`

**Issue:** `confirmation_status` is computed independently per candidate
(`reads > minRead & percent_gt_4_int > minCov`). The workflow then maps every candidate whose status
is `'pass'`. A rank-2 candidate can therefore be `'pass'` and get fully mapped even when the rank-1
major is `below_threshold`. The inline comment at `workflows/hcvtyper.nf:428-431` asserts "rank 1 ==
the major gate", but the gate is not actually conditional on the major passing at the mapping stage —
the "never report a minor on a failed major" rule is enforced only later in `summarize.R`
(`gate_flag`, `minor_typable`, and the GATE-03 secondary gate). The reporting outcome is likely
preserved, but the comment overstates the guarantee and the pipeline does wasted mapping work for
minors that will be suppressed. If a future edit trusts the comment and drops the downstream gate,
the core invariant breaks.

**Fix:** Either gate the workflow filter on the major's pass state (carry a `major_pass` flag in each
candidate row and require it for rank ≥ 2), or soften the comment to state plainly that the
major-dependency is enforced only in `summarize.R`, not at the fan-out filter.

### WR-03: `n_candidates ≥ 3` candidates are silently dropped from GLUE and from the summary

**File:** `conf/modules_hcv.config:130-131`, `bin/GLUE_json_parser.R:9`, `bin/summarize.R:226`

**Issue:** The config comment (`conf/modules_hcv.config:128-129`) acknowledges that `summarize.R`
only consumes the `major`/`minor` slots, but the consequence is broader than noted: for rank ≥ 3 the
slot becomes `.cand{k}.`, and (a) `GLUE_json_parser.R` matches only `paste0(major_minor,
".nodup.json$")` for literal `major`/`minor`, so rank-3 GLUE JSONs are never parsed; (b)
`summarize.R`'s `first_major_minor <- str_split(...)[[3]]` yields `cand3`, which falls through every
`case_when(first_major_minor == "major"/"minor" ~ …)` to NA and is then collapsed away by the
per-sample `slice(1)`. So even setting aside CR-01, a rank-3 candidate produces orphaned, unreported
output files. Combined with CR-01 this means `n_candidates > 2` is not safely usable.

**Fix:** Hard-fail or warn at workflow start if `n_candidates > 2` until the `cand{k}` slot is wired
through GLUE and `summarize.R` (the config already flags this as Open Q1 — make it enforced, not just
documented).

### WR-04: `id_files` empty-guard and `1:length()` loop can crash SUMMARIZE (pre-existing)

**File:** `bin/summarize.R:543-549`

**Issue:** `if (length(id_files > 0))` evaluates `id_files > 0` (an element-wise comparison over a
character vector, producing a logical vector of the same length) and then takes its `length`, instead
of `length(id_files) > 0`. It happens to be falsy only when `id_files` is empty, so the `id_df`
creation is skipped in exactly the empty case — but then `for (i in 1:length(id_files))` at line 549
runs `1:0 == c(1, 0)`, dereferences `id_files[1]` (NA) and an undefined `id_df`, and errors. This is a
classic `1:length()` antipattern. It is **pre-existing** (not introduced by this phase) but lives in a
reviewed file and would abort the whole summary if the `id/` staging dir is ever empty.

**Fix:** `if (length(id_files) > 0) { … }` and wrap the loop in the same guard, or use
`for (i in seq_along(id_files))`.

### WR-05: `kraken_df` / stats loops also use `1:length(...)` (pre-existing)

**File:** `bin/summarize.R:139, 168, 217, 281, 348`

**Issue:** Same `for (i in 1:length(x))` antipattern as WR-04 across the kraken, first-mapping,
withdup-stats, nodup-stats and coverage loops. Each crashes (or iterates a spurious `i = 0`) if its
corresponding staging dir is empty. In normal runs these dirs are non-empty, but a skip/edge-case run
(e.g. all samples filtered before mapping) can hit them. Pre-existing; called out because this phase
changes what flows into these dirs (single fan-out instead of two alias routes).

**Fix:** Replace each `1:length(...)` with `seq_along(...)`.

## Info

### IN-01: Stub `versions.yml` indentation differs from real script output

**File:** `modules/local/parsefirstmapping/main.nf:49-54` vs `86-89`

**Issue:** The real script emits `versions.yml` with two-space-indented keys via a `<<-` heredoc;
the stub emits them via `echo '  r-base: stub'`. The values intentionally differ (`stub`), but any
snapshot test that captures `versions` from a real run and a stub run will see different content. The
in-file comment already explains the choice; flagged only so reviewers of the snapshot regen are aware
the two paths are deliberately divergent.

**Fix:** None required; ensure the nf-test snapshots for the real and stub runs are kept separate
(they already are).

### IN-02: `ch_denovo` mix-then-collect ordering is non-deterministic

**File:** `workflows/hcvtyper.nf:486`

**Issue:** `BLASTPARSE.out.csv.collect({it[1]}).mix(BLASTPARSE.out.blast_res.collect({it[1]})).collect()`
stages both file sets into one channel for SUMMARIZE. The `mix` order between the two collected lists
is not guaranteed. SUMMARIZE consumes them by `list.files(pattern = ...)` (filename-globbed, not
order-dependent), so this is harmless today, but it is fragile if a future consumer relies on channel
order.

**Fix:** No action needed given the glob-based consumer; consider a comment noting order-independence.

### IN-03: Magic gate constants embedded in defaults

**File:** `bin/summarize_mapping_to_all_references.R:32`, `bin/summarize.R:25-29`

**Issue:** Fallback defaults (`n_candidates <- 2L`, `denovo_min_contig_length <- 1000`,
`denovo_min_kmer_cov <- 2.0`, `denovo_min_blast_identity <- 90`) are duplicated as literals in both R
scripts and again in `conf/modules_hcv.config`. Drift between the config defaults and the R fallbacks
would be silent. Not a bug today (config always passes the args), but a maintainability risk.

**Fix:** Keep the config the single source of truth and treat the R fallbacks purely as crash-guards
(as the comments already intend); optionally assert the arg is present in production rather than
silently defaulting.

---

_Reviewed: 2026-06-12_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
