# Phase 6: Neutral Candidate Selection - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-06-12
**Phase:** 6-neutral-candidate-selection
**Areas discussed:** N-candidate output shape, Candidate ranking & dedup rule, Where 1a/1b + 2k1b exceptions live, Phase 6/9 boundary (filenames + gate logic)

---

## N-candidate output shape

| Option | Description | Selected |
|--------|-------------|----------|
| Long: one row per candidate | sample, candidate_rank, candidate_ref, reads, cov, confirmation_status — N rows/sample; natural channel fan-out, scales to N | ✓ |
| Numbered wide slots | cand1_ref, cand2_ref… single row/sample; ragged for variable N | |
| Keep wide 2-slot, defer N>2 | Minimal change; punts REFSEL-02's parameterizable-N | |

**User's choice:** Long: one row per candidate (Recommended)
**Notes:** Accepted the larger one-time change to hcvtyper.nf routing + summarize.R join in exchange for clean N-scaling and uniform per-candidate mapping.

---

## Candidate ranking & dedup rule

| Option | Description | Selected |
|--------|-------------|----------|
| Rank by reads, one ref per distinct subtype | Top ref per subtype, top N subtypes by reads; preserves 1a/1b | ✓ |
| Rank by reads, one ref per distinct genotype | Always different-genotype candidates; loses 1a/1b | |
| Top-N references, no dedup | Can pick redundant same-subtype refs | |

**User's choice:** Rank by reads, one ref per distinct subtype (Recommended)
**Notes:** Distinct-subtype dedup keeps the 1a/1b co-infection case alive without a special exception at selection. Ranking metric changes from today's coverage-breadth minor tiebreak to read recruitment (REFSEL-01) — an intended behaviour change.

---

## Where 1a/1b + 2k1b exceptions live

| Option | Description | Selected |
|--------|-------------|----------|
| Move to Phase 8 classification | Selection neutral; validity rules become co-infection-role rules in Phase 8 | ✓ |
| Keep at selection in Phase 6 | Behaviour-preserving; criterion #5 holds standalone; bakes interpretation into selection | |
| Split | Structural dedup at selection, validity to Phase 8 | |

**User's choice:** Move to Phase 8 classification (Recommended)
**Notes:** Truest to "interpretation only at summary." Explicitly accepted the consequence that Phase 6 in isolation may map an extra/different second reference (2k1b or same-genotype second subtype), so success-criterion #5 cannot hold standalone — captured in CONTEXT D-05 and flagged for planner/verifier. True reproduction gated at Phase 9 (COMPAT-01).

---

## Phase 6/9 boundary (filenames + gate logic)

| Option | Description | Selected |
|--------|-------------|----------|
| Additive shim in Phase 6, cut over in Phase 9 | Keep _major.fa/_minor.fa + legacy columns (cand_1→major, cand_2→minor); confirmation_status added but minor_call/gate_flag retained; full cutover Phase 9 | ✓ |
| Cut over to cand1/cand2 in Phase 6 | Rename + summarize.R parsing now; pulls COMPAT-02 forward; Phase 6 becomes large | |

**User's choice:** Additive shim in Phase 6, cut over in Phase 9 (Recommended)
**Notes:** Keeps every intermediate phase shippable; the deliberately lockstepped filename-migration + summarize.R parsing + one-release column aliasing stay together in Phase 9.

## Claude's Discretion

- Exact long-format CSV column names/types (subject to legacy-column retention).
- `confirmation_status` vocabulary in Phase 6 (real gating/roles land Phase 8).
- Candidate-count parameter name + validation (default 2, integer ≥ 1), following typed-param conventions.
- Single-candidate / no_mapping / empty-stats behaviour — preserve existing safe defaults, generalize to candidate set.

## Deferred Ideas

- Combined dominance score + dominant/co-infection/background roles → Phase 8 (SCORE/CLASS).
- Re-homed different-genotype / 1a-1b / 2k1b validity rules → Phase 8.
- Filename slot migration + summarize.R parsing cutover + legacy column aliasing → Phase 9 (COMPAT-02/03).
- Per-genotype assembly support + genotype-level join → Phase 7 (ASUP).
- Reviewed-not-folded todos: HCVGLUE-parallel (paused v2.0 Phase 5), Tanoti removal (already shipped).
