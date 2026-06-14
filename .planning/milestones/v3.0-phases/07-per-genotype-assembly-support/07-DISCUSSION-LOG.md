# Phase 7: Per-Genotype Assembly Support - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-06-12
**Phase:** 07-per-genotype-assembly-support
**Areas discussed:** Signal vs verdict split, Table grain, Best-contig definition, Additive vs cutover

---

## Signal vs Verdict Split

| Option | Description | Selected |
|--------|-------------|----------|
| Raw metrics + status string | Raw metrics only (best contig length, BLAST % identity, BLAST aln length, k-mer cov); threshold verdict → Phase 8 (CLASS-02); no-support = explicit `assembly_support = "none"` + NA metrics, no row loss | ✓ |
| Metrics + substantiality bool | Phase 7 also applies the threshold floors and emits a per-candidate `substantial` boolean | |
| NA-only, no status string | Raw metrics only; no-support inferred from NA across all columns, no explicit status field | |

**User's choice:** Raw metrics + status string (recommended)
**Notes:** Mirrors Phase 6 D-05 — mechanical signal in Phase 7, interpretation (calibration floors) in Phase 8 reusing existing `denovo_min_*` params.

---

## Table Grain

| Option | Description | Selected |
|--------|-------------|----------|
| Subtype-grain, collapse at join | `blast_parse.R` emits one row per subtype (carrying a genotype key via `genotype_from_subtype`); `summarize.R` collapses to `denovo_match_level` at join time | ✓ |
| Genotype-grain, pre-aggregated | `blast_parse.R` collapses to one row per genotype directly; `denovo_match_level='subtype'` could not be honored from this table | |

**User's choice:** Subtype-grain, collapse at join (recommended)
**Notes:** Keeps `denovo_match_level = "subtype"` usable without recomputation; criterion #1's "per genotype" satisfied at default match level by roll-up.

---

## Best-Contig Definition

| Option | Description | Selected |
|--------|-------------|----------|
| Coherent row, best by length | Pick single best contig by `sc_length`; report that contig's identity, aln length, k-mer cov — one coherent contig backs all four numbers | ✓ |
| Coherent row, best by bitscore | Pick best contig by BLAST bitscore; report its four metrics coherently | |
| Per-metric maxima | Independent max of each metric across the group — numbers may come from different contigs | |

**User's choice:** Coherent row, best by length (recommended)
**Notes:** Matches "substantial contig" framing (ERR1810453's 2,949 bp 2b); most direct route to reproducing today's length-based `minor_contig_length` (criterion #4).

---

## Additive vs Cutover

| Option | Description | Selected |
|--------|-------------|----------|
| Additive shim | Add neutral per-genotype table + candidate join; keep `classify_minor_denovo` / `minor_denovo_status` / `review_flag` working unchanged; legacy path retired in Phase 8 | ✓ |
| Full cutover now | Replace `classify_minor_denovo` and rewire `review_flag` immediately in Phase 7 | |

**User's choice:** Additive shim (recommended)
**Notes:** Mirrors Phase 6 D-06 — keeps Phase 7 shippable and non-breaking; full cutover would change reporting before Phase 8's role logic exists.

---

## Claude's Discretion

- Join site: `bin/summarize.R`, reusing the existing `df_denovo` / `df_blast_out` left-join plumbing (keyed by `sampleName`).
- Candidate-side genotype key: derive from `candidate_ref` via `genotype_from_subtype()` — identical helper on both join sides.
- Exact column names/types of the support table and per-candidate support columns; whether the support table is a new `BLASTPARSE` CSV or folded into existing outputs.
- Preserve typed-empty-tibble / NA-fill behavior on skip-assembly runs (T-03-01 guard).

## Deferred Ideas

- Substantiality threshold as the corroboration verdict + "guilty until corroborated" + role assignment — Phase 8 (CLASS-02/SCORE).
- Retiring the legacy minor-coupled confirmation path once the new join feeds classification — Phase 8.
- Filename slot migration `.major.`/`.minor.` → `.cand1.`/`.cand2.` + `summarize.R` parsing cutover + legacy column aliasing — Phase 9 (COMPAT-02/03).
