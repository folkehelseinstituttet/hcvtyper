# Phase 6: Neutral Candidate Selection - Context

**Gathered:** 2026-06-12
**Status:** Ready for planning

<domain>
## Phase Boundary

Make reference selection during the run **dominance-neutral**: after first-pass mapping, emit a ranked set of N candidate references (`cand_1..cand_n`, default 2) ordered by read recruitment, map each candidate through **one uniform path**, and replace the `gate_flag`/`minor_call` plumbing with per-candidate `confirmation_status` — with no "major"/"minor" semantics baked into the run.

**Covers requirements:** REFSEL-01 (neutral ranked candidates), REFSEL-02 (candidate count parameter, default 2), REFSEL-03 (uniform per-candidate mapping + `confirmation_status`).

**Explicitly NOT in this phase:**
- The combined dominance score and dominant/co-infection/background role assignment → **Phase 8** (SCORE/CLASS).
- The co-infection-validity rules (different-genotype requirement, 1a/1b allow, 2k1b block) → **moved to Phase 8** (see D-05).
- The output-filename slot migration `.major.`/`.minor.` → `.cand1.`/`.cand2.` and `bin/summarize.R` parsing change, plus legacy column aliasing → **Phase 9** (COMPAT-02/03), kept lockstepped there.

</domain>

<decisions>
## Implementation Decisions

### Candidate output shape (N-candidate data contract)
- **D-01:** The selection step emits a **long-format** candidate table — **one row per candidate** (e.g. `sample, candidate_rank, candidate_ref, candidate_reads, candidate_cov, confirmation_status`), N rows per sample. Chosen over numbered wide slots (`cand1_ref`, `cand2_ref`…) and over keeping the fixed 2-slot wide row.
- **Why:** A per-candidate row makes REFSEL-03's "map each candidate through one uniform path" a natural Nextflow channel fan-out (one channel element per row) and scales to any N without ragged columns. Accepted tradeoff: a larger one-time change to `workflows/hcvtyper.nf` routing and to `bin/summarize.R`'s consumption of the selection CSV.

### Candidate ranking & de-duplication
- **D-02:** Candidates are chosen by **read recruitment, one reference per distinct subtype**: pick the top reference (most reads) within each subtype, then take the top N subtypes by total reads. (Not genotype-level dedup, not raw top-N references.)
- **Why:** Distinct-subtype dedup preserves the 1a/1b co-infection case (both genotype 1, different subtypes can both be candidates) which genotype-level dedup would collapse — the exact case the current code special-cases. Raw top-N could pick two near-identical refs of one subtype (redundant mapping).
- **D-03:** Ranking is by **read recruitment** (REFSEL-01), replacing today's asymmetric metric where `major` is by reads but the `minor` is chosen by **coverage breadth** among valid different-genotype refs. This is an intended behaviour change — the second candidate's mapping target may differ from today's minor.

### Uniform mapping path
- **D-04:** Replace the two aliased subworkflow calls `MAJOR_MAPPING`/`MINOR_MAPPING` in `workflows/hcvtyper.nf` with **one uniform per-candidate mapping fan-out** over the long-format candidate channel — every candidate gets identical mapping/stats treatment. `confirmation_status` is carried per-candidate through channels into the summary inputs without index misalignment (use explicit metadata-key joins, per the existing `TARGETED_MAPPING` join discipline).

### Where co-infection-validity rules live
- **D-05:** The **different-genotype requirement, the 1a/1b allow rule, and the 2k1b block** (today inside `is_valid_minor()` in `bin/summarize_mapping_to_all_references.R`) **move OUT of selection and INTO Phase 8 classification.** Phase 6 selection picks top-N distinct subtypes by reads with **no genotype/recombinant filtering**; in Phase 8 a same-genotype-as-dominant candidate or a 2k1b pair resolves to `background`.
- **Why:** Truest realization of the milestone thesis — selection is mechanical, co-infection validity is interpretation and belongs at the summary. COMPAT-04 (exceptions preserved) is satisfied by Phase 8 reproducing the same allow/block outcomes, not by keeping them in selection.
- **CONSEQUENCE (must be honored by planner + verifier):** Because of D-05 + D-03, **Phase 6 in isolation can select/map a different or extra second reference than today** (e.g. a same-genotype second subtype, or a 2k1b pair, that `is_valid_minor()` currently excludes). Therefore ROADMAP Phase 6 **success-criterion #5** ("reproduce the current two-slot selection on the regression fixtures") **cannot hold standalone** and must NOT be treated as a Phase-6 blocker. True non-breaking reproduction is only meaningful once Phase 8 re-suppresses those candidates as `background`, and the golden-baseline reproduction is formally gated at **Phase 9 / COMPAT-01**. Interpret criterion #5 as "reproduces at default settings across Phase 6+8, verified at Phase 9."

### Phase 6 / Phase 9 compatibility boundary
- **D-06:** Phase 6 is **additive** with a legacy compatibility shim. It builds neutral candidate selection internally but **keeps the legacy on-disk contract** so nothing downstream breaks before Phase 9:
  - Targeted-mapping outputs keep the existing `_major.fa` / `_minor.fa` filename tags (`cand_1` → major slot, `cand_2` → minor slot).
  - The legacy `major_*` / `minor_*` / `minor_call` / `gate_flag` CSV columns that `bin/summarize.R` reads are retained (emit alongside / reconstructable from the new long-format candidate table).
  - `confirmation_status` is introduced as a per-candidate field, but legacy `minor_call`/`gate_flag` stay until Phase 9.
- **Why:** Keeps every intermediate phase shippable and runnable. The deliberately lockstepped filename-migration + `summarize.R` parsing change + one-release column aliasing all stay in Phase 9 (COMPAT-02/03) as the roadmap intended.
- **Note on default N=2:** With candidate count = 2 and the shim, `cand_1`/`cand_2` map onto the existing two mapping slots, so routing topology is unchanged at the default.

### Claude's Discretion
- Exact column names/types of the long-format candidate CSV (subject to the legacy-column-retention constraint in D-06) — planner/researcher to specify.
- The neutral value `confirmation_status` carries in Phase 6 (e.g. a provisional `selected` / threshold outcome) given the real gating + role logic lands in Phase 8. Keep current major-pass threshold behaviour observable so nothing regresses, but the field's vocabulary is open.
- The candidate-count parameter name and validation (e.g. `params.n_candidates` / `params.max_candidates`, default 2, integer ≥ 1). Follow existing typed-param conventions (`conf/modules_hcv.config`, schema).
- Single-candidate and `no_mapping` / empty-stats behaviour — preserve the existing safe defaults (`gate_flag = "no_mapping"`, no crash on length-zero frames) while generalizing to the candidate set.

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Phase scope & requirements
- `.planning/ROADMAP.md` §"Phase 6: Neutral Candidate Selection" — goal, success criteria (note criterion #5 caveat in D-05), dependencies.
- `.planning/REQUIREMENTS.md` — REFSEL-01/02/03 (Phase 6), plus ASUP/CLASS/SCORE (Phases 7/8) and COMPAT (Phase 9) for cross-phase awareness.
- `.planning/PROJECT.md` — Core Value, constraints (non-breaking, genotype-level match, preserved exceptions), milestone framing.

### Code touchpoints (the files this phase changes)
- `bin/summarize_mapping_to_all_references.R` — the selection script; today produces the 10-column `major_*/minor_*/minor_call/gate_flag` row and writes `_major.fa`/`_minor.fa`. Primary site for neutral candidate ranking (D-01..D-03) and the validity-rule removal (D-05).
- `bin/genotype_utils.R` — canonical `genotype_from_subtype()` helper (2k1b-aware); subtype/genotype derivation used by ranking + dedup.
- `modules/local/parsefirstmapping/main.nf` — PARSEFIRSTMAPPING process; emits `major_mapping`/`minor_mapping`/`csv` channels + stub contract. Long-format output + per-candidate emits land here.
- `workflows/hcvtyper.nf` — `MAJOR_MAPPING`/`MINOR_MAPPING` aliases of `TARGETED_MAPPING` and the read-count/coverage filter routing (~lines 379–435); collapse to uniform per-candidate fan-out (D-04).
- `subworkflows/local/targeted_mapping/` — the uniform mapping path each candidate runs through; honor its metadata-key join discipline.
- `bin/summarize.R` — consumes the selection CSV (`major_reads`/`minor_reads`/`gate_flag`) and parses the mapping-output filename slot (3rd dot-field = major/minor). Legacy contract preserved in Phase 6 per D-06; full cutover in Phase 9.

### Conventions / patterns
- `.planning/codebase/CONVENTIONS.md`, `.planning/codebase/ARCHITECTURE.md` — module/script naming, channel patterns, R helper conventions.

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `genotype_from_subtype()` (`bin/genotype_utils.R`) — already staged into the PARSEFIRSTMAPPING workdir as a `path(genotype_utils)` input; reuse for subtype→genotype in ranking/dedup. Do NOT hardcode an absolute path (container portability).
- The empty-frame / NA-safe guards already in `summarize_mapping_to_all_references.R` (length-zero rowwise guard, `as.numeric` coercion of `minRead`/`minCov`, default `gate_flag="no_mapping"`) — generalize these to the candidate set rather than reinventing.
- `TARGETED_MAPPING` subworkflow already encapsulates the per-reference mapping+stats+consensus path — the single uniform path candidates fan out into.

### Established Patterns
- Selection logic emits decisions as CSV columns that Nextflow then routes on (R-emits, Nextflow-routes) — the v1.0 major-gate pattern. The long-format candidate table continues this: R ranks/emits rows, Nextflow fans out per row.
- Explicit `.join()` by metadata key in `TARGETED_MAPPING` (avoids index misalignment in parallel execution) — required for carrying `confirmation_status` per candidate (D-04).
- Typed params declared across schema + `conf/modules_hcv.config` + profiles (v1.0 `denovo_*` precedent) — model the candidate-count param the same way.

### Integration Points
- PARSEFIRSTMAPPING output → mapping fan-out → GET_MAPPING_STATS → SUMMARIZE. The long-format change ripples through the routing in `hcvtyper.nf` and the selection-CSV read in `summarize.R` (kept working via the D-06 shim).
- Downstream of mapping, the `_major.fa`/`_minor.fa` filename tag is the parsing key in `summarize.R` — the reason the rename is deferred to Phase 9.

</code_context>

<specifics>
## Specific Ideas

- Naming is locked at the milestone level: `cand_1..cand_n` for run-time candidates; `dominant` / `co-infection` / `background` for Phase 8 roles.
- Default candidate count = 2 (reproduces today's two-slot topology under the shim).
- The "guilty until corroborated" stance and the combined breadth-evenness-weighted dominance score are Phase 8 — Phase 6 only needs to make candidates first-class and neutrally ranked, and to stop filtering on co-infection validity.

</specifics>

<deferred>
## Deferred Ideas

- Combined dominance score (reads + k-mer coverage + breadth, breadth-evenness weighted) and strain-role classification — **Phase 8** (SCORE/CLASS).
- Re-homing the different-genotype / 1a-1b / 2k1b validity rules as classification rules — **Phase 8** (per D-05).
- Filename slot migration `.major.`/`.minor.` → `.cand1.`/`.cand2.`, `summarize.R` parsing cutover, and one-release legacy column aliasing — **Phase 9** (COMPAT-02/03).
- Per-genotype assembly-support evidence and its genotype-level join to candidates — **Phase 7** (ASUP).

### Reviewed Todos (not folded)
- `2026-06-06-refactor-hcvglue-parallel-docker.md` — matched only on generic keywords (modules/local/parallel); it is the **paused v2.0 Phase 5** work, unrelated to candidate selection. Not folded.
- `2026-06-06-remove-tanoti-mapper-completely.md` — matched on generic keywords (modules/local/mapping); Tanoti removal already shipped (per PROJECT.md). Not folded.

</deferred>

---

*Phase: 06-neutral-candidate-selection*
*Context gathered: 2026-06-12*
