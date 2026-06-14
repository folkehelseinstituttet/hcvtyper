# Phase 8: Dominance Scoring + Strain-Role Classification - Context

**Gathered:** 2026-06-13
**Status:** Ready for planning

<domain>
## Phase Boundary

At the **summary step only** (`bin/summarize.R`), assign each neutral Phase-6 candidate a combined, breadth-evenness-weighted **dominance score**, classify each candidate into exactly one strain role — **dominant / co-infection / background** ("guilty until corroborated") — surface background/artefact candidates explicitly with a reason (never silently drop), and derive one **overall sample call** (monoinfection / co-infection / indeterminate) per sample. This is where the index-hopping / cross-mapping false co-infections get caught.

**Covers requirements:** SCORE-01, SCORE-02 (combined breadth-evenness-weighted dominance score), CLASS-01 (per-candidate role), CLASS-02 (co-infection = abundance floor AND assembly support, else background), CLASS-03 (background surfaced with explicit reason), CLASS-04 (overall sample call from roles).

**Explicitly NOT in this phase:**
- Filename-slot migration `.major.`/`.minor.` → `.cand1.`/`.cand2.` and the `summarize.R` filename-parsing cutover → **Phase 9** (COMPAT-02).
- Legacy `Major_*`/`Minor_*` **output-column** aliasing for one release → **Phase 9** (COMPAT-03). (Phase 8 retires the legacy confirmation *logic*; the legacy wide columns remain and are aliased in Phase 9.)
- The CI regression suite covering ranking/roles/score/golden reproduction → **Phase 9** (TEST-01). Phase 8 adds unit-level coverage of the new functions; the full golden-baseline reproduction is gated at Phase 9.

</domain>

<decisions>
## Implementation Decisions

### Dominance score + dominant determination (SCORE-01, SCORE-02)
- **D-01 (dominant = top score, gated):** Exactly one candidate is `dominant` = the highest-scoring candidate **that also passes the major-gate** (minRead/minCov). If no candidate passes the gate, there is **no dominant** → sample call `indeterminate`. Preserves the v1.0 major-gate (handles ERR1810469's failed 3a major).
- **D-02 (score form = normalized weighted sum):** Dominance score = weighted sum of normalized/transformed components — `log10(reads)`, breadth fraction (0–1), a CV-evenness factor (0–1), and a `log(k-mer cov)` term — with the **largest weight on breadth-evenness**. Read counts alone are contamination-prone (the false 4g had more reads than many genuine minors), so evenness must dominate.
- **D-03 (evenness metric = CV over full reference):** Breadth-evenness is the **coefficient of variation (sd/mean) of per-position depth INCLUDING zero positions**, across the full reference length, mapped to a 0–1 factor (e.g. `1/(1+CV)`). Penalizes both patchiness and low breadth in one number — a 67%-breadth spiky 4g scores low. Per-position depth is available from `depth/*.tsv` (`cov$X3` in `summarize.R` ~L372).
- **D-04 (weights = calibrated params):** Expose the weights (and the evenness-transform constant) as **named params** (e.g. `score_weight_evenness` / `score_weight_reads` / `score_weight_kmercov`) with **defaults calibrated by the researcher against the handoff evidence table** so the benchmark outcomes reproduce. Matches the `denovo_*` param precedent and the manuscript's before/after toggling need.
- **D-05 (k-mer cov = bonus-only):** The k-mer-cov term contributes **only positively** (a capped boost); `assembly_support = "none"` / NA k-mer cov yields **no boost, never a penalty**. A genuine low-yield dominant that failed de novo is not demoted — honors the v1.0 fall-back rule.
- **D-06 (deterministic tie-break):** On near-equal top scores (balanced co-infection, e.g. sim1 1a:1b), highest score wins with a **deterministic tie-break** (then by mapped reads, then ref name) for reproducibility. The overall call is co-infection either way, so which candidate is labeled `dominant` is cosmetic. (No co-dominant flag.)

### Co-infection abundance floor (CLASS-02)
- **D-07 (floor = reuse minRead/minCov):** A non-dominant candidate's abundance floor = the **existing minRead/minCov** applied to that candidate (same as the dominant gate). Thresholds **cannot** separate artefact from genuine (the false 4g cleared them easily); the floor is just a minimal presence gate, and **de novo corroboration — not a stricter cutoff — is the discriminator** ("guilty until corroborated").
- **D-08 (coverage source = targeted/second mapping):** The floor (and the breadth-evenness score) read each candidate's **targeted (second) mapping** coverage — `cov_breadth` + `avg_depth` + per-position depth, and the `min_targeted_read`/`min_targeted_cov` gate. Floor and score therefore share one coverage source.
- **D-09 (one threshold set):** The dominant major-gate and the co-infection floor use the **same minRead/minCov** (one threshold concept; reproduces v1.0 which gated both on the same numbers). No separate lenient co-infection floor.

### Corroboration verdict + HCV exceptions (CLASS-02, CLASS-03, COMPAT-04)
- **D-10 (corroboration = reuse denovo_* floors):** "Has genotype-level assembly support" = the candidate's joined Phase-7 assembly-support metrics clear the **existing** `denovo_min_contig_length` (1000) AND `denovo_min_kmer_cov` (2.0) AND `denovo_min_blast_identity` (90), at `denovo_match_level` (default genotype). Already calibrated (03-RESEARCH); confirms ERR1810453's 2,949 bp partial 2b, ignores ~300 bp noise contigs.
- **D-11 (asymmetric refute rule):** A non-dominant candidate that **clears the floor but has no assembly support of its own** is classified `background` **only if de novo produced a substantial contig for the DOMINANT** (de novo demonstrably worked, so absence is real evidence — the 4g case). If de novo assembled nothing substantial for the dominant either, de novo is **inconclusive**: keep the candidate as `co-infection` flagged **`uncorroborated_kept`**, do not suppress. Preserves genuine low-yield/IVT minors. Matches the v1.0 asymmetric-refute Key Decision. ⚠️ This **tempers the literal CLASS-02 wording** ("no support → background") — verifier must check against the asymmetric rule, not strict suppression.
- **D-12 (exceptions = port `is_valid_minor()` verbatim):** Lift the existing `is_valid_minor()` logic into the Phase-8 role classifier as **post-scoring special-cases**: a co-infection requires a **different genotype** than the dominant, **except 1a/1b** (same gt1, different subtype, allowed as co-infection), and **2k1b** recombinant pairs are **blocked**. A same-genotype (non-1a/1b) candidate or a 2k1b pair → `background`. Preserves COMPAT-04 exactly. (D-05/Phase-6 moved this filtering out of selection into here.)
- **D-13 (role reason = coded field + derived sentence):** Each candidate carries a machine-stable **`role_reason`** with a controlled vocabulary — e.g. `corroborated` / `below_floor` / `refuted_denovo` / `same_genotype_as_dominant` / `recombinant_2k1b` / `uncorroborated_kept` — **plus** a human-readable sentence **derived** from it (reusing the `review_flag` sentence style from quick-task 260610-8t0). Stable for CI assertions and analyst-readable. (Final vocabulary subject to planner refinement.)

### Overall sample call + legacy migration (CLASS-04)
- **D-14 (3-value call, nuance in the flag):** Overall sample call is exactly the 3 spec values. **≥1 co-infection role** (including `uncorroborated_kept`) → `co-infection`; **1 dominant + only background** → `monoinfection`; **no candidate passes the major-gate** → `indeterminate`. Background candidates do not change the call (surfaced, not counted). The uncorroborated nuance lives in the per-candidate `role_reason`, not a 4th sample-call value (matches criterion #5).
- **D-15 (retire legacy confirmation logic now):** Phase 8 owns the cutover. The new corroboration verdict **replaces** `classify_minor_denovo` / `minor_denovo_status`, and `review_flag` is **rewired onto the new roles + `role_reason`**. No two confirmation systems coexisting. Phase 9 stays focused on filename migration + `Major_*`/`Minor_*` *column* aliasing + the regression suite. ⚠️ Departs from the Phase 6/7 additive precedent deliberately — avoids contradictory confirmation columns in `Summary.csv`.
- **D-16 (wide Summary + long candidate detail):** `Summary.csv` stays **one-row-per-sample (wide)**: dominant + corroborated co-infection fill the existing `Major_*`/`Minor_*` slots (aliased in Phase 9) plus the new overall sample-call column. The **full per-candidate role + dominance score + `role_reason`, including every background candidate**, lands in the long-format `*.candidates.csv` (Phase 6 emit) so backgrounds are surfaced not dropped (CLASS-03). Preserves the wide-Summary contract Phase 9 column aliasing assumes.

### Claude's Discretion
- Exact column names/types for the new score, role, and `role_reason` fields (follow `denovo_*` / `cand_*` naming precedent and the legacy-retention constraints).
- Exact normalization constants/transforms within D-02 (e.g. log base, breadth source `min_1` vs `min_10`, how the k-mer-cov bonus is capped) — researcher to fit against the evidence table.
- Whether the role classifier is a new sourceable R helper (e.g. `bin/classify_roles.R`) or inline in `summarize.R` — pick the least-disruptive shape consistent with the existing `denovo_confirm.R` / `assembly_support_join.R` pattern.
- Precise measurement of "de novo worked for the dominant" in D-11 (reuse the dominant's joined assembly-support metrics against the same `denovo_*` floors).
- Final `role_reason` controlled vocabulary (D-13) and the exact derived-sentence wording.

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Phase scope & requirements
- `.planning/ROADMAP.md` §"Phase 8: Dominance Scoring + Strain-Role Classification" — goal, 5 success criteria, dependency on Phase 6 AND Phase 7.
- `.planning/REQUIREMENTS.md` — SCORE-01/02, CLASS-01..04 (Phase 8); plus COMPAT-01..04 / TEST-01 (Phase 9) for awareness of what Phase 8 must NOT finish.
- `.planning/PROJECT.md` — Core Value, the "substantial contig" calibration constraint (≥~1,000 bp / identity / k-mer floor), genotype-level match constraint, non-breaking constraint, and the Key Decisions table (asymmetric refute rule, genotype-level match).
- `hcvtyper_handoff_denovo_informed_selection.md` (repo root) — **§2 evidence table is the calibration anchor** (false 4g → refute/background; ERR1810447 full 2b → co-infection; ERR1810453 2,949 bp partial → confirm; IVT extreme-ratio genuine minors must NOT be suppressed). §4 asymmetric fall-back rule, §5 substantial-contig thresholds + genotype-vs-subtype match.

### Upstream phase decisions this phase builds on
- `.planning/phases/06-neutral-candidate-selection/06-CONTEXT.md` — D-01 (long-format `*.candidates.csv` = the classification input), D-03 (2nd candidate ranked by reads, may differ from old minor), **D-05 (the different-genotype/1a-1b/2k1b validity rules were moved OUT of selection INTO Phase 8 — this phase re-applies them as role logic)**, D-06 (legacy shim retained until Phase 9).
- `.planning/phases/07-per-genotype-assembly-support/07-CONTEXT.md` — D-01 (assembly support = raw metrics only; **the substantiality floor verdict is Phase 8's job**), D-02/D-03 (per-subtype support rolled to `denovo_match_level` at join; single best contig by `sc_length`), and the deferred list explicitly handing the corroboration verdict + role assignment + legacy-path retirement to Phase 8.

### Code touchpoints (the files this phase changes)
- `bin/summarize.R` — primary site. Reads the candidate CSV + assembly-support join (Phase 7) + coverage (`df_coverage`, ~L348-437, per-position depth at ~L372) + confirmation layer (~L797-815) + builds `review_flag` (~L938+). Add: dominance score, role classification, overall sample call; **retire** the `classify_minor_denovo` / `minor_denovo_status` path and rewire `review_flag` (D-15).
- `bin/denovo_confirm.R` — `classify_minor_denovo()` is the legacy minor-coupled verdict to be retired (D-15); its threshold floors (`denovo_min_*`) are reused as the corroboration verdict (D-10).
- `bin/summarize_mapping_to_all_references.R` — contains `is_valid_minor()` (different-genotype / 1a-1b / 2k1b). Port that logic into the Phase-8 role classifier verbatim (D-12).
- `bin/genotype_utils.R` — `genotype_from_subtype()` (2k1b-aware); used to compare candidate genotype vs dominant genotype for D-12 and the corroboration match level.
- The long-format `*.candidates.csv` (Phase 6 emit) — gains per-candidate role / dominance-score / `role_reason` columns (D-16).
- `conf/modules_hcv.config` (~L328) — already passes `denovo_min_* denovo_match_level minRead minCov min_targeted_read min_targeted_cov` to `summarize.R`. Add the new `score_weight_*` params here + the nf-schema + profiles (D-04), mirroring the `denovo_*` param wiring.

### Conventions / patterns
- `.planning/codebase/CONVENTIONS.md`, `.planning/codebase/ARCHITECTURE.md` — R helper conventions, R-emits/Nextflow-routes pattern, typed-empty-tibble guards.

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- **Per-position depth is already in hand:** `summarize.R` reads `depth/*.tsv` into `cov` and uses `cov$X3` (per-position depth) at ~L372-407 to compute `avg_depth` and breadth (positions ≥1/≥5/≥10 ÷ ref_length). The CV-evenness factor (D-03) is computed directly from this same `cov$X3` vector — no new process needed.
- **`is_valid_minor()`** in `bin/summarize_mapping_to_all_references.R` — the validated 1a/1b allow + 2k1b block + different-genotype rule to port verbatim (D-12).
- **`denovo_*` floors + `classify_minor_denovo()`** in `bin/denovo_confirm.R` — calibration-validated thresholds reused as the corroboration verdict (D-10); the function body itself is retired (D-15).
- **`genotype_from_subtype()`** (`bin/genotype_utils.R`, 2k1b-aware) — already sourced/staged in `summarize.R`; reuse for genotype comparisons.
- **`review_flag` sentence builder** (`summarize.R` ~L938+, refined by quick-tasks 260609/260610) — rewire onto roles for the derived human sentence (D-13, D-15).
- **Typed-empty-tibble guards (PLUMB-02 / T-03-01)** — generalize to the new score/role columns so skip-assembly / no_mapping runs NA-fill rather than abort.

### Established Patterns
- **R-emits decisions, Nextflow-routes** — classification stays in R (`summarize.R`); no new Nextflow routing needed (this is the terminal summary step).
- **`left_join` + typed-empty fallback** — the Phase-7 assembly-support join is already the classifier's evidence input; roles compute over the joined frame.
- **Params declared across nf-schema + `conf/modules_hcv.config` + profiles** (v1.0 `denovo_*` precedent) — model `score_weight_*` the same way (D-04).

### Integration Points
- Phase 6 `*.candidates.csv` (long) + Phase 7 assembly-support join + `df_coverage` (per-candidate targeted-mapping breadth/depth) → **Phase 8 role classifier in `summarize.R`** → wide `Summary.csv` (sample call + Major_/Minor_ slots) + enriched long `*.candidates.csv` (all roles incl. background). This is the read path roles slot into.

</code_context>

<specifics>
## Specific Ideas

- **Evidence-table calibration targets (handoff §2) are the verification anchor:** false 4g (53,279 reads, 67.7% breadth, no 4g contig) → `background`/`refuted_denovo`; ERR1810447 full 9,207 bp 2b → `co-infection`/`corroborated`; ERR1810453 2,949 bp partial 2b → `co-infection`/`corroborated`; sim1 1a:1b + sim2 2a:3a + ERR1810505 + ERR1810511 (true co-infections) → preserved; IVT extreme-ratio genuine minors → not suppressed (`uncorroborated_kept` if de novo failed).
- **Naming locked at milestone level:** roles are exactly `dominant` / `co-infection` / `background`; sample call is exactly `monoinfection` / `co-infection` / `indeterminate`.
- **Breadth-evenness is the headline discriminator** (SCORE-02) — it must out-weigh raw read count strongly enough that the spiky high-read 4g loses to a genuine even minor.

</specifics>

<deferred>
## Deferred Ideas

- A 4th `possible_coinfection` sample-call value for the uncorroborated case — rejected to keep criterion #5's 3-value contract; nuance carried in `role_reason` instead. Revisit only with roadmap sign-off.
- A co-dominant / balanced-co-infection flag for near-tie top scores — deferred (D-06 deterministic tie-break chosen); revisit if calibration shows balanced cases are common.
- Config-driven `coinfection_allow_pairs` / `block_pairs` generalization of the HCV exceptions — deferred in favor of the verbatim `is_valid_minor()` port (D-12); a future extensibility option.
- Filename-slot migration, `Major_*`/`Minor_*` column aliasing, and the full CI golden-baseline regression suite — **Phase 9** (COMPAT-02/03, TEST-01).

</deferred>

---

*Phase: 08-dominance-scoring-strain-role-classification*
*Context gathered: 2026-06-13*
