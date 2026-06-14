# Phase 8: Dominance Scoring + Strain-Role Classification - Research

**Researched:** 2026-06-13
**Domain:** R data-classification logic inside `bin/summarize.R` (HCV strain-role assignment); no new libraries, no new Nextflow routing
**Confidence:** HIGH (all code touchpoints read in-session; calibration anchor is the in-repo handoff evidence table)

## Summary

This is a **brownfield, R-logic-only** phase. Every input it consumes already exists and was verified in-session: the Phase-6 long-format `*.candidates.csv` (8 typed columns), the Phase-7 `join_assembly_support()` output (wide `cand_<rank>_assembly_support_*` block, already left-joined into `final`), and per-candidate per-position depth (`cov$X3`) inside the existing coverage loop at `summarize.R:360-413`. No new process, no new channel, no new container. The work is: compute a per-candidate dominance score, classify each candidate into exactly one role, surface backgrounds, derive one sample call, **retire** the legacy `classify_minor_denovo`/`minor_denovo_status`/`apply_denovo_layer` path (D-15), and rewire `review_flag` onto the new roles.

The phase has an unusually complete `<decisions>` block (D-01..D-16) — research did not need to choose an architecture, it needed to **verify the decisions are buildable against the real code** and **calibrate the score weights against the handoff evidence table** (the one genuine open research task, D-04). All decisions are buildable as written. Two latent risks surfaced that the planner must address explicitly (see Common Pitfalls): a **calibration-default mismatch** between `nextflow.config` (denovo floors ship as 500 / 10.0) and the 03-RESEARCH-validated 1000 / 2.0 that D-10 assumes; and the **positional-arg coupling** in the SUMMARIZE `ext.args` string that any new `score_weight_*` param must slot into without shifting existing positions.

**Primary recommendation:** Implement the role classifier and dominance score as **new pure sourced R helpers** (`bin/classify_roles.R` + a small `bin/dominance_score.R`, or one combined file) mirroring the `assembly_support_join.R` / `denovo_layer.R` convention (pure functions, no I/O, no `commandArgs`, defensive tidyverse guard), source them into `summarize.R`, and add `bin/tests/test_*.R` unit tests in the existing self-contained harness (`source()` + in-memory tibbles + `fail()`/`ok()`, run by `bin/tests/run_all.sh`). Compute the CV-evenness factor inside the existing `summarize.R:360-413` cov loop from `cov$X3` (no new process). Calibrate `score_weight_*` defaults so the handoff §2 evidence table reproduces (false 4g → background; ERR1810447/ERR1810453 2b → co-infection; true co-infections preserved).

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

**Dominance score + dominant determination (SCORE-01, SCORE-02)**
- **D-01 (dominant = top score, gated):** Exactly one candidate is `dominant` = the highest-scoring candidate **that also passes the major-gate** (minRead/minCov). If no candidate passes the gate, there is **no dominant** → sample call `indeterminate`. Preserves the v1.0 major-gate (handles ERR1810469's failed 3a major).
- **D-02 (score form = normalized weighted sum):** Dominance score = weighted sum of normalized/transformed components — `log10(reads)`, breadth fraction (0–1), a CV-evenness factor (0–1), and a `log(k-mer cov)` term — with the **largest weight on breadth-evenness**. Read counts alone are contamination-prone (the false 4g had more reads than many genuine minors), so evenness must dominate.
- **D-03 (evenness metric = CV over full reference):** Breadth-evenness is the **coefficient of variation (sd/mean) of per-position depth INCLUDING zero positions**, across the full reference length, mapped to a 0–1 factor (e.g. `1/(1+CV)`). Penalizes both patchiness and low breadth in one number — a 67%-breadth spiky 4g scores low. Per-position depth is available from `depth/*.tsv` (`cov$X3` in `summarize.R` ~L372).
- **D-04 (weights = calibrated params):** Expose the weights (and the evenness-transform constant) as **named params** (e.g. `score_weight_evenness` / `score_weight_reads` / `score_weight_kmercov`) with **defaults calibrated by the researcher against the handoff evidence table** so the benchmark outcomes reproduce. Matches the `denovo_*` param precedent and the manuscript's before/after toggling need.
- **D-05 (k-mer cov = bonus-only):** The k-mer-cov term contributes **only positively** (a capped boost); `assembly_support = "none"` / NA k-mer cov yields **no boost, never a penalty**. A genuine low-yield dominant that failed de novo is not demoted — honors the v1.0 fall-back rule.
- **D-06 (deterministic tie-break):** On near-equal top scores (balanced co-infection, e.g. sim1 1a:1b), highest score wins with a **deterministic tie-break** (then by mapped reads, then ref name) for reproducibility. The overall call is co-infection either way, so which candidate is labeled `dominant` is cosmetic. (No co-dominant flag.)

**Co-infection abundance floor (CLASS-02)**
- **D-07 (floor = reuse minRead/minCov):** A non-dominant candidate's abundance floor = the **existing minRead/minCov** applied to that candidate (same as the dominant gate). Thresholds **cannot** separate artefact from genuine; the floor is just a minimal presence gate, and **de novo corroboration — not a stricter cutoff — is the discriminator** ("guilty until corroborated").
- **D-08 (coverage source = targeted/second mapping):** The floor (and the breadth-evenness score) read each candidate's **targeted (second) mapping** coverage — `cov_breadth` + `avg_depth` + per-position depth, and the `min_targeted_read`/`min_targeted_cov` gate. Floor and score therefore share one coverage source.
- **D-09 (one threshold set):** The dominant major-gate and the co-infection floor use the **same minRead/minCov** (one threshold concept; reproduces v1.0 which gated both on the same numbers). No separate lenient co-infection floor.

**Corroboration verdict + HCV exceptions (CLASS-02, CLASS-03, COMPAT-04)**
- **D-10 (corroboration = reuse denovo_* floors):** "Has genotype-level assembly support" = the candidate's joined Phase-7 assembly-support metrics clear the **existing** `denovo_min_contig_length` (1000) AND `denovo_min_kmer_cov` (2.0) AND `denovo_min_blast_identity` (90), at `denovo_match_level` (default genotype). Already calibrated (03-RESEARCH); confirms ERR1810453's 2,949 bp partial 2b, ignores ~300 bp noise contigs.
- **D-11 (asymmetric refute rule):** A non-dominant candidate that **clears the floor but has no assembly support of its own** is classified `background` **only if de novo produced a substantial contig for the DOMINANT** (de novo demonstrably worked, so absence is real evidence — the 4g case). If de novo assembled nothing substantial for the dominant either, de novo is **inconclusive**: keep the candidate as `co-infection` flagged **`uncorroborated_kept`**, do not suppress. Preserves genuine low-yield/IVT minors. Matches the v1.0 asymmetric-refute Key Decision. ⚠️ This **tempers the literal CLASS-02 wording** ("no support → background") — verifier must check against the asymmetric rule, not strict suppression.
- **D-12 (exceptions = port `is_valid_minor()` verbatim):** Lift the existing `is_valid_minor()` logic into the Phase-8 role classifier as **post-scoring special-cases**: a co-infection requires a **different genotype** than the dominant, **except 1a/1b** (same gt1, different subtype, allowed as co-infection), and **2k1b** recombinant pairs are **blocked**. A same-genotype (non-1a/1b) candidate or a 2k1b pair → `background`. Preserves COMPAT-04 exactly. (D-05/Phase-6 moved this filtering out of selection into here.)
- **D-13 (role reason = coded field + derived sentence):** Each candidate carries a machine-stable **`role_reason`** with a controlled vocabulary — e.g. `corroborated` / `below_floor` / `refuted_denovo` / `same_genotype_as_dominant` / `recombinant_2k1b` / `uncorroborated_kept` — **plus** a human-readable sentence **derived** from it (reusing the `review_flag` sentence style from quick-task 260610-8t0). Stable for CI assertions and analyst-readable. (Final vocabulary subject to planner refinement.)

**Overall sample call + legacy migration (CLASS-04)**
- **D-14 (3-value call, nuance in the flag):** Overall sample call is exactly the 3 spec values. **≥1 co-infection role** (including `uncorroborated_kept`) → `co-infection`; **1 dominant + only background** → `monoinfection`; **no candidate passes the major-gate** → `indeterminate`. Background candidates do not change the call (surfaced, not counted). The uncorroborated nuance lives in the per-candidate `role_reason`, not a 4th sample-call value (matches criterion #5).
- **D-15 (retire legacy confirmation logic now):** Phase 8 owns the cutover. The new corroboration verdict **replaces** `classify_minor_denovo` / `minor_denovo_status`, and `review_flag` is **rewired onto the new roles + `role_reason`**. No two confirmation systems coexisting. Phase 9 stays focused on filename migration + `Major_*`/`Minor_*` *column* aliasing + the regression suite. ⚠️ Departs from the Phase 6/7 additive precedent deliberately — avoids contradictory confirmation columns in `Summary.csv`.
- **D-16 (wide Summary + long candidate detail):** `Summary.csv` stays **one-row-per-sample (wide)**: dominant + corroborated co-infection fill the existing `Major_*`/`Minor_*` slots (aliased in Phase 9) plus the new overall sample-call column. The **full per-candidate role + dominance score + `role_reason`, including every background candidate**, lands in the long-format `*.candidates.csv` (Phase 6 emit) so backgrounds are surfaced not dropped (CLASS-03). Preserves the wide-Summary contract Phase 9 column aliasing assumes.

### Claude's Discretion
- Exact column names/types for the new score, role, and `role_reason` fields (follow `denovo_*` / `cand_*` naming precedent and the legacy-retention constraints).
- Exact normalization constants/transforms within D-02 (e.g. log base, breadth source `min_1` vs `min_10`, how the k-mer-cov bonus is capped) — researcher to fit against the evidence table.
- Whether the role classifier is a new sourceable R helper (e.g. `bin/classify_roles.R`) or inline in `summarize.R` — pick the least-disruptive shape consistent with the existing `denovo_confirm.R` / `assembly_support_join.R` pattern.
- Precise measurement of "de novo worked for the dominant" in D-11 (reuse the dominant's joined assembly-support metrics against the same `denovo_*` floors).
- Final `role_reason` controlled vocabulary (D-13) and the exact derived-sentence wording.

### Deferred Ideas (OUT OF SCOPE)
- A 4th `possible_coinfection` sample-call value for the uncorroborated case — rejected to keep criterion #5's 3-value contract; nuance carried in `role_reason` instead.
- A co-dominant / balanced-co-infection flag for near-tie top scores — deferred (D-06 deterministic tie-break chosen).
- Config-driven `coinfection_allow_pairs` / `block_pairs` generalization of the HCV exceptions — deferred in favor of the verbatim `is_valid_minor()` port (D-12).
- Filename-slot migration, `Major_*`/`Minor_*` column aliasing, full CI golden-baseline regression suite — **Phase 9** (COMPAT-02/03, TEST-01).
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|------------------|
| SCORE-01 | Combined score over mapped read count, k-mer coverage, mapping coverage breadth | All three inputs verified present: `candidate_reads`/`Reads_nodup_mapped_*` (reads), `cand_<rank>_assembly_support_best_contig_kmer_cov` (k-mer cov, from Phase-7 join), `cov_breadth_*` + `cov$X3` (breadth). Score = D-02 weighted sum computed in a new sourced helper. |
| SCORE-02 | Breadth evenness weighted strongly vs raw read count | CV-evenness factor (D-03) computed from `cov$X3` (per-position depth incl. zeros) inside the existing `summarize.R:360-413` cov loop; `score_weight_evenness` default > other weights, calibrated against handoff §2. |
| CLASS-01 | Each candidate classified dominant / co-infection / background | Classifier reads the Phase-7-joined long candidate frame; emits one role per candidate into the long `*.candidates.csv` (D-16) + wide `Major_*`/`Minor_*` slots. |
| CLASS-02 | Co-infection only if abundance floor AND assembly support, else background ("guilty until corroborated") | Floor = `minRead`/`minCov` reused (D-07/D-09); support = `assembly_support` from the Phase-7 join cleared against `denovo_*` floors (D-10); asymmetric refute (D-11) tempers literal wording. |
| CLASS-03 | Background surfaced with explicit reason, never dropped | `role_reason` coded field + derived sentence (D-13); every candidate incl. background written to long `*.candidates.csv` (D-16). Candidates are the LEFT side of every join → no row loss (verified in `assembly_support_join.R`). |
| CLASS-04 | Overall sample call from roles | 3-value derivation (D-14) over the per-candidate roles; one new wide column on `Summary.csv`. |
</phase_requirements>

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Dominance score computation | R logic (`bin/summarize.R` + new helper) | — | Terminal summary step; all evidence already joined into `final`/`candidate_support` here |
| CV-evenness from per-position depth | R logic (`summarize.R` cov loop) | — | `cov$X3` only exists inside the L360-413 loop; must compute CV there before discarding |
| Role classification | R logic (new sourced helper) | — | Pure function over the joined candidate frame; mirrors `assembly_support_join.R` |
| HCV exception rules (1a/1b, 2k1b) | R logic (ported `is_valid_minor()`) | — | Verbatim port (D-12); `genotype_from_subtype()` already sourced |
| Overall sample call | R logic (`summarize.R`) | — | Aggregation over per-candidate roles (D-14) |
| Param wiring (`score_weight_*`) | Nextflow config + nf-schema | R arg parse | `nextflow.config` defaults + `nextflow_schema.json` + `conf/modules_hcv.config` ext.args + `summarize.R` `commandArgs` (mirror `denovo_*`) |
| Routing / channels | **None — no change** | — | "R-emits decisions, Nextflow-routes"; this is the terminal step, nothing routes downstream |

## Standard Stack

No new libraries. Phase 8 is pure R inside the existing SUMMARIZE container.

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| R (base) | 4.x (SUMMARIZE container) | Scripting, `sd`/`mean`/`log10` for the score | Already the summarize runtime [VERIFIED: modules/local/summarize/main.nf container] |
| tidyverse | 3536dd50a17de0ab wave build | dplyr/tidyr/readr/purrr for the classifier joins + `pmap_chr` sentence builder | Every existing `bin/*.R` helper uses it [VERIFIED: bin/summarize.R:3] |
| seqinr | (in container) | FASTA read/write (unchanged; not touched by Phase 8) | Already present [VERIFIED: container env] |

**Installation:** None. The SUMMARIZE container (`community.wave.seqera.io/library/r-gridextra_r-png_r-seqinr_r-tidyverse:3536dd50a17de0ab`) already ships R + tidyverse + seqinr + gridExtra + png. [VERIFIED: modules/local/summarize/main.nf]

### Supporting
| Asset | Purpose | When to Use |
|-------|---------|-------------|
| `bin/genotype_utils.R::genotype_from_subtype()` | 2k1b-aware genotype extraction | Comparing candidate genotype vs dominant genotype (D-12); already sourced in summarize.R:10 |
| `bin/tests/run_all.sh` | R unit-test runner (globs `test_*.R`, exits non-zero on any fail) | New `test_classify_roles.R` / `test_dominance_score.R` are auto-picked-up |

**Alternatives considered:** None — the decisions block locks the approach. Do not introduce `testthat` (the existing harness is plain `Rscript` + `fail()`/`ok()`; matching it keeps `run_all.sh` working).

## Package Legitimacy Audit

> Not applicable — Phase 8 installs **no** external packages. All R libraries are pre-provisioned in the existing SUMMARIZE container. No `npm`/`pip`/`cargo`/`conda` install occurs in this phase.

## Architecture Patterns

### System Architecture Diagram (data flow into the Phase-8 classifier)

```
PARSEFIRSTMAPPING (Phase 6)                 TARGETED_MAPPING (per candidate)
  ├─ *.candidates.csv  ────────┐              └─ depth/<id>.<ref>.{major|minor}.nodup….tsv
  │   (8 typed cols, long)     │                   │ (per-position depth, X3)
  └─ *.parsefirstmapping.csv   │                   │
      (legacy wide shim)       │                   ▼
                               │            summarize.R cov loop (L360-413)
SPADES+BLAST (Phase 7)         │              ├─ breadth (cov_breadth_min_1/5/10)
  └─ *.assembly_support.csv ───┤              └─ NEW: CV-evenness from cov$X3  ◄── D-03
      (per-subtype contig)     │                   (compute HERE, before X3 discarded)
                               ▼
                  join_assembly_support()  (Phase 7, ASUP-02, genotype level)
                               │  candidates = LEFT side → no row loss
                               ▼
              candidate_support  (long) ── pivot_wider ──► cand_<rank>_assembly_support_* (wide)
                               │                                  │
                               ▼                                  ▼
        ┌──────────────────────────────────────────────────────────────────┐
        │  NEW Phase-8 role classifier (sourced helper, pure)               │
        │   1. dominance score per candidate (D-02): log10(reads),          │
        │      breadth, CV-evenness factor, log(kmer_cov) bonus-only        │
        │   2. dominant = top score AND passes minRead/minCov gate (D-01)   │
        │   3. each non-dominant: floor(minRead/minCov) AND support? →      │
        │      co-infection ; asymmetric refute (D-11) → background ;       │
        │      no-support-but-dominant-also-unassembled → uncorroborated_kept│
        │   4. HCV exceptions ported from is_valid_minor() (D-12)           │
        │   5. role_reason coded vocab + derived sentence (D-13)            │
        │   6. overall sample call (D-14)                                   │
        └──────────────────────────────────────────────────────────────────┘
                               │                                  │
                               ▼                                  ▼
         enriched *.candidates.csv (long)            Summary.csv (wide, one row/sample)
         every candidate incl. background            dominant→Major_*, co-inf→Minor_*,
         + role + score + role_reason (CLASS-03)     + overall_sample_call (CLASS-04)
                                                     + rewired review_flag (D-15)
  RETIRE: classify_minor_denovo / minor_denovo_status / apply_denovo_layer (D-15)
```

### Recommended File Structure
```
bin/
├── summarize.R                      # consumes the new helpers; retires the legacy path; rewires review_flag
├── dominance_score.R                # NEW pure helper: score_candidates(df, weights, evenness_const) — OR fold into classify_roles.R
├── classify_roles.R                 # NEW pure helper: classify_roles(scored_df, minRead, minCov, denovo_floors, match_level)
├── genotype_utils.R                 # UNCHANGED (sourced; genotype_from_subtype)
├── assembly_support_join.R          # UNCHANGED (Phase 7; produces the support input)
├── denovo_confirm.R                 # RETIRE consumption (D-15); keep file or delete per planner
├── denovo_layer.R                   # RETIRE consumption (D-15)
└── tests/
    ├── run_all.sh                   # UNCHANGED runner (auto-globs new tests)
    ├── test_dominance_score.R       # NEW: score monotonicity + evenness-dominates-reads (the 4g case)
    ├── test_classify_roles.R        # NEW: all 3 roles, asymmetric refute, exceptions, sample call
    └── fixtures/                    # NEW evidence-table fixtures (4g refute, 2b confirm, true co-infections)
```

### Pattern 1: Pure Sourced R Helper (the house pattern)
**What:** A `bin/*.R` file that defines functions only — no `commandArgs`, no file I/O, no global mutation — guarded by a defensive `if (!exists("group_by")) library(tidyverse)`. Staged into the SUMMARIZE task workdir as a declared `path()` input and `source()`d relatively.
**When to use:** For every new piece of Phase-8 decision logic, so it is unit-testable without a pipeline run.
**Example:**
```r
# Source: bin/assembly_support_join.R:36-50 (verbatim house pattern)
if (!exists("group_by")) {
  library(tidyverse)
}
classify_roles <- function(candidates_df, minRead, minCov,
                           denovo_min_contig_length = 1000,
                           denovo_min_kmer_cov = 2.0,
                           denovo_min_blast_identity = 90,
                           match_level = "genotype",
                           score_weights = list(...)) {
  stopifnot(match_level %in% c("genotype", "subtype"))   # WR-04 guard, copy from line 50
  # ... pure logic; return the candidates frame + role/score/role_reason columns ...
}
```

### Pattern 2: CV-Evenness Inside the Existing Cov Loop
**What:** The per-position depth vector `cov$X3` only exists transiently inside `summarize.R:360-413`; it is read per candidate (one tsv per `.major.`/`.minor.` slot) to compute breadth, then discarded. The CV-evenness factor (D-03) MUST be captured inside that loop.
**When to use:** Add a column to `tmp_df` (e.g. `cv_evenness`) computed as `sd(cov$X3)/mean(cov$X3)` over the **full** reference (zeros included, which `cov$X3` already contains because samtools depth at `-aa` / the existing read produces one row per position), then map to a 0–1 factor `1/(1+CV)`.
**Example:**
```r
# Source: derived from summarize.R:372-413 (the existing cov loop)
cov <- read_tsv(cov_files[i], col_names = FALSE)
ref_length <- nrow(cov)
tmp_df$avg_depth[i] <- mean(cov$X3)
# NEW (D-03): CV over per-position depth incl. zeros, mapped to 0–1 evenness factor.
cv_raw <- if (ref_length > 0 && mean(cov$X3) > 0) sd(cov$X3) / mean(cov$X3) else NA_real_
tmp_df$cv_evenness[i] <- if (!is.na(cv_raw)) 1 / (1 + cv_raw) else 0
```
**✅ VERIFIED:** SAMTOOLS_DEPTH runs with `-aa` (`conf/modules_hcv.config:259` — `--count-orphans --no-BAQ --max-depth 0 --min-BQ 0 -aa`), so `cov$X3` **does** contain one row per reference position including zero-depth positions. The CV-over-full-reference (D-03) is therefore correct as designed — a low-breadth spiky candidate has many zero rows → high CV → low evenness factor. No breadth-folding workaround needed. [VERIFIED: conf/modules_hcv.config:259]

### Pattern 3: pmap_chr Sentence Builder (rewire, don't replace the style)
**What:** `review_flag` is built at `summarize.R:1096-1124` via `pmap_chr` over per-row trigger columns, emitting full human sentences joined with `" | "`. D-13/D-15 rewire this onto `role_reason` + roles.
**When to use:** Derive the human `role_reason` sentence and the rewired `review_flag` from the coded vocabulary using the same `pmap_chr` idiom (consistency with quick-task 260610-8t0).

### Anti-Patterns to Avoid
- **Resurrecting a suppressed candidate.** The legacy layer was downgrade-only by construction (`denovo_layer.R:62-66`). The new classifier assigns roles from scratch, so this invariant no longer holds automatically — the planner must ensure the asymmetric-refute (D-11) and exception (D-12) rules can only demote to `background`, never promote a below-floor candidate to co-infection.
- **Computing CV after the cov loop.** `cov$X3` is gone after L413. Computing CV from the surviving `avg_depth`/breadth scalars is impossible — it must be in the loop.
- **Inline re-implementing the helper in the test.** The whole `bin/tests/` convention asserts on the REAL sourced function (see `test_assembly_support_join.R` header). Do the same.
- **Two confirmation systems.** D-15 mandates retiring `minor_denovo_status`; do not leave it alongside the new `role`/`role_reason` columns (that is the explicit Phase-6/7-precedent departure).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Genotype from subtype (2k1b-aware) | New `substr`/`if` logic | `genotype_from_subtype()` (genotype_utils.R) | Already sourced, 2k1b-correct, single-sourced canonical rule |
| 1a/1b allow + 2k1b block + different-genotype | New rule from scratch | Port `is_valid_minor()` verbatim (summarize_mapping_to_all_references.R) | D-12 mandates verbatim port; preserves COMPAT-04 exactly |
| Genotype-level assembly-support match | Re-join contigs to candidates | Reuse `join_assembly_support()` output (already in `candidate_support` / wide `cand_*` cols) | Phase 7 already did the genotype-level join; classifier reads `assembly_support` |
| Substantial-contig floors | New threshold constants | Reuse `denovo_min_contig_length`/`kmer_cov`/`blast_identity` (D-10) | Already calibration-validated (03-RESEARCH) |
| Per-candidate breadth/depth | New samtools process | Existing `cov_breadth_*` + `cov$X3` in the cov loop (D-08) | Targeted-mapping depth already staged in `depth/` |
| Human-readable reason sentences | New formatter | `pmap_chr` sentence idiom (summarize.R:1096-1124) | Established style; MultiQC orange-highlight already wired for `review_flag` |
| Test harness | testthat / new runner | `bin/tests/*.R` + `run_all.sh` (source + fail/ok) | The CI guard globs `test_*.R`; matching keeps it green |

**Key insight:** Phases 6 and 7 deliberately produced the *evidence streams* (neutral candidates, genotype-level assembly support, per-candidate coverage) so Phase 8 is **pure combination logic** — almost nothing here is novel infrastructure; the risk is entirely in the *score calibration* and the *legacy-retirement surface*.

## Runtime State Inventory

> Phase 8 is logic-only (R inside one terminal process). It writes no datastore keys, registers no OS state, and renames no live-service config. The one rename-adjacent action is the **legacy-column retirement** (D-15), which is a code edit, not a data migration.

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| Stored data | None — pipeline is stateless per-run; outputs are files written fresh each run. | None |
| Live service config | None — no external service holds Phase-8 strings. | None |
| OS-registered state | None — no scheduled tasks / daemons. | None |
| Secrets/env vars | None. | None |
| Build artifacts | The SUMMARIZE container is pinned by digest; no rebuild needed (no new deps). New `bin/*.R` helpers must be added as `path()` inputs to `modules/local/summarize/main.nf` (mirror `denovo_confirm`/`assembly_support_join` staging) and to the `workflows/hcvtyper.nf` SUMMARIZE call site. | Code edit: stage new helper files |
| Legacy output columns (D-15) | `minor_denovo_status`, `coinfection_flag` (and the `denovo_*_subtype_match` columns feeding `review_flag`) are emitted in `Summary.csv` / `summary_mqc.tsv` and **hard-coded in the SUMMARIZE stub block** (`modules/local/summarize/main.nf` stub CSV header). The nf-test snapshot (`tests/default.nf.test.snap`) pins these. | Code edit: update stub header + regen snapshot. ⚠️ See MEMORY "Phase 7 SUMMARIZE snapshot regen deferred" — a snapshot regen is already pending; fold Phase-8 column changes into the same regen (needs docker/disk). |

**Nothing found in categories Stored data / Live service config / OS-registered state / Secrets — verified by reading the SUMMARIZE module, modules_hcv.config, and confirming the pipeline writes only run-local files.**

## Common Pitfalls

### Pitfall 1: Calibration-default mismatch (nextflow.config ships 500 / 10.0, D-10 assumes 1000 / 2.0)
**What goes wrong:** D-10 reuses `denovo_min_contig_length` / `denovo_min_kmer_cov` as the corroboration floor and cites the **validated** 1000 / 2.0. But `nextflow.config:46-47` ships `denovo_min_contig_length = 500` and `denovo_min_kmer_cov = 10.0`, `conf/modules_hcv.config` does **not** override them, and the SUMMARIZE `ext.args` passes the *param* values. The R-helper defaults (1000 / 2.0) only fire when the arg is empty — which it never is. So the **effective runtime floors are 500 / 10.0**, not the validated 1000 / 2.0.
**Why it happens:** 03-RESEARCH validated 1000 / 2.0 / 90 against real contig data, but the config defaults were left at the earlier 500 / 10.0 (Phase-1/2 values) and never corrected.
**How to avoid:** The researcher/planner must decide explicitly: either (a) correct `nextflow.config` (and `nextflow_schema.json`) to 1000 / 2.0 so the validated floors are the shipped behaviour, or (b) re-validate that 500 / 10.0 still produces the correct evidence-table verdicts. Note 10.0 is a *stricter* k-mer floor than the validated 2.0 — it could refute the genuine low-k-mer-cov partial 2b (ERR1810453, k-mer cov ~5 per the test fixture) that D-10 must confirm. **This is a real correctness risk for CLASS-02 and must be resolved before calibrating the score.** [VERIFIED: nextflow.config:46-47, conf/modules_hcv.config:328, bin/denovo_confirm.R:37]
**Warning signs:** ERR1810453's 2b candidate classifies as `background` instead of `co-infection` in a calibration run.

### Pitfall 2: Positional ext.args coupling — a new score_weight_* param shifts existing positions
**What goes wrong:** SUMMARIZE `ext.args` is a **space-joined positional string** (`conf/modules_hcv.config:328`); `summarize.R` reads args by index (`args[4]`..`args[11]`). The mapping is subtle: config passes `… minRead minCov n_candidates` (positions 9/10/11 after samplesheet/version/name), and `summarize.R` reads `args[9]`/`args[10]` as **`min_targeted_read`/`min_targeted_cov`** (fed from `minRead`/`minCov` — which is exactly D-09's "one threshold set"), and `args[11]` as `n_candidates`. Appending `score_weight_*` anywhere except the **end** silently re-maps every later index.
**Why it happens:** No named-arg parsing; pure positional coupling between Groovy string order and R `args[n]`.
**How to avoid:** Append new `score_weight_*` (and any evenness-transform constant) **after** `n_candidates` in both the `ext.args` string and the `summarize.R` `commandArgs` block, with defensive defaults (`if (length(args) >= N && nchar(args[N]) > 0) … else <default>`) exactly like the existing `denovo_*` / `min_targeted_*` parse. [VERIFIED: conf/modules_hcv.config:328, bin/summarize.R:29-45]
**Warning signs:** `n_candidates` parses as a weight, or weights parse as NA → score collapses to read-count-only.

### Pitfall 3: CV-evenness needs full-reference per-position depth incl. zeros — RESOLVED, but guard the zero-coverage edge case
**What goes wrong:** D-03's CV "penalizes both patchiness and low breadth in one number." That only holds if `cov$X3` contains a row for **every** reference position (zeros included). **This is satisfied:** SAMTOOLS_DEPTH runs with `-aa` (`conf/modules_hcv.config:259`), so every reference position is emitted. The residual edge case is a candidate with `mean(cov$X3) == 0` (no reads at all) → `sd/mean` = NaN/Inf.
**How to avoid:** Guard the zero-mean / zero-length case explicitly (`if (ref_length > 0 && mean(cov$X3) > 0) … else evenness factor = 0`), as in the Pattern 2 snippet. A no-coverage candidate gets evenness 0 (correctly scored to the bottom). [VERIFIED: conf/modules_hcv.config:259 `-aa`]
**Warning signs:** NaN/Inf dominance scores for no-coverage candidates; a candidate with one spike scoring oddly high.

### Pitfall 4: Snapshot / stub drift on Summary.csv schema change
**What goes wrong:** `modules/local/summarize/main.nf` hard-codes the full `Summary.csv` header in its **stub** block, and `tests/default.nf.test.snap` pins the output. Adding `overall_sample_call` / removing `minor_denovo_status` (D-15) breaks both.
**How to avoid:** Update the stub header in lockstep, then regen the nf-test snapshot. ⚠️ MEMORY notes a **SUMMARIZE snapshot regen is already deferred** (Phase 7, pending docker/disk) — fold Phase-8 schema changes into that same regen, and budget for the host-disk constraint (MEMORY "Host disk near-full"). The nf-test invocation is the pinned `~/.nf-test` 0.9.3 + `--profile docker` (MEMORY "nf-test local install").
**Warning signs:** `nf-test` fails on snapshot mismatch; CI red.

### Pitfall 5: "De novo worked for the dominant" measurement (D-11)
**What goes wrong:** The asymmetric refute hinges on "did de novo assemble a substantial contig for the **dominant**?" If measured against the wrong candidate's support row, the 4g refute fails (or a genuine minor gets wrongly refuted).
**How to avoid:** Reuse the dominant candidate's own joined assembly-support metrics (`cand_<dominant_rank>_assembly_support_best_contig_*`) against the same `denovo_*` floors. The dominant is known only *after* scoring (D-01), so the asymmetric-refute step must run **after** the dominant is chosen. Order: score → choose dominant → evaluate each non-dominant's corroboration **relative to the dominant's de-novo success**. [VERIFIED: D-11 + join structure]

## Code Examples

### Reading the per-candidate evidence the classifier consumes (already in summarize.R)
```r
# Source: bin/summarize.R:579 (Phase-7 join) + :620-646 (wide pivot) + :894 (left_join into final)
candidate_support <- join_assembly_support(candidates_long, support_df, denovo_match_level)
# long: one row/candidate with assembly_support ("supported"/"none") + best_contig_* metrics
# This long frame (NOT the wide pivot) is the natural classifier input — it keeps one row
# per candidate, which the role classifier and the enriched *.candidates.csv both need.
```

### Porting is_valid_minor() exception logic (D-12)
```r
# Source: bin/summarize_mapping_to_all_references.R (the 1a/1b allow + 2k1b block) + genotype_utils.R
# As a post-scoring special case in classify_roles():
same_geno <- genotype_from_subtype(cand_subtype) == genotype_from_subtype(dom_subtype)
is_1a1b   <- (cand_subtype %in% c("1a","1b")) && (dom_subtype %in% c("1a","1b")) && cand_subtype != dom_subtype
is_2k1b   <- "2k1b" %in% c(cand_subtype, dom_subtype)
# co-infection blocked (-> background) when: (same_geno & !is_1a1b) OR is_2k1b
```

### Evidence-table calibration assertion (the headline test, SCORE-02)
```r
# Source: handoff §2 evidence table — the false 4g must lose to a genuine even minor.
# Synthetic: 4g = 53279 reads, breadth 0.677, spiky depth (high CV) ; genuine 2b = fewer reads, even.
# Assert: dominance_score(genuine_2b) classification path != background-by-score, and the
# 4g candidate (no 4g contig, dominant 1a HAS a contig) -> role == "background", role_reason == "refuted_denovo".
```

## State of the Art

| Old Approach (pre-Phase-8) | Current Approach (Phase 8) | When Changed | Impact |
|----------------------------|----------------------------|--------------|--------|
| `minor_typable` YES/NO + downgrade-only `minor_denovo_status` (confirmed/refuted/unconfirmed/not_evaluated) | Per-candidate `role` (dominant/co-infection/background) + `role_reason` + dominance score | Phase 8 (this) | Replaces the minor-coupled 2-strain model with an N-candidate role model |
| Read-count-driven major/minor selection | Breadth-evenness-weighted dominance score | Phase 8 | Catches high-read spiky cross-mapping artefacts (the 4g) |
| Backgrounds silently dropped | Backgrounds surfaced in long `*.candidates.csv` with reason | Phase 8 | CLASS-03 |

**Deprecated/retired in this phase (D-15):**
- `classify_minor_denovo()` consumption (denovo_confirm.R) — verdict replaced by the corroboration step.
- `apply_denovo_layer()` consumption (denovo_layer.R) — downgrade layer replaced by role assignment.
- `minor_denovo_status`, `coinfection_flag` columns — replaced by `role` / `role_reason` / `overall_sample_call`.
- `review_flag` triggers rewired from `minor_denovo_status`/`coinfection_flag`/`gate_flag` onto roles.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | ~~`cov$X3` contains every reference position incl. zeros~~ — **VERIFIED, no longer an assumption**: SAMTOOLS_DEPTH runs `-aa` (conf/modules_hcv.config:259). CV-over-full-reference is correct. | Pattern 2 / Pitfall 3 | Resolved. |
| A2 | ERR1810453's genuine partial 2b has k-mer cov ≈ 5 (above 2.0, below 10.0), so the 10.0 nextflow.config floor would wrongly refute it. | Pitfall 1 | If its k-mer cov is actually ≥10, the config mismatch is benign for this sample (but still a latent correctness gap). Validate against the real contig numbers / 03-RESEARCH fixture. |
| A3 | Calibrating `score_weight_*` defaults so evenness out-weighs reads is sufficient to make the 4g lose to genuine even minors across the whole evidence table (not just the named samples). | D-04 / Calibration | Under-calibration → false co-infection survives (CLASS-02 fails); over-calibration → genuine low-yield minor demoted. The manuscript's before/after claim depends on this. |
| A4 | The validated denovo floors are 1000 / 2.0 / 90 and the planner will reconcile nextflow.config to them (or re-validate 500 / 10.0). | Pitfall 1 / D-10 | Wrong floor → wrong corroboration verdict for partial-contig genuine minors. |
| A5 | ~~`is_valid_minor()` still in the file~~ — **VERIFIED FALSE**: `is_valid_minor()` is **gone** from `bin/` (grep finds it only in `bin/tests/test_candidate_selection.R:24,183` as comments noting its Phase-6 removal). The D-12 "port verbatim" must recover the rule from **git history** (`git log -p -- bin/summarize_mapping_to_all_references.R`) or the Phase-6 SUMMARY, NOT a live file. | D-12 / Don't Hand-Roll / Open Q3 | Planner must add a recovery step; the canonical rule text is in history. The rule itself (different-genotype; allow 1a/1b; block 2k1b) is well-documented in D-12, handoff §5, and STATE — reconstructable even if the exact source is lost. |

## Open Questions

1. **Does the depth tsv include zero-depth positions?** (A1) — **RESOLVED**
   - SAMTOOLS_DEPTH runs with `-aa` (`conf/modules_hcv.config:259`). `cov$X3` includes every reference position incl. zeros. CV-evenness (D-03) is correct as designed. Only guard the zero-mean edge case.

2. **nextflow.config denovo floors: correct to 1000/2.0 or re-validate 500/10.0?** (Pitfall 1, A2, A4) — **OPEN, highest-priority**
   - What we know: 03-RESEARCH validated 1000/2.0/90; config ships 500/10.0/90; modules_hcv.config does not override.
   - What's unclear: whether shipping 500/10.0 was intentional or an un-synced default. 10.0 k-mer is *stricter* and could refute the genuine partial 2b.
   - Recommendation: resolve before score calibration — the corroboration verdict (D-10) is an *input* to the role classifier. Likely a one-line config correction + schema update; surface it as an explicit planned task.

3. **Where does the canonical `is_valid_minor()` rule live now?** (A5) — **RESOLVED: it is gone from `bin/`**
   - Confirmed removed in Phase 6; grep finds it only in test comments. Recover the verbatim rule from `git log -p -- bin/summarize_mapping_to_all_references.R` (or the Phase-6 SUMMARY). The rule (different-genotype required; allow 1a/1b; block 2k1b) is independently documented in D-12 + handoff §5, so it is reconstructable even without the exact source line.
   - Recommendation: add an explicit "recover is_valid_minor from history" step to the plan; do not assume a live file.

4. **One combined helper or two (`dominance_score.R` + `classify_roles.R`)?** (Claude's Discretion)
   - Recommendation: one file `classify_roles.R` exposing two functions (`score_candidates()`, `classify_roles()`) keeps staging simple (one new `path()` input) while keeping the score independently testable. Planner's call.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| R + tidyverse + seqinr | summarize.R + new helpers | ✓ (in SUMMARIZE container) | wave `3536dd50a17de0ab` | — |
| Rscript (local) | `bin/tests/run_all.sh` unit tests | ✓ assumed on dev host | — | run tests inside the container if absent |
| nf-test 0.9.3 | snapshot regen (Pitfall 4) | ✓ (`~/.nf-test`, MEMORY) | 0.9.3 pinned | invoke via PATH prefix + `--profile docker` |
| Docker + free disk | nf-test snapshot regen | ⚠️ disk near-full (MEMORY) | — | `sudo rm work/` + `docker volume prune` before regen |
| Manuscript mounts (`/mnt/N/.../HCVTyper_*_data/`) | real-contig calibration (D-04) for score weights | ⚠️ read-only mounts, may not be present this session | — | use the committed evidence-table numbers in handoff §2 + the 03-RESEARCH fixtures as the calibration source |

**Missing dependencies with no fallback:** None block code authoring. Score calibration (D-04) ideally uses the real per-contig/per-position numbers; if the mounts are unavailable, calibrate against the handoff §2 table values + `bin/tests/fixtures` and flag the defaults for an analyst validation pass (consistent with the milestone's "analyst validates after" out-of-scope note in REQUIREMENTS.md).

## Validation Architecture

> `workflow.nyquist_validation: true` — section included.

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Plain `Rscript` unit tests (`bin/tests/test_*.R`, source + `fail()`/`ok()`), NOT testthat; pipeline-level nf-test 0.9.3 |
| Config file | `bin/tests/run_all.sh` (R) ; `tests/default.nf.test` + `tests/nextflow.config` (nf-test) |
| Quick run command | `bash bin/tests/run_all.sh` |
| Full suite command | `bash bin/tests/run_all.sh` then (gated) `PATH=~/.nf-test:$PATH nf-test test tests/default.nf.test --profile docker` |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| SCORE-01 | Score combines reads + kmer_cov + breadth | unit | `Rscript bin/tests/test_dominance_score.R` | ❌ Wave 0 |
| SCORE-02 | Evenness out-weighs reads (4g loses to even minor) | unit | `Rscript bin/tests/test_dominance_score.R` | ❌ Wave 0 |
| CLASS-01 | Each candidate gets exactly one of 3 roles | unit | `Rscript bin/tests/test_classify_roles.R` | ❌ Wave 0 |
| CLASS-02 | co-infection = floor AND support; asymmetric refute → background | unit | `Rscript bin/tests/test_classify_roles.R` | ❌ Wave 0 |
| CLASS-03 | background surfaced with role_reason, no row loss | unit | `Rscript bin/tests/test_classify_roles.R` | ❌ Wave 0 |
| CLASS-04 | sample call derived from roles (3-value) | unit | `Rscript bin/tests/test_classify_roles.R` | ❌ Wave 0 |
| (integration) | Summary.csv schema + no crash on no-mapping/skip-assembly | nf-test snapshot | `nf-test test tests/default.nf.test --profile docker` | ⚠️ exists, snapshot regen deferred |

### Sampling Rate
- **Per task commit:** `bash bin/tests/run_all.sh` (R units, < 30s, no docker).
- **Per wave merge:** `bash bin/tests/run_all.sh` (all R tests).
- **Phase gate:** R suite green + (gated, when docker/disk available) nf-test snapshot regenerated and green. Full golden-baseline reproduction is **Phase 9 / TEST-01**, not a Phase-8 gate (per CONTEXT scope).

### Wave 0 Gaps
- [ ] `bin/tests/test_dominance_score.R` — covers SCORE-01/02 (monotonicity + evenness-dominates the 4g case)
- [ ] `bin/tests/test_classify_roles.R` — covers CLASS-01..04 (all roles, asymmetric refute D-11, exceptions D-12, sample call D-14)
- [ ] `bin/tests/fixtures/` evidence-table rows — 4g refute, ERR1810447/ERR1810453 2b confirm, sim1 1a:1b + sim2 2a:3a + ERR1810505/ERR1810511 preserved, IVT `uncorroborated_kept`
- [ ] No framework install needed (Rscript + tidyverse already present; `run_all.sh` auto-globs new tests)
- [ ] (Deferred-but-track) SUMMARIZE stub header update + nf-test snapshot regen — coordinate with the already-pending Phase-7 regen

## Security Domain

> `security_enforcement: true`, ASVS level 1. This phase processes **trusted internal data** (pipeline-generated CSV/TSV from earlier stages, params from tracked config). There is no auth, session, network, or external-user-input surface.

### Applicable ASVS Categories
| ASVS Category | Applies | Standard Control |
|---------------|---------|-----------------|
| V2 Authentication | no | — (no auth surface) |
| V3 Session Management | no | — |
| V4 Access Control | no | — |
| V5 Input Validation | yes (low) | `as.numeric`/`as.integer` coercion of `score_weight_*` args yields NA on garbage → defensive default fallback (existing pattern, summarize.R:25-45); `stopifnot(match_level %in% ...)` guard (assembly_support_join.R:50). Pin `read_csv` `col_types` for new columns (CR-01/CR-02 precedent) to avoid type-inference aborts. |
| V6 Cryptography | no | — (no secrets/crypto) |

### Known Threat Patterns for R-in-Nextflow summary step
| Pattern | STRIDE | Standard Mitigation |
|---------|--------|---------------------|
| Malformed/empty upstream CSV (no-mapping, skip-assembly) aborts SUMMARIZE (DoS) | Denial of Service | Typed-empty-tibble guards (generalize PLUMB-02 / T-03-01 to the new score/role columns) — NA-fill, never `stop()` |
| Positional `ext.args` misparse silently changes scoring behaviour | Tampering (config) | Append new args at the end + defensive index-guarded parse (Pitfall 2) |
| Numeric-looking genotype key (`"1"`→double) aborts the join | Denial of Service | Pin `col_types` / `as.character` coercion (already done in assembly_support_join.R:90) |

No high/critical security findings — this is internal trusted-data processing. `security_block_on: high` is satisfied.

## Sources

### Primary (HIGH confidence — read in-session)
- `bin/summarize.R` (full; key regions L1-120, L330-441, L509-646, L860-1005, L1075-1188) — the primary edit site, confirmation layer, review_flag builder, cov loop, joins.
- `bin/summarize_mapping_to_all_references.R` — neutral candidate selection + (historical) `is_valid_minor()` site.
- `bin/denovo_confirm.R`, `bin/denovo_layer.R` — the legacy path to retire (D-15) + the calibrated floors to reuse (D-10).
- `bin/assembly_support_join.R` + `bin/tests/test_assembly_support_join.R` — Phase-7 join (classifier input) + the house unit-test pattern.
- `bin/genotype_utils.R` — `genotype_from_subtype()` (2k1b-aware).
- `bin/tests/run_all.sh`, `bin/tests/test_coinfection.R`, `bin/tests/test_summarize_denovo.R`, `bin/tests/fixtures/` — test harness + fixture conventions.
- `conf/modules_hcv.config:13-22,119-129,327-336`, `nextflow.config:46-53`, `nextflow_schema.json:191-225`, `modules/local/summarize/main.nf` — param wiring + staging + stub schema.
- `hcvtyper_handoff_denovo_informed_selection.md` §2 (evidence table), §4 (asymmetric rule), §5 (substantial-contig thresholds) — the calibration anchor.
- `.planning/phases/0{3,6,7}-*/0*-CONTEXT.md` + 03-RESEARCH.md — upstream decisions + validated denovo floors.
- `.planning/REQUIREMENTS.md`, `.planning/STATE.md`, `.planning/phases/08-*/08-CONTEXT.md` — scope + decisions.

### Secondary (MEDIUM)
- MEMORY notes: SUMMARIZE snapshot regen deferred; nf-test local install; host disk near-full; SUMMARIZE channel arity.

### Tertiary (LOW)
- None — no WebSearch needed; this is an internal-codebase logic phase.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — no new deps; container verified.
- Architecture / data flow: HIGH — every input traced through the real code (candidates → join → cov loop → classifier → outputs).
- Decisions buildability: HIGH — all D-01..D-16 confirmed implementable as written.
- Calibration (score weights): MEDIUM — depends on real per-contig/per-position numbers (manuscript mounts) which may be unavailable this session; handoff §2 table + 03 fixtures are the fallback anchor (A3).
- Pitfalls: HIGH — the calibration-default mismatch and positional-arg coupling were verified directly in config/schema/R.

**Research date:** 2026-06-13
**Valid until:** 2026-07-13 (stable internal codebase; re-check only if `summarize.R` / the config arg order changes)
