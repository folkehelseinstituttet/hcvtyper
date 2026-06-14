# Phase 7: Per-Genotype Assembly Support - Context

**Gathered:** 2026-06-12
**Status:** Ready for planning

<domain>
## Phase Boundary

Reframe de novo/BLAST evidence from a **dominance-coupled, minor-only confirmation** into an **independent per-genotype "assembly support" signal**, then join that signal to the Phase 6 long-format candidate set at genotype level.

Today (`bin/denovo_confirm.R::classify_minor_denovo` + `bin/blast_parse.R` §7) the evidence is computed *relative to the major/minor call* — `blast_parse.R` picks a `minor_ref` as "best hit of a different genotype than major," and `classify_minor_denovo()` returns confirmed/refuted/unconfirmed *for the reported minor*. Phase 7 makes the evidence neutral: for every (sub)type seen in the de novo contigs, summarize best contig length, BLAST % identity, BLAST alignment length, and k-mer coverage — computed with no reference to which candidate is dominant — and left-join it to candidates at a parameterized match level (default genotype).

**Covers requirements:** ASUP-01 (per-genotype assembly support, mapping-independent), ASUP-02 (genotype-level join to candidates, parameterized match level).

**Explicitly NOT in this phase:**
- The "substantial enough to corroborate?" threshold verdict and dominant/co-infection/background role assignment → **Phase 8** (CLASS-02/SCORE). Phase 7 emits raw signal; Phase 8 applies the calibration floors that already exist as params.
- Retiring the legacy minor-coupled confirmation path (`classify_minor_denovo` / `minor_denovo_status` / `review_flag`) → kept working as a shim this phase (D-04); consumed/retired in Phase 8.
- Filename slot migration `.major.`/`.minor.` → `.cand1.`/`.cand2.` and the `summarize.R` parsing cutover → **Phase 9** (COMPAT-02/03).

</domain>

<decisions>
## Implementation Decisions

### Signal vs verdict — Phase 7/8 boundary
- **D-01:** The per-genotype support table emits **raw metrics only** — best contig length, BLAST % identity, BLAST alignment length, k-mer coverage. The "substantial enough to corroborate?" threshold verdict (the `sc_length ≥ min_len & kmer_cov ≥ min_kmer & pident ≥ min_pid` floor) is **NOT applied in Phase 7** — it stays a Phase 8 (CLASS-02) concern, reusing the existing `denovo_min_contig_length` / `denovo_min_kmer_cov` / `denovo_min_blast_identity` params.
- **Why:** Mirrors the Phase 6 D-05 split — the mechanical signal is computed neutrally; interpretation ("is this support real?") lives at the summary/classification step. Keeps Phase 7 additive and threshold-agnostic.
- **No-support representation (criterion #3):** A candidate with no genotype-level evidence resolves to an **explicit status** (e.g. `assembly_support = "none"`) with **NA metrics**, via `left_join` + NA-fill — **no row loss**. Consistent with the v1.0 `left_join` plumbing and the PLUMB-02 typed-empty-tibble guard already in `summarize.R`.

### Table grain — taxonomic granularity of the support table
- **D-02:** `bin/blast_parse.R` emits the support table at **subtype grain** (one row per subtype seen in the contigs), each row carrying a **genotype key** derived via `genotype_from_subtype()` (2k1b-aware). `bin/summarize.R` **collapses to the level set by `denovo_match_level`** at join time (default `genotype`).
- **Why:** Keeps `denovo_match_level = "subtype"` honorable without recomputing the table — the finer grain is always available, collapse happens at the join. Criterion #1's "for every genotype" is satisfied at the default match level (subtype rows roll up to genotype). Chosen over pre-aggregating to genotype in `blast_parse.R`, which would strand the subtype match level.

### Best-contig definition — how the four metrics are chosen
- **D-03:** Within each (sub)type group, pick the **single best contig by contig length (`sc_length`)** and report **that contig's** identity, alignment length, and k-mer coverage — one coherent contig backs all four numbers. NOT independent per-metric maxima (which could mix unrelated noise contigs), NOT best-by-bitscore.
- **Why:** Matches the milestone's "substantial contig" framing (confirm ERR1810453's 2,949 bp 2b while ignoring ~300 bp noise) and is the most direct route to reproducing today's length-based `minor_contig_length` for criterion #4.

### Additive vs cutover — migration strategy
- **D-04:** Phase 7 is **additive with a legacy shim** (mirrors Phase 6 D-06). It ADDS the neutral per-genotype support table + candidate join, while keeping `classify_minor_denovo()` / `minor_denovo_status` / `review_flag` working **unchanged**, so Phase 7 ships non-breaking on its own. The legacy minor-coupled path is consumed/retired when Phase 8 classification lands on the new join.
- **Why:** Every intermediate phase stays shippable and runnable. Rewiring `review_flag` now (full cutover) would change reporting behavior before Phase 8's role logic exists, risking regression with nothing to validate against.
- **CONSEQUENCE for planner + verifier:** Phase 7 success-criterion #4 ("the assembly-support join reproduces the de novo evidence currently attached to the minor slot, default flags, N=2") is checked against the **new join's output for `cand_2`**, which must equal today's `denovo_minor_*` for the regression fixtures. Because the legacy path still runs in parallel (D-04), both the old `denovo_minor_*` columns and the new per-candidate support columns appear in `Summary.csv` this phase — that duplication is intentional and removed in Phase 8/9.

### Claude's Discretion
- **Join site:** `bin/summarize.R`, reusing the existing `df_denovo` / `df_blast_out` left-join plumbing (already keyed by `sampleName`, already sources `genotype_from_subtype()` and `denovo_confirm.R`). No new Nextflow module expected.
- **Candidate-side genotype key:** derive from `candidate_ref` (the Phase 6 long-format candidate column) via `genotype_from_subtype()` — same helper used on the support-table side, so the join keys are computed identically on both sides.
- Exact column names/types of the new support table and the per-candidate support columns in `Summary.csv` (subject to the legacy-column-retention constraint in D-04). Follow the existing `denovo_*` naming precedent.
- Whether the support table is a new CSV emitted by `BLASTPARSE` (e.g. `*.assembly_support.csv`) or folded into the existing `*.blastparse.csv` / `*_blast_out.csv` outputs — planner/researcher to pick the least-disruptive shape given `summarize.R` already reads `*_blast_out.csv` into `df_blast_out`.
- Preserve existing safe defaults on empty/skip-assembly runs (typed-empty tibble → NA support columns, never abort — the T-03-01 DoS guard).

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Phase scope & requirements
- `.planning/ROADMAP.md` §"Phase 7: Per-Genotype Assembly Support" — goal, the four success criteria (note criterion #4 interpretation in D-04), dependency on Phase 6.
- `.planning/REQUIREMENTS.md` — ASUP-01, ASUP-02 (Phase 7); plus CLASS-02/SCORE (Phase 8) and COMPAT (Phase 9) for cross-phase awareness of what Phase 7 must NOT do.
- `.planning/PROJECT.md` — Core Value, the "substantial contig" calibration constraint (≥~1,000 bp / identity / k-mer floor), genotype-level match constraint, non-breaking constraint.
- `.planning/phases/06-neutral-candidate-selection/06-CONTEXT.md` — D-01 (long-format candidate table = the join target), D-06 (additive + legacy shim pattern this phase mirrors).

### Code touchpoints (the files this phase changes)
- `bin/blast_parse.R` — currently §7 picks `major_ref`/`minor_ref` coupled to dominance and writes `*.blastparse.csv` (`sample, major_ref, major_contig_length, minor_ref, minor_contig_length`); §6 filters contigs ≥500 bp. Primary site for the neutral per-subtype support table (D-01..D-03). The legacy `*.blastparse.csv` shape is retained this phase (D-04).
- `bin/denovo_confirm.R` — `classify_minor_denovo(blast_out_df, major_geno, minor_geno, min_len, min_kmer, min_pid, match_level)`: the current minor-coupled confirmed/refuted/unconfirmed verdict. Kept working as a shim (D-04); its threshold floors are the calibration that moves to Phase 8.
- `bin/summarize.R` — reads `*.blastparse.csv` → `df_denovo` (~L441-465) and `*_blast_out.csv` → `df_blast_out` (~L476-492); runs the confirmation layer (~L797-815) and builds `review_flag` (~L938+). Site of the new genotype-level candidate join + NA-fill (D-01, D-02) and the candidate-side genotype key (Discretion).
- `bin/genotype_utils.R` — canonical `genotype_from_subtype()` (2k1b-aware); used on BOTH the support-table side and the candidate side so join keys match.
- `modules/local/blastparse/main.nf` — BLASTPARSE process; declared outputs (`blast_res`, `csv`, `contigs`, `major_fasta`/`minor_fasta`, `png`) + stub contract. If a new support CSV is added, declare it here and in the stub.
- `conf/modules_hcv.config` (~L328) — passes `denovo_min_contig_length denovo_min_kmer_cov denovo_min_blast_identity denovo_match_level denovo_confirm_minor minRead minCov` as args to `summarize.R`. Params already exist; Phase 7 reuses, does not add.

### Conventions / patterns
- `.planning/codebase/CONVENTIONS.md`, `.planning/codebase/ARCHITECTURE.md` — R helper conventions, R-emits/Nextflow-routes pattern, channel/join discipline.

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `genotype_from_subtype()` (`bin/genotype_utils.R`) — already sourced in `summarize.R` and staged into workdirs; 2k1b-aware. Reuse for BOTH the support-table genotype key and the candidate-side key (D-02, Discretion).
- The de novo confirmation params (`denovo_min_contig_length` 1000, `denovo_min_kmer_cov` 2.0, `denovo_min_blast_identity` 90, `denovo_match_level` "genotype") and the threshold-floor logic in `classify_minor_denovo()` — already calibration-validated (03-RESEARCH). Phase 7 leaves them in place; Phase 8 reuses the floors as the corroboration verdict.
- `summarize.R`'s typed-empty-tibble guards for `df_denovo` / `df_blast_out` (PLUMB-02 / T-03-01) — generalize to the new support columns so skip-assembly runs NA-fill rather than abort.
- `blast_parse.R` already extracts `sc_length` (contig length, from `NODE_..._length_<len>`) and `kmer_cov` (`..._cov_<cov>`) and `subtype` per BLAST hit, plus `scaf_top` (best hit per scaffold). The per-subtype support roll-up builds directly on `scaf_top`.

### Established Patterns
- **R-emits, Nextflow-routes:** selection/evidence logic emits CSV columns; Nextflow routes on them. The support table continues this — `blast_parse.R` emits per-subtype support rows, `summarize.R` consumes and joins.
- **`left_join` + typed-empty fallback** for de novo evidence in `summarize.R` — the exact plumbing the new candidate join reuses (criterion #3 NA-fill, no row loss).
- **Additive + legacy-column retention** (Phase 6 D-06) — keep legacy `denovo_minor_*` columns alongside the new per-candidate support columns until Phase 8/9.

### Integration Points
- BLASTPARSE → `*.blastparse.csv` / `*_blast_out.csv` → `summarize.R` (`df_denovo` / `df_blast_out`) → confirmation layer + `review_flag` → `Summary.csv`. Phase 7 inserts the per-genotype support table into this read path and adds the candidate join; the long-format candidate rows (Phase 6) are the new join's left side.
- The join keys candidates (`candidate_ref` → genotype) to support rows (subtype → genotype) at `denovo_match_level` — both sides via `genotype_from_subtype()`.

</code_context>

<specifics>
## Specific Ideas

- Four support metrics are fixed by ASUP-01 / criterion #1: **best contig length, BLAST % identity, BLAST alignment length, k-mer coverage** — note BLAST alignment length (the `length` column) is reported in addition to full contig length (`sc_length`); they are distinct and both required.
- Default `denovo_match_level = "genotype"` — a 3b candidate is corroborated by 3a assembly support (criterion #2, the established genotype-match rule).
- Criterion #4 is the reproduction anchor: at default flags + N=2, the new join's `cand_2` support must equal today's `denovo_minor_*` evidence on the regression fixtures.

</specifics>

<deferred>
## Deferred Ideas

- Applying the substantiality threshold floors as the corroboration verdict, the "guilty until corroborated" stance, and dominant/co-infection/background role assignment — **Phase 8** (CLASS-02/SCORE).
- Retiring `classify_minor_denovo` / `minor_denovo_status` / the legacy `review_flag` minor-coupled path once the new join feeds classification — **Phase 8**.
- Filename slot migration `.major.`/`.minor.` → `.cand1.`/`.cand2.`, `summarize.R` filename parsing cutover, legacy `Major_*`/`Minor_*` column aliasing — **Phase 9** (COMPAT-02/03).

None — discussion stayed within phase scope.

</deferred>

---

*Phase: 07-per-genotype-assembly-support*
*Context gathered: 2026-06-12*
