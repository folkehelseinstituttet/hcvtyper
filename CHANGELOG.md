# folkehelseinstituttet/hcvtyper: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### `Added`

### `Changed`

### `Fixed`

### `Dependencies`

### `Deprecated`

## 2.0.0 - 2026.09.11

Major release. The major/minor strain model is replaced end to end by a neutral, evidence-first candidate model: references are ranked without dominance semantics, each candidate is corroborated against an independent _de novo_ assembly, and roles (`dominant` / `co-infection` / `background`) are assigned only at the reporting step. Output filenames and the `Summary.csv` schema change accordingly — see `Breaking changes`.

**Upgrading:** the last official release was **v1.1.7**. The `v1.2.0` and `1.3.0` tags exist in the repository but were never published as releases, so their changes are listed under their own headings below and are also included here — anyone coming from v1.1.7 receives all three sets of changes at once.

### `Breaking changes`

- **Per-candidate output filenames changed: `.major.` / `.minor.` → `.cand1.` / `.cand2.`** across every per-candidate output (BAM, mapping stats, consensus, depth, variation plots). Any downstream script globbing `*.major.*` or `*.minor.*` must be updated. `bin/summarize.R` recovers `candidate_rank` by joining against `*.candidates.csv` rather than parsing a hard-coded filename field.
- **`Summary.csv` schema changed.** New columns: `overall_sample_call`, the `Major_role_*` / `Minor_role_*` family, `Major_evidence_state` / `Minor_evidence_state`, the per-candidate `cand_N_*` metrics, and the evidence-contribution columns. `Major_genotype` / `Minor_genotype` are now derived from the role-corrected subtypes, which **changes reported values** on samples whose roles are reversed relative to first mapping (see `Fixed`). Consumers should key on column names rather than positions.
- **`--skip_assembly` no longer exists.** It was deleted from the code in `b074dbd` — where it had been silently accepted while doing nothing — and is now removed from `nextflow_schema.json`. Assembly always runs.
- **`evidence_state` is three-valued** (`confirmed` / `probable` / `weak`). The former fourth state `refuted` was structurally unreachable and has been removed — see `Removed`.
- **Release tags are bare version numbers from 2.0.0 onward** (`2.0.0`, not `v2.0.0`), matching `manifest.version` and nf-core convention; the `v` belongs only in the displayed release name. Existing tags are untouched — `v1.0`–`v1.2.0` keep their prefix and `-r v1.1.7` keeps working, while `1.3.0` was already tagged bare before the convention was written down. New releases are pinned as `-r 2.0.0`.

### `Added`

- **Competitive joint mapping (`JOINT_MAPPING`), replacing independent per-candidate mapping.** Reads are mapped **once** against a combined index of all selected candidate references, deduplicated on the combined BAM, then split per candidate — so each read is assigned to the reference it fits best instead of being counted once per candidate. This is what makes cross-mapping artefacts visible to the dominance score rather than inflating every candidate equally. Emits per-candidate deduplicated BAMs (meta carrying `candidate_nodup_reads`), combined pre- and post-dedup idxstats, per-candidate depth, consensus and variation plots.
- **_De novo_ subtype rescue (`RESCUE_EVALUATION`).** Where a high-quality contig contradicts a mapped candidate's subtype, the candidate's reference is **reassigned before mapping** rather than surfacing as an unresolved conflict afterwards. Triggering thresholds: `--rescue_min_length` (3000 bp), `--rescue_min_pident` (85), `--rescue_min_aln_length` (3000 bp), `--rescue_min_kmer_cov` (2.0), plus `--rescue_1a1b_length` (5000 bp), a stricter contig-length floor when both the candidate and the de novo subtype fall in {1a, 1b}. Two guards can block a replacement: `--rescue_kmer_cov_ratio` (10) blocks when the candidate's own-subtype contig exceeds the rescue target's k-mer coverage by that factor (cross-mapping-noise guard), and `--rescue_dominant_protect_cov` (90) protects a candidate at or above that coverage whose own subtype is itself assembled. Every rescue, nomination and block decision is written to a per-sample `*.rescue_audit.csv`, so no reference reassignment is silent.
- **Consensus distance (`CONSENSUS_DISTANCE`).** A new `*.consensus_distance.tsv` per sample, reporting similarity percentage and difference count between candidate consensus sequences from a Needleman–Wunsch global alignment (`pwalign`). Gap columns count as differences; `N` and zero-coverage columns are excluded from the denominator. The separation between two called strains is now measured rather than assumed.
- **Neutral candidate selection (`--n_candidates`, default 2):** Reference selection ranks up to N candidates neutrally by read recruitment — no major/minor dominance semantics during the run. The top reference per distinct subtype is selected; `is_valid_minor()` validity filtering is removed from selection and moved to classification. A new `*.candidates.csv` (long-format, one row per candidate) is emitted by PARSEFIRSTMAPPING alongside the legacy wide CSV.
- **Per-genotype assembly support (`*.assembly_support.csv`):** `blast_parse.R` emits a per-subtype assembly-support roll-up (best contig by `sc_length`) carrying four metrics: best contig length, BLAST % identity, BLAST alignment length, and k-mer coverage. These are joined to the neutral candidates at genotype level (parameterised by `--denovo_match_level`, default `genotype`) via the new `bin/assembly_support_join.R` helper, replacing the former major/minor-specific de novo evidence pass-through. Candidates with no matching assembly evidence resolve to `assembly_support = "none"` with NA metrics (no row loss).
- **Dominance scoring:** Each candidate receives a combined dominance score over log10(mapped reads), coverage breadth fraction, CV-of-depth evenness, and log10(1 + k-mer coverage). Breadth evenness is weighted 3× raw read count by default (`--score_weight_evenness 3.0`, `--score_weight_reads 1.0`, `--score_weight_kmercov 0.5`), making the score resistant to index-hopping artefacts that have high read count but uneven breadth.
- **Strain-role classification (`dominant` / `co-infection` / `background`):** The new `bin/classify_roles.R` helper classifies each candidate at the summary step using a "guilty until corroborated" rule: a non-dominant candidate is reported as co-infection only when it (a) clears the abundance floor (`minRead` + `minCov`) **and** (b) has genotype-level assembly support; otherwise it is classified `background`. Background/artefact candidates are surfaced explicitly in `candidates.csv` with a `role_reason` (e.g. `no_own_assembly`, `weak_own_assembly_below_floor`, `below_floor`) and never silently dropped.
- **`overall_sample_call` column in `Summary.csv`:** Derived from candidate roles — `monoinfection`, `co-infection`, or `indeterminate` — reported once per sample.
- **New `candidates.csv` output:** Per-sample long-format file listing every candidate (including background) with `candidate_rank`, `dominance_score`, `role`, `role_reason`, and all assembly-support metrics. Published to `summary/candidates/`.
- **New scoring parameters:** `--score_weight_evenness`, `--score_weight_reads`, `--score_weight_kmercov`, `--score_evenness_k`.
- **`--review_min_offgenotype_contig_length` (default 1000 bp):** Contig-length floor for the monoinfection "de novo assembly found a different-genotype contig" review sentence. Deliberately separate from — and higher than — `--denovo_min_contig_length` (500 bp): that floor confirms a minor that mapping already supports, whereas this trigger fires where the contig is the _only_ evidence. See `Changed` below.
- **Per-candidate evidence report** in `candidates.csv`/`Summary.csv`: `evidence_summary` (a human-readable contig-corroboration sentence per candidate) plus three per-metric contribution columns (`contig_identity_contribution`, `contig_length_contribution`, `contig_kmer_contribution`), so a reader can reconstruct why a candidate landed in its evidence state and what drove its `assembly_support_score` without re-running the pipeline.
- **`Summary.csv` regression differ (`bin/tests/compare_summary_regression.R`):** Old-vs-new `Summary.csv` comparison that fails, with named examples, if any column outside an allowed set changed. Reads both files as all-character so a ragged schema or differing type inference cannot abort the comparison, and aligns on `sampleName` so row order cannot register as a difference. Not part of `run_all.sh` (it takes two file paths).
- **Threshold sweep tool (`bin/tests/offgeno_flag_sweep.R`):** Re-derives the off-genotype-contig review floor on a cohort of result directories, sweeping both candidate length definitions against the 2k1b exclusion and reporting must-keep retention and the `provisional`→`high` impact per cell. Aborts on a cohort-size mismatch or missing must-keep samples rather than emitting a plausible-looking table.
- **Regression suite extension (`bin/tests/test_compat.R`):** New auto-discovered test file covering golden strain-call reproduction (COMPAT-01), cand-slot filename lockstep (COMPAT-02), legacy + role column co-presence (COMPAT-03), and the 1a/1b co-infection and 2k/1b recombinant suppression exceptions (COMPAT-04). Runs automatically with `bash bin/tests/run_all.sh` (CI `r-regression` job, no YAML change required).
- **CI now asserts pipeline output, not just exit status.** A dedicated `nf-test` job runs the minimal profile and checks the published file tree and stable-content md5s against `tests/default.nf.test.snap`; the snapshot had drifted 92 commits before this job existed.
- **User documentation refreshed to the current classification model.** `README.md` gains a "How the pipeline calls strains" walkthrough of the eight-stage analysis, parameter tables for candidate selection, de novo corroboration, the rescue family, dominance-score weights, the review-flag floor and the contamination check, a rewritten `Summary.csv` column reference, and documentation of `summary/candidates.csv` and the per-sample evidence files under `parsefirstmapping/` and `blastparse/`. `docs/output_interpretation.md` gains the evidence-annotated review-sentence behaviour, the off-genotype-contig substantiality floor and its rationale, the `denovo_minor_*` contig-coherence note, and a threshold-tuning index. Deprecated columns are now labelled as such rather than described as live.
- **Three shipped-but-undocumented parameters declared in `nextflow_schema.json`:** `--review_min_offgenotype_contig_length`, `--rescue_kmer_cov_ratio` and `--rescue_dominant_protect_cov`. All three were live in `nextflow.config` and consumed by the R helpers, but absent from the schema — so they were missing from `--help`, from the launch interfaces built from the schema, and from any parameter validation.
- **`cand_cov_breadth` column in `summary/candidates.csv`** — the deduplicated targeted breadth@>=5x per candidate, the value the dominance score's breadth term now reads. `candidate_cov` is unchanged and still reports the first-pass all-reference breadth.
- **`bin/tests/test_dominance_breadth_source.R`** (DBS-1..9) — covers the per-row source preference, the units trap, rescued-candidate scoring and eligibility, the deliberate non-change that keeps eligibility from demoting anyone, and that `untypable` stays reachable for genuinely uncovered samples. Verified to fail against the pre-fix code.
- **A breadth-source wiring guard in `tests/default.nf.test`** asserting `summary/candidates.csv` carries `cand_cov_breadth`. No function-level test can catch "nothing populates the column the score prefers" — the resolver stays green either way — so this assertion is the only thing standing between a silent revert and a score that quietly returns to first-pass breadth.

### `Changed`

- **Continuous assembly-support scoring replaces the binary `own_substantial` threshold:** each candidate's own de novo/BLAST assembly support is now expressed as a continuous `assembly_support_score` (weighted length + identity + k-mer coverage, with a logistic identity term centered at 88%) instead of an ANDed length ≥ 500 bp / k-mer ≥ 2 / identity ≥ 90% pass/fail gate. A candidate one point below the old 90% identity cliff with a near-full-length, concordant contig is no longer silently discarded to `background`.
- **Per-candidate evidence state (`confirmed` / `probable` / `weak`)** replaces the binary per-candidate `own_substantial` flag. A candidate's state is computed entirely from its own score and concordance, and is never forced by another candidate's dominance.
- **`overall_sample_call` is derived from the count of `confirmed`/`probable` candidates**, independent of dominance ordering — a strong non-dominant candidate is reported as co-infection rather than gated out by which candidate happens to be dominant.
- **`Major_role_*` / `Minor_role_*` columns added to `Summary.csv`** alongside the new overall call: `Major_role_reference`, `Major_role_subtype`, `Major_role_dominance_score`, `Minor_role_reference`, `Minor_role_subtype`, `Minor_role_dominance_score`, `overall_sample_call`.
- The legacy `apply_denovo_layer` / `minor_denovo_status` / `coinfection_flag` classification path in `bin/summarize.R` is retired; `bin/denovo_layer.R` is still staged and unit-tested but is no longer called from the main reporting path. `review_flag` is now set from the per-sample role roll-up.
- **The "Major subtype conflict" review sentence now carries the measured evidence.** This trigger sets `call_confidence = "review"` on its own — the hardest sample-level signal — and named two subtypes with no numbers behind either. The conflicting contig's metrics were not obtainable from `Major_best_contig_*`, which describes the _candidate's_ genotype group, while the contig causing the disagreement belongs to `denovo_major_subtype`'s group. A short-anchor conflict now reads _"…mapping (3a) vs contig (1a) — 1816 bp contig, 120 bp aligned (7%), 89.0% identity, k-mer cov 1.3; only 7% of the contig aligns to any reference in the panel, so the contig subtype is weakly supported and the conflict may be an artefact of a short anchor rather than a real discrepancy."_, while a full-length conflict keeps the original "possible reference mismatch or highly divergent strain" reading. A 69 bp anchor cannot support a subtype call, so it cannot support a _disagreement_ with one either.
- **The "different-genotype contig" review sentence now carries the measured evidence.** It previously named a subtype and nothing else, then asked for a manual review — while the aligned length, identity and k-mer coverage of that contig existed only in `blastparse/*.assembly_support.csv`, never reaching `Summary.csv` (assembly support is joined at candidate grain, and an off-genotype contig is not a candidate). The sentence now reads e.g. _"…different-genotype contig (6i) — 1620 bp contig, 69 bp aligned (4%), 91.3% identity, k-mer cov 1.0; only 4% of the contig aligns to any reference in the panel, so the subtype assignment is weakly supported — the contig may be largely non-HCV, chimeric, or too divergent to type."_, while a contig aligning over its full length keeps the co-infection wording. The **aligned fraction**, not the contig length, is what separates the two: genuine minors align over 99–100% of their contigs. No flag is suppressed by this — the analyst decides, with the numbers present.
- **The monoinfection "different-genotype contig" review sentence is now gated on contig length and the 2k/1b pair rule.** It previously applied no substantiality floor at all: across five routine runs it fired on 72 of 140 samples (51%), was the only sentence on every one of them, and so accounted for 92% of the cohort's `provisional` calls — median triggering contig 606 bp, shortest 142 bp, i.e. below even `--denovo_min_contig_length`. Two gates now apply: the contig must reach `--review_min_offgenotype_contig_length` (1000 bp), and the pair must not be a 2k/1b recombinant against a genotype 1 or 2 major (the same `is_valid_minor()` rule 2 already applied to candidate promotion, which the review trigger never consulted). Length is the only leg used: the three legacy-typable minors this build demotes to monoinfection sit at k-mer coverage 1.42–1.97, _below_ `--denovo_min_kmer_cov`, so a k-mer leg would suppress exactly the samples most worth reviewing. Policy lives in a new pure helper `offgenotype_contig_reviewable()`, replacing a hand-rolled `substr(x, 1, 1)` comparison. **Affects `review_flag` and `call_confidence` only** — no subtype, typability, `overall_sample_call` or resistance value changes.
- **`review_flag` messages now name the specific candidate** they refer to (rank/slot, reference, or subtype) and carry the concrete measured value alongside the floor or expectation it missed (e.g. "candidate 2 contig matched 2c but mapping says 3a"; "identity 89.0 below 90 floor"), rather than only naming the trigger — including flags raised by a non-dominant candidate in an otherwise-monoinfection sample.

### `Fixed`

- **The dominance score's breadth term read the first-pass mapping, not the targeted one** (`260810-dbs`). `score_candidates()` prefers a `cand_cov_breadth` column and falls back to `candidate_cov` — but **no module has ever emitted `cand_cov_breadth`**, so the fallback was the only live path, and `candidate_cov` is the first-pass all-reference breadth while every other term in the same score is a second-pass quantity. This contradicts decision **D-08** ("the floor and the breadth-evenness score read each candidate's targeted (second) mapping coverage"), and the score's own calibration fixtures were built from targeted breadths. `bin/summarize.R`'s coverage loop now carries the deduplicated targeted breadth@>=5x (`cov_breadth_min_5`) as `cand_cov_breadth`, and the score resolves the two axes **per row** — targeted where it exists, first-pass where the candidate was never targeted-mapped, never 0. Worked example (ERR1810469, Thomson 2016): the 3a candidate was scored on 46% first-pass breadth against a real targeted breadth of 18.17%, and its 1a partner on 69% against 96.02%.
- **Rescued and nominated candidates scored as if none of their reference were covered.** `rescue_evaluation.R` correctly blanks `candidate_reads`/`candidate_cov` after a reference replacement (those numbers described the _displaced_ reference), but nothing recomputed breadth for the score, so the breadth term silently fell to 0 — a systematic 3.0-point penalty on exactly the candidates the de novo layer exists to surface. In the Thomson 2016 run this hit all five de-novo-surfaced co-infections, whose real targeted breadths are 62–99%.
- **A rescued or nominated candidate could never be selected as dominant, and a sole one was reported `untypable`.** The dominance eligibility test read the same blanked `candidate_cov`, so `NA` was treated as "no coverage" rather than "not measured on this axis". A sample whose only candidate had been rescued was reported with `overall_sample_call = untypable` — documented as "no usable coverage on any candidate" — while simultaneously carrying a `co-infection` role, regardless of how strong its targeted coverage and contig evidence were. Eligibility now resolves coverage across both axes; it can only ever _add_ candidates to the eligible pool, never demote one.
- **Coverage breadth below 1% scored as if it were 100%.** The percent-to-fraction coercion used a `> 1` units heuristic, so any value in `[0, 1]` was read as an already-fractional breadth. Because the first-pass breadth is rounded to an integer, every candidate at 0.5–1.5% breadth landed on exactly `1` and collected the **full** 3.0-point breadth award: at the default weights, 1% breadth scored 6.599 against 2% breadth's 3.659. Breadth is now coerced as a 0–100 percent unconditionally. (Raised as a warning during the Phase-8 review and dismissed because such candidates "always fail the gate" — no longer true once D-07/D-09 made the floor informational.)
- **`Summary.csv` reported no role statistics at all for `co-infection (indeterminate dominance)` samples.** When the indeterminate-dominance trigger fires, both candidates take `role = "indeterminate"`, which matched neither the `dominant` nor the `co-infection` filter that fills the wide slots — so `Major_role_reference`, `Minor_role_reference`, both subtypes, **both dominance scores**, both role reasons, both evidence states and the contig-support tokens in `Major_evidence`/`Minor_evidence` all read `NA`, on exactly the samples where a reader most wants the numbers. `indeterminate` candidates are now slot-eligible: the Major slot takes the lowest-ranked one and the Minor slot the next. Ordering is by **`candidate_rank`, not by dominance score**, deliberately — every other `Major_*`/`Minor_*` column in the row (reference, read counts, coverage breadth, average depth, consensus similarity, and the GLUE resistance profile) is keyed to rank, so ordering the role slots by score would have put `Major_role_reference` on one strain while the coverage and resistance columns described the other. Since an indeterminate call means precisely that the two cannot be ranked, rank order is equally defensible and keeps the row internally consistent. The uncertainty is not lost: `role` stays `indeterminate` in `candidates.csv`, both `role_reason` columns read `indeterminate_dominance_conflict`, and `overall_sample_call` is unchanged.
- **The `indeterminate` review sentence named the wrong mechanism.** It read _"No candidate passed the major-gate — overall sample call indeterminate."_, but `classify_roles()` has not applied a `minRead`/`minCov` gate to dominance since D-07/D-09 turned that floor into an informational annotation. `overall_sample_call = "indeterminate"` is reached when coverage exists but `dom_idx` is `NA` — that is, no candidate is both covered on some axis **and** `concordance_ok` — so in practice every covered candidate's mapping genotype conflicts with the genotype HCV-GLUE assigned to the same BAM, which is the only live source of `discordant`. The old wording sent the reader to `--minRead`/`--minCov` and to the first-pass statistics, neither of which has anything to do with why the sample was not called. Now reads _"No candidate eligible to be called dominant — no candidate combines measurable coverage with a mapping identity that agrees with GLUE. Overall sample call indeterminate. Please review."_ **This changes an emitted `review_flag` string**; consumers matching on the old text must be updated. The same stale wording is corrected in `docs/output_interpretation.md` (the `overall_sample_call` table and the `review_flag` trigger list).
- **The `below_floor` annotation in `candidates.csv` ignored the targeted mapping.** It compared the first-pass read count (with duplicates) and first-pass breadth against thresholds named `min_targeted_read`/`min_targeted_cov`, so it shared neither source with the score nor semantics with its own parameters, and reported every rescued candidate as failing a floor it clears by a wide margin. It now reads `targeted_reads_nodup` and `cand_cov_breadth`, falling back per row.
- **De novo "dominance" disagreement raised `review` on clean co-infections.** `summarize.R` escalated to `call_confidence = review` whenever `denovo_major_subtype_match == "NO"`, reporting that "de novo and mapping disagree on which strain is dominant". They had never been compared on the same axis: `denovo_major_ref` comes from a BLAST frame sorted by **bitscore** (`blast_parse.R:301`), which ranks strains by proximity to the reference _panel_, not by abundance. On the designed 7:3 `sim2` mixture it named 3a the de novo major on a 100% panel match despite 2.4× less k-mer coverage and 2.2× fewer mapped reads than 2a. The message is removed and both `denovo_*_subtype_match` legs are scoped to non-co-infection calls; dominance remains owned by the abundance-based D2 trigger. The `denovo_*_ref` fields deliberately keep their bitscore ordering — they also feed the off-genotype review trigger, which wants substantiality, not abundance.
- **Coverage and depth silently `NA` when a candidate had no depth file.** `df_coverage` derives both reference columns from the staged `depth/` filenames while the left side of the join takes them from `candidates.csv`. A sample with a rank-2 candidate that never reached JOINT_MAPPING carried `Minor_reference = NA` on one side and the real name on the other; dplyr matches NA-to-NA, so 1-candidate and 2-candidate/2-file samples joined while **2-candidate/1-file samples matched nothing and lost every coverage and depth column** (15 of 93 samples in the validation cohort). Because the typability gate reads that coverage, it also forced `major_typable = NO` on samples with ~692k reads at 100% breadth. `Minor_reference` is dropped from the join key.
- **`major_contig_length` could describe a different contig than `major_ref`.** The major slot derived its reference, contig and length three different ways — row 1 of `scaf_top`, the longest contig in `scaf_top` unfiltered, and the longest contig in the _full_ BLAST table filtered to `major_ref`. On `sim1` that reported `major_ref = 1a_HQ850279` (own contig 9,076 bp) beside `major_contig_length = 9,339 bp`, the length of the 1b contig. The slot now uses the same one-row discipline as the minor slot. `denovo_major_contig_length` has no logic consumers, so this corrects a reported number and cannot change a call.
- **GLUE aggregation empty after the cand-slot rename (`bin/GLUE_json_parser.R`):** The parser was globbing for `*.major.nodup.json` / `*.minor.nodup.json`, but this release renamed those slots to `cand1` / `cand2`, leaving both `GLUE_collected_report_*.tsv` files header-only and every GLUE column NA in `Summary.csv` for all samples. The glob is now a regex alternation `(cand1|major)` / `(cand2|minor)` so both new and legacy filenames are matched.
- **`daclatasvir` column-name typo in the GLUE-absent fallback (`bin/summarize.R`):** The fallback placeholder block used `daclasvir*` (missing `ta`), disagreeing with the parser's `daclatasvir*` column name. Renamed to `daclatasvir` / `daclatasvir_mut` / `daclatasvir_mut_short` so both code paths produce a consistent `Summary.csv` header.
- **De novo confirmation floor defaults settled against the SRA validation cohort:** `--denovo_min_contig_length` **500 bp**, `--denovo_min_kmer_cov` **2.0×**, `--denovo_min_blast_identity` **90%**. The previous 10.0× k-mer floor was too strict and would have refuted ERR1810453's genuine partial 2b co-infection (~5× k-mer coverage); the 500 bp length floor keeps genuine short minor contigs such as ERR1810507's (829 bp / 30× / 92%). Not to be confused with `--review_min_offgenotype_contig_length` (1000 bp), a different floor answering a different question — see `Added`.
- **`denovo_minor_contig` named a different contig than `denovo_minor_ref` and `denovo_minor_contig_length` described.** The reference was selected in `blast_parse.R` from `scaf_top` — one row per contig, keeping only contigs whose _own_ top hit is off-genotype — while the contig name was then re-derived in `summarize.R` as the top-bitscore hit to that reference across the _full_ hit table, an unrestricted search that could return a contig the selection had excluded. In one production sample: `denovo_minor_ref = 6i_DQ835770` and `denovo_minor_contig_length = 1620` (both `NODE_3`, the genuine off-genotype contig) reported alongside `denovo_minor_contig = NODE_2_length_3232` — a 1a contig whose conserved 5′UTR/core region hit the same 6i reference at 12× the bitscore. Three fields, two contigs, and an analyst sent to the wrong sequence. `blast_parse.R` now captures the minor selection as one row and emits `minor_contig` alongside `minor_ref` and `minor_contig_length`; `summarize.R` consumes it and the minor-slot re-derivation is removed. The major slot keeps its re-derivation and was never affected — `denovo_major_ref` is the globally best hit, so its row is necessarily the top-bitscore row for that reference. **`blastparse.csv` gains a `minor_contig` column**, and `denovo_minor_contig_length` now reports the selected contig rather than the longest contig sharing that reference.
- **`Major_genotype` / `Minor_genotype` were stale legacy GLUE-slot columns.** They were absent entirely from 2 of 5 runs of the same pipeline version (134 vs 136 columns), because they were only created inside the `gt_check` block, which is guarded on GLUE reports being present — so a downstream consumer reading `Major_genotype` broke on some runs and silently reported no genotype on others. They were also swapped against the subtype columns on the one sample whose roles are reversed relative to first mapping (`Major_genotype=2` alongside `Major_subtype=3a`), because `Major_subtype`/`Minor_subtype` are corrected to the role-based assignment but the genotype fields never were. Both columns are now derived from the already-corrected subtypes via `genotype_from_subtype()`, so they are 2k1b-aware, always present, and in lockstep with the subtype columns by construction. **This does change two reported values** on the affected sample; `Major`/`Minor` and `Major_subtype`/`Minor_subtype` are unchanged.
- **Variation-plot grid no longer silently empty after the cand-slot rename:** `bin/summarize.R` now splits variation-plot PNGs by `_cand1.` / `_cand2.` patterns (was `"major"` / `"minor"`, which never matched after the filename migration).
- **`plot_bam_variation.R` cand-slot extraction:** The BAM basename is now parsed by searching for a field matching `^cand[0-9]+$` (with a legacy fallback to position 3), replacing the hard-coded position-3 `str_split` that read `"nodup"` instead of the cand slot after the rename.
- **MultiQC silently merged the two FastQC read mates into one sample.** MultiQC's default `fn_clean_exts` truncates a sample name at the literal substring `.trim` wherever it appears, not only as a trailing extension. `FASTQC_TRIM` published as `${meta.id}.trimmed_1` / `_2`, so both mates truncated to bare `${meta.id}`, became the same sample, and whichever mate MultiQC parsed last silently overwrote the first — non-deterministically, since parse order varied per run. 13 files in the report inherited the flip. The prefix is now `${meta.id}_trimmed`, which carries no `.trim` substring, with the FastQC (trimmed) `path_filters` in `assets/multiqc_config.yml` following the rename.
- **Release tagging never fired the release workflow.** `auto-release.yml` triggered only on `v*.*.*`, so the bare `1.3.0` tag produced no GitHub release. The filter now matches both forms and the version-extraction step strips any leading `v`.

### `Validation`

- The dominance-score fixes listed at the top of `Fixed` were checked by re-running `summarize.R` pre- and post-fix over two existing runs (Thomson 2016 accessions, simulated co-infections, and the IVT dilution series) — **20 samples, 30 candidates, 11 multi-candidate samples**. The pre-fix rerun reproduces both runs' shipped `candidates.csv` byte-for-byte, so the comparison isolates this change. Result: **no dominance ordering flips, no `role` changes, and no `overall_sample_call` changes**; `dominance_score` moves on every candidate (−1.16 to +2.88) and `below_floor` on 5 of 30. Both anchor cases reproduce their predicted values exactly — ERR1810469's 3a/1a pair at 5.9479→5.1130 and 7.5571→8.3677, and ERR1810447's rescued 2b at 5.2307→8.1086. Re-running an existing dataset is therefore expected to change the scores and the floor annotation without changing which strains are reported.

### `Re-running pre-release results`

This applies only to the dominance-score fixes at the top of `Fixed`, and only if you already ran a pre-release `dev` build — 2.0.0 changes the pipeline end to end, so there is no in-place upgrade path from a 1.x run. Those scoring fixes are confined to the `SUMMARIZE` step: the targeted breadth the score now reads (`cov_breadth_min_5`) is computed inside `summarize.R`'s own coverage loop from depth files that already exist — so **no re-mapping, re-assembly, BLAST, GLUE or Kraken is needed** to bring an existing run up to date.

- **`nextflow run ... -resume`** re-runs exactly one process (plus MultiQC downstream of it); `SUMMARIZE` takes well under a minute on a routine run. This is the recommended route and the only one that also refreshes the MultiQC report.
- **⚠️ Editing a script in `bin/` does NOT invalidate the Nextflow cache** (verified on 24.10.1). These fixes are picked up on `-resume` only because `classify_roles.R` is passed to `SUMMARIZE` as an explicit `path()` input, and input files _are_ content-hashed. Both changed files must therefore be updated together — and a future patch touching only `bin/summarize.R` would be silently ignored by `-resume`.
- **Without the `work/` directory**, `summarize.R` can be run standalone: all eleven of its staged input directories reconstruct from the published output directory, and doing so reproduces the pipeline's `Summary.csv` and `candidates.csv` byte-for-byte. Note that the per-candidate depth files are published to **`samtools/` as `*.nodup.tsv`**, not to a `depth/` directory; that `*.nodup.idxstats` also matches the first-mapping stats, which are not staged; that `blastparse/*.csv` also matches the rescue audits, which are not staged; and that the variation directory holds two cohort-level PNGs that are not staged. Take the argument list verbatim from the run's own `.command.sh` so the run's `--minRead`/`--minCov` and score weights are preserved. This route does not regenerate the MultiQC report.
- **What moves:** `dominance_score` on every candidate, `below_floor` on some, the new `cand_cov_breadth` column, and the previously-empty role columns on indeterminate-dominance samples. `bin/tests/compare_summary_regression.R` will flag `Major_dominance_score`/`Minor_dominance_score` as disallowed differences — that is the intended change, not a regression.

### `Deprecated`

- **`Major_*` / `Minor_*` summary columns are aliased in 2.0.0** alongside the new `Major_role_*` / `Minor_role_*` equivalents. These legacy columns will be removed in the next release (COMPAT-03 drop).
- **`minor_denovo_status` and `coinfection_flag`** are no longer populated by the main reporting path; they remain in `Summary.csv` as NA-filled stubs in 2.0.0 and will be removed in the next release.

### `Removed`

- **`evidence_state == "refuted"` and the `refuted_denovo` role reason — structurally unreachable.** Contig attribution and contig contradiction were keyed on the same field: `assembly_support_join()` attaches a contig to a candidate only when their genotypes are **equal**, while the refuted predicate required them to **differ**. The band could never be taken in either `--denovo_match_level` mode. Confirmed across 358 candidate rows from 10 result directories — zero refuted, zero rows where the attributed contig's genotype differs from the candidate's — and removing the machinery leaves `Summary.csv` and `candidates.csv` **byte-identical** on a 10-sample replay. `Major_evidence_state` / `Minor_evidence_state` are now three-valued, and `any_refuted_denovo` is gone from the `call_confidence` triggers. It is not being repaired: the band required the contradicting contig to _fail_ the substantiality floors, which is the profile of assembly noise — the nine controlled IVT mixtures carry six off-genotype contigs at 226–320 bp / 0.79–1.22× k-mer that would each have refuted a correct dominant candidate. Genuine mapping-vs-assembly contradiction is owned by `RESCUE_EVALUATION`, which is sample-scoped and explicitly seeks a different-subtype contig. `concordance_status == "discordant"` and the `discordant_identity` role reason are **retained** — their GLUE leg is unaffected and live; only their de novo leg was dead.
- **`skip_assembly` removed from `nextflow_schema.json`.** The parameter was deleted from the code in `b074dbd` (assembly now always runs) but remained in the schema, where it was still offered to users, accepted on the command line, and silently did nothing.
- **`TARGETED_MAPPING` subworkflow**, superseded by `JOINT_MAPPING`.
- **Dead code removed from the repository:** `lib/NfcoreSchema.groovy`, `lib/WorkflowCommons.groovy` and `lib/WorkflowNanopore.groovy` (unreferenced; `lib/*.groovy` is auto-compiled into every run's classpath), and the unused helper scripts `bin/join_glue_report_with_summary.R` and `bin/summarize_depth.R`. Verified to leave pipeline output byte-identical.
- **`.github/workflows/release.yml`** — a manual release workflow that could not run to completion: its manifest-bump `sed` matched nothing, which failed the step before tagging. Releases are handled by `auto-release.yml` on tag push.

### `Known limitations`

- **`below_floor` is inverted in name**: the column is assigned `clears_floor`, so `TRUE` means the candidate **clears** the floor. The value is now correct; the name is not. Renaming a published column is deferred as an outward-facing contract change.
- **HCV-GLUE requires the `docker` profile.** The GLUE step is the one step not run inside its Nextflow-managed container: it drives the **host's** container runtime to start `cvrbioinformatics/gluetools-mysql` and run `cvrbioinformatics/gluetools` against each BAM, so the task needs the host Docker socket bind-mounted. It is therefore unsupported under `-profile singularity` or `-profile conda`, and unavailable on hosts that forbid non-root Docker socket access or run-time pulls from Docker Hub — many HPC sites. Run those with `--skip_hcvglue`; the rest of the pipeline is unaffected and the GLUE-derived columns are reported as `NA`, with the genotype/subtype call falling back to the mapping-based call. See the Requirements section of `README.md`.

## 1.3.0 - 2026.06.11

### `Added`

- **Subtype concordance columns in `Summary.csv`:** `denovo_major_subtype` and `denovo_minor_subtype` (extracted from the de novo/BLAST top-hit reference), plus `denovo_major_subtype_match` and `denovo_minor_subtype_match` (`YES`/`NO`/`NA`), cross-comparing the mapping-selected reference subtype against the de novo/BLAST subtype for the major and minor strain.
- **`review_flag` column in `Summary.csv`** highlighting samples that need manual inspection, written as human-readable sentences (e.g. de novo/mapping major-subtype conflict, minor candidate refuted by de novo, possible co-infection suppressed by the quality gate, or uncertain major/minor assignment). Multiple reasons are joined with `|` and the value is `NA` when nothing needs review. MultiQC now highlights any non-`NA` `review_flag` cell orange in the Results summary table via conditional formatting in `assets/multiqc_config.yml`.

### `Fixed`

- MultiQC Results summary table no longer disappears when a `review_flag` sentence contains a comma: the per-sample table handed to MultiQC is now written as TSV (`summary_mqc.tsv`) instead of CSV, sidestepping MultiQC's non-RFC-4180 comma-splitting.
- Corrected `nextflow_schema.json`, which had drifted from `nextflow.config`: the `denovo_min_contig_length` (1000 → 500) and `denovo_min_kmer_cov` (2.0 → 10.0) defaults now match the values shipped in v1.2.0, and the five `contamination_*` parameters (`hop_rate`, `min_dir_ratio`, `genome_size`, `kmer_size`, `min_aln_cov`) are now documented in the schema.
- Restored green CI on `dev`: the pipeline test now runs to completion (it was failing because the declared minimum Nextflow could not load `nf-schema@2.1.0`), and the Prettier/Black/EditorConfig linting jobs pass again. GLUE nf-test JSON fixtures (which contain GLUE `DEBUG` log lines before their JSON payload and are parsed accordingly) are excluded from Prettier via `.prettierignore`.

### `Changed`

### `Removed`

### `Dependencies`

- Raised the minimum Nextflow version to **24.04.0** (`nextflow.config` and the CI test matrix). The previous `>=23.04.0` floor could not actually run the pipeline: `nf-schema@2.1.0` requires `>=23.10.0`, and the `resourceLimits` process directive used by the test profile requires `>=24.04.0`.
- Updated deprecated GitHub Actions in the CI/linting workflows (`actions/upload-artifact` v3 → v4, `dawidd6/action-download-artifact` v2 → v6, `actions/checkout` v3 → v4, `actions/setup-node` v3 → v4, `nf-core/setup-nextflow` v1 → v2, `actions/setup-python` v4 → v5) and switched the nf-core lint step to the restructured 4.x CLI (`nf-core pipelines lint`).

### `Deprecated`

## v1.2.0 - 2026.06.09

De novo-informed strain selection: de novo/BLAST evidence and a major-gate now drive minor-strain reporting, behind `--denovo_confirm_minor` (default ON; setting it `false` reproduces pre-v1.1.7 output).

### `Added`

- **Major-gate (Change 1):** a candidate minor strain is only evaluated/reported when its major passes both `minRead` and `minCov`. A failed major reports first-mapping stats plus a `gate_flag` reason, with no minor call and no major genotype call.
- **De novo confirmation of the candidate minor (Change 2):** after first-mapping selection, the candidate minor is cross-checked against de novo/BLAST evidence and classified `confirmed_by_denovo` / `refuted` / `unconfirmed`, using a calibrated substantial-contig test (length + k-mer-coverage + BLAST identity) and genotype-level (not subtype) matching with an asymmetric refute rule. A refuted minor is downgraded to single-infection while its `Minor_*` columns are preserved for QC.
- **`minor_denovo_status` column** in `Summary.csv` and `summary_mqc.csv`, surfacing the basis for each minor call.
- Added `coinfection_flag` column to Summary.csv. When `minor_typable = NO` but `minor_denovo_status = confirmed_by_denovo` (gate suppressed the minor call while de novo assembly still confirms the minor genotype), the flag reads `possible_multiple_strains` to prompt manual review of QC plots. All other cases are `NA`.
- **New parameters:** `--denovo_confirm_minor` (default `true`), `--denovo_min_contig_length`, `--denovo_min_kmer_cov`, `--denovo_min_blast_identity`, `--denovo_match_level` (default `genotype`).
- **R regression guard** (`bin/tests/run_all.sh`) covering the major-gate, the confirm/refute/fall-back branches, genotype-level matching, flag-OFF legacy reproduction (against a committed golden baseline), and non-suppression of genuine co-infections — wired into CI as the `r-regression` job.

### `Changed`

- Lowered default `--denovo_min_contig_length` from 1000 to **500 bp** and raised default `--denovo_min_kmer_cov` from 2.0 to **10.0×**, calibrated against the SRA validation cohort: the 829 bp / ~20× ERR1810507 minor contig is now considered substantial evidence, while short spurious contigs (~300 bp / ~1–5×) remain below threshold.
- Single-sourced the 2k1b-aware genotype helper into `bin/genotype_utils.R` (`genotype_from_subtype()`), staged as a process input and used by both the selection and confirmation sides.
- BLASTPARSE per-contig CSVs (`*.blastparse.csv`, `*_blast_out.csv`) are now consumed by `SUMMARIZE` via `left_join` on `sampleName` (NA-fill on missing samples, no row loss).

### `Fixed`

- **Fixed per-sample BLAST filter no-op in `denovo_layer.R`:** `filter(sampleName == .data$sampleName)` inside `rowwise()` was comparing the column to itself (`.data` refers to the data frame, not the current row). The filter was a no-op, causing `classify_minor_denovo()` to receive pooled BLAST data from all samples. Replaced with a local variable captured before the pipe. This caused incorrect `minor_denovo_status` values — e.g. `sim11asingle` (1a single-infection with a 4g minor candidate) was falsely reported `confirmed_by_denovo` because other samples in the cohort had substantial 4g contigs.
- **Fixed secondary major-gate using first-mapping stats:** `summarize_mapping_to_all_references.R` gates on first-mapping idxstats which can be inflated by cross-mapping reads (e.g. ERR1810469: 1a reads cross-map to 3a, giving >499 first-mapping reads and >29% coverage, while targeted 3a mapping yields only 248 deduplicated reads / 24% coverage). Added a secondary gate in `summarize.R` that re-checks `Reads_nodup_mapped_major` and `Major_cov_breadth_min_5` from targeted mapping against the same `minRead`/`minCov` thresholds.
- Fixed the always-truthy `length(minor_ref > 0)` predicate in `summarize_mapping_to_all_references.R` (was `length(minor_ref) > 0`).
- Moved the minor-gate decision into the R layer (`minor_call` / `gate_flag` columns), eliminating the `NA.toInteger()` crash class in the Nextflow minor branch.
- Fixed two latent `--skip_assembly` plumbing bugs (undefined `BLASTPARSE.out`; de novo columns vanishing instead of NA-filling).
- Fixed sample mix-up risk in `TARGETED_MAPPING` subworkflow: the `reference` key is now added to the meta map before the `multiMap` split, ensuring all branches (`build`, `fasta`, `reads`) share the same meta key throughout the subworkflow.
- Fixed potential index/sample mismatch in `TARGETED_MAPPING` (bowtie2 path): `BOWTIE2_ALIGN` now receives reads, index, and fasta joined by meta key rather than positionally.
- Fixed `ggsave()` crash in `contamination_report.R` when running cohorts with more than ~53 samples. Plot cell size now scales down proportionally for large N so dimensions stay within ggplot2's 50-inch limit.
- Contamination check heatmap limited to 30 samples.

### `Removed`

- Removed the dead, non-functional `strategy == "denovo"` reference-selection branch and the undeclared `params.minDenovoLength`; the `strategy` parameter is removed from `nextflow_schema.json` and all config profiles. Reference selection now runs a single mapping-based path.
- **Breaking:** Removed the TANOTI mapper and the `--mapper` / `tanoti_stringency_1` / `tanoti_stringency_2` parameters entirely. `bowtie2` is now the only supported mapper; the mapper-selection branch and the bespoke `docker.io/jonbra/viral_haplo:1.3` image are gone. Configurations that set `--mapper tanoti` (or the stringency parameters) will no longer work. This is a non-backwards-compatible change and warrants a major-version bump.

### `Dependencies`

### `Deprecated`

## v1.1.7 - 2026.06.01

### `Added`

- Added contamination check reporting with a TSV of cross-sample contig pairs, a heatmap PNG, and a MultiQC-compatible JSON table.

### `Fixed`

- Fixed sample mix-up risk in `TARGETED_MAPPING` subworkflow: the `reference` key is now added to the meta map before the `multiMap` split, ensuring all branches (`build`, `fasta`, `reads`) share the same meta key throughout the subworkflow. Previously the enrichment happened inside the `BOWTIE2_ALIGN` input map after the split, causing `ch_aligned` to carry a different meta key than `ch_input.build` / `ch_input.fasta`, which could silently pair the wrong reference with the wrong sample in `SAMTOOLS_SORMADUP`, `STATS_WITHDUP`, `STATS_MARKDUP`, and `IVAR_CONSENSUS` during parallel multi-sample runs.
- Fixed the `reads` branch of the `multiMap` in `TARGETED_MAPPING` to emit `[meta, reads]` instead of `[meta, fasta, reads]`. The extra `fasta` element was silently bundled into the reads input of `TANOTI_ALIGN` (which expects a 2-element tuple), potentially causing alignment failures or wrong reference use in the tanoti mapper path.
- Fixed potential index/sample mismatch in `TARGETED_MAPPING` (bowtie2 path): `BOWTIE2_ALIGN` now receives reads, index, and fasta joined by meta key rather than positionally. Previously, `BOWTIE2_BUILD.out.index` was passed as a separate positional channel; since build tasks complete in non-deterministic order under parallel execution, sample A's reads could be aligned against sample B's index. The fix joins all three channels by meta key before calling `BOWTIE2_ALIGN`.

### `Dependencies`

### `Deprecated`

## v1.1.6 - 2026.02.25

### `Added`

- Do not publish fastq files from FASTP by default
- Added major coverage and major reference from the first mapping to Summary.csv

### `Fixed`

### `Dependencies`

### `Deprecated`

## v1.1.5 - 2025.11.10

### `Added`

- Added check in `bam_coverage.R` to ensure that the reference name extracted from the depth filename matches the reference name found in the depth file itself. If they do not match, the script will stop and print an error message.
- Added validation to ensure that the sample IDs from the metadata match those in the CSV file when joining channels before the MAJOR_MAPPING and MINOR_MAPPING processes. If there is a mismatch, the workflow will fail with an informative error message.
- Renamed `script_name_stringency` to `pipeline_version`

### `Fixed`

### `Dependencies`

### `Deprecated`

## v1.1.4 - 2025.10.27

### `Added`

### `Fixed`

- Updated wrong config references to KRAKEN2 database process names in the server config file.
- Limit blastparse dot plot to top 100 contigs
- Pipeline name and version are now correctly passed to the summary process and included in the final summary file.

### `Dependencies`

### `Deprecated`

## v1.1.3 - 2025.10.21

### `Added`

- Replaced custom dumpsoftwareversions module with built-in softwareVersionsToYAML functionality in the main workflow.
- Moved test datasets to a dedicated branch `test-datasets` to reduce repository size.

### `Fixed`

- Fixed wrong path name to KRAKEN2_KRAKEN2 process in the server config file.
- Correct samplesheet input now available when running the minimal test profile.

### `Dependencies`

### `Deprecated`

## v1.1.2 - 2025.10.16

### `Added`

### `Fixed`

Pipeline version is now fetched from the manifest block of `nextflow.config` and passed to the SUMMARIZE process and included in the final `Summary.csv` file.

### `Dependencies`

### `Deprecated`

## v1.1.0 and v1.1.1 - 2025.10.15

### `Added`

- **Major nf-core compliance update**: Achieved 100% nf-core lint compliance (0 failed tests) with comprehensive modernization
  - Migrated from nf-validation to nf-schema plugin (v2.1.0)
  - Updated JSON schemas to draft-2020-12 format
  - Implemented resourceLimits instead of deprecated max\_\* parameters
  - Added proper nf-test infrastructure with default.nf.test and .nftignore files
  - Created external script (`bin/run_hcvglue.sh`) to resolve Docker template string issues
  - Fixed parameter type consistency (`hcvglue_threshold` as integer)
  - Updated validation configuration for nf-schema compatibility
  - Updated all nf-core modules to latest versions
- Updated main.nf to latest nf-core template standards with proper workflow hierarchy and PIPELINE_INITIALISATION/COMPLETION subworkflows.
- Implemented nf-core standard workflow completion handlers following template patterns.
- Changed name of pipeline to `HCVTyper`.
- Changed profile name for the minimal test from `test_illumina` to `test`.
- Added plotting of variation per site in the bam files from the targeted mapping.
- The Kraken2 database PlusPFP-8 by Ben Langmead will be downloaded automatically if not overridden.
- Cleaned up various publish directories.
- Added visual representation of the denovo BLAST output.
- Removed deprecated workflows for HBV and ROV.
- No need to specify the `agens` or `platform` parameters.
- Rewrote local modules using the nf-core module create tool.
- Use Seqera containers to create conda, docker, and singularity environments.
- Changed the names of the parameters `minAgensRead` and `minAgensCov` to `minRead` and `minCov`.
- Renamed summary directory to `summary`.
- Analyzing Spades contigs instead of scaffolds. For low coverage samples, sometimes scaffolds are not produced.
- Publish blast plots in the `QC` directory.
- Allow for co-infection of 1a and 1b subtypes. Otherwise, co-infections must belong to subtypes to be identified.
- Renamed "abundance_major" and "abundance_minor" to "percent_mapped_reads_major" and "percent_mapped_reads_minor" in the summary file.
- Added total mapped reads (with duplicates) to the summary file and a calculation per sample of the fraction of mapped reads compared to the median for the entire batch.
- Samples that are filtered out during the workflow, for example due to empty fastq files, will be included in the final summary file with NA values.
- Added option to choose between cutadapt or fastp for read trimming. Default is cutadapt.

### `Fixed`

- Resolved all configuration warnings by updating process selectors to match nf-core workflow naming conventions.
- Fixed extensive linting errors in workflows/hcvtyper.nf including variable declarations and parameter naming conflicts.
- Corrected workflow.onComplete handler placement and implementation following nf-core patterns.
- The summarize R script can handle cases when GLUE report is missing. GLUE columns will all be NA.
- Fixed bug in the making of consensus sequence in cases of co-infection. The filenames would not separate between the two strains and only a single consensus would be written.
- Spades may produce empty contigs.fa file. Filter out these instances.
- Parsing the blast output and plotting may sometimes fail if there are a lot (hundreds) of contigs. Only use maximum 100 contigs for plotting (sorted by evalue and bitscore).
- Fixed a bug in the plotting of the denovo blast results where the ordering of the contigs was lost. First, the blast hits per contig are sorted first by evalue and then by bitscore and the top blast hit per contig is retained. Then, if there are more than 100 contigs, only the top 100 contigs sorted by bitscore are retained for plotting.
- Handling cases where the de novo assembled contigs produced no blast hits.

### `Dependencies`

- Updated pipeline structure to comply with latest nf-core template standards and DSL2 best practices.

### `Deprecated`

## v1.0.6 - 2025.02.12

### `Added`

Filter empty idxstats files prior to PARSEFIRSTMAPPING in the HCV workflow
Run GLUE and create json and html files for all potential major and minor strains.
GLUE is run as one single process on all bam files. To avoid conflicts with running docker images.
Compare GLUE genotypes and mapping genotypes for minor strains.
Updated tidyverse version in GLUE_PARSER.

### `Fixed`

GLUE json parser does not fail on corrupt GLUE json files.

### `Dependencies`

### `Deprecated`

## v1.0.5 - 2025.01.15

### `Added`

### `Fixed`

Joining the GLUE summary file and the sequencing summary uses tsv-files and not csv

### `Dependencies`

### `Deprecated`

## v1.0.4 - 2025.01.07

### `Added`

Ignoring errors in the SAMTOOLS_SORMADUP module after the first mapping. When using the Tanoti mapper many bam files fails in this step for some reason.

### `Fixed`

SUMMARIZE module expects tsv and not csv as output.
Reverted back to running GLUE outside of Nextflow. Some bugs in the module.

### `Dependencies`

### `Deprecated`

## v1.0.3 - 2024.12.17

### `Added`

### `Fixed`

HCV_GLUE process now runs with the docker profile.
Write final summary file as tsv and not csv

### `Dependencies`

### `Deprecated`

## v1.0.2 - 2024.12.09

### `Added`

### `Fixed`

GLUE json parser script can handle lines beginning with "DEBUG"

### `Dependencies`

### `Deprecated`

## v1.0.1 - 2024.11.24

### `Added`

Adhere versioning to Semantic Versioning.

### `Fixed`

Renamed niph to folkehelseinstituttet
Updated repo name and versions throughout

### `Dependencies`

### `Deprecated`

## v1.0 - 2024.11.08

## v1.0dev - [date]

Initial release of niph/viralseq, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
