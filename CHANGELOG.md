# folkehelseinstituttet/hcvtyper: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### `Added`

- **Neutral candidate selection (`params.n_candidates`, default 2):** Reference selection now ranks up to N candidates neutrally by read recruitment — no major/minor dominance semantics during the run. The top reference per distinct subtype is selected; `is_valid_minor()` validity filtering is removed from selection and moved to classification. A new `*.candidates.csv` (long-format, one row per candidate) is emitted by PARSEFIRSTMAPPING alongside the legacy wide CSV.
- **Per-genotype assembly support (`*.assembly_support.csv`):** `blast_parse.R` now emits a per-subtype assembly-support roll-up (best contig by `sc_length`) carrying four metrics: best contig length, BLAST % identity, BLAST alignment length, and k-mer coverage. These are joined to the neutral candidates at genotype level (parameterised by `--denovo_match_level`, default `genotype`) via the new `bin/assembly_support_join.R` helper, replacing the former major/minor-specific de novo evidence pass-through. Candidates with no matching assembly evidence resolve to `assembly_support = "none"` with NA metrics (no row loss).
- **Dominance scoring:** Each candidate receives a combined dominance score over log10(mapped reads), coverage breadth fraction, CV-of-depth evenness, and log10(1 + k-mer coverage). Breadth evenness is weighted 3× raw read count by default (`--score_weight_evenness 3.0`, `--score_weight_reads 1.0`, `--score_weight_kmercov 0.5`), making the score resistant to index-hopping artefacts that have high read count but uneven breadth.
- **Strain-role classification (`dominant` / `co-infection` / `background`):** The new `bin/classify_roles.R` helper classifies each candidate at the summary step using a "guilty until corroborated" rule: a non-dominant candidate is reported as co-infection only when it (a) clears the abundance floor (`minRead` + `minCov`) **and** (b) has genotype-level assembly support; otherwise it is classified `background`. Background/artefact candidates are surfaced explicitly in `candidates.csv` with a `role_reason` (e.g. `refuted_denovo`, `uncorroborated_kept`, `below_floor`) and never silently dropped.
- **`overall_sample_call` column in `Summary.csv`:** Derived from candidate roles — `monoinfection`, `co-infection`, or `indeterminate` — reported once per sample.
- **New `candidates.csv` output:** Per-sample long-format file listing every candidate (including background) with `candidate_rank`, `dominance_score`, `role`, `role_reason`, and all assembly-support metrics. Published to `summary/candidates/`.
- **New parameters:** `--score_weight_evenness`, `--score_weight_reads`, `--score_weight_kmercov`, `--score_evenness_k`.
- **`--review_min_offgenotype_contig_length` (default 1000 bp):** Contig-length floor for the monoinfection "de novo assembly found a different-genotype contig" review sentence. Deliberately separate from — and higher than — `--denovo_min_contig_length`: that floor confirms a minor that mapping already supports, whereas this trigger fires where the contig is the *only* evidence. See `Changed` below.
- **Summary.csv regression differ (`bin/tests/compare_summary_regression.R`):** Old-vs-new `Summary.csv` comparison that fails, with named examples, if any column outside an allowed set changed. Reads both files as all-character so a ragged schema or differing type inference cannot abort the comparison, and aligns on `sampleName` so row order cannot register as a difference. Not part of `run_all.sh` (it takes two file paths).
- **Threshold sweep tool (`bin/tests/offgeno_flag_sweep.R`):** Re-derives the off-genotype-contig review floor on a cohort of result directories, sweeping both candidate length definitions against the 2k1b exclusion and reporting must-keep retention and the `provisional`→`high` impact per cell. Aborts on a cohort-size mismatch or missing must-keep samples rather than emitting a plausible-looking table.
- **Regression suite extension (`bin/tests/test_compat.R`):** New auto-discovered test file covering golden strain-call reproduction (COMPAT-01), cand-slot filename lockstep (COMPAT-02), legacy + role column co-presence (COMPAT-03), and the 1a/1b co-infection and 2k/1b recombinant suppression exceptions (COMPAT-04). Runs automatically with `bash bin/tests/run_all.sh` (CI `r-regression` job, no YAML change required).

### `Fixed`

- **GLUE aggregation empty after cand-slot rename (`bin/GLUE_json_parser.R`):** The parser was globbing for `*.major.nodup.json` / `*.minor.nodup.json` but v3.0 renamed those slots to `cand1` / `cand2`, leaving both `GLUE_collected_report_*.tsv` files header-only and every GLUE column NA in `Summary.csv` for all samples. The glob on line 9 is now a regex alternation `(cand1|major)` / `(cand2|minor)` so both new and legacy filenames are matched. Driver calls on lines 259/262 updated to pass `"cand1|major"` / `"cand2|minor"`.
- **`daclatasvir` column-name typo in GLUE-absent fallback (`bin/summarize.R`):** The fallback placeholder block (lines 1220–1222) used `daclasvir*` (missing `ta`), disagreeing with the parser's `daclatasvir*` column name. Renamed to `daclatasvir` / `daclatasvir_mut` / `daclatasvir_mut_short` so both code paths produce a consistent `Summary.csv` header.
- **Reconciled de novo confirmation floor defaults** to the values validated against the SRA cohort: `--denovo_min_contig_length` 500 → **1000** bp; `--denovo_min_kmer_cov` 10.0 → **2.0**×. The prior 10.0× k-mer floor was incorrectly strict and would refute ERR1810453's genuine partial 2b co-infection (k-mer coverage ~5×).
- **`Major_genotype` / `Minor_genotype` were stale legacy GLUE-slot columns.** They were absent entirely from 2 of 5 runs of the same pipeline version (134 vs 136 columns), because they were only created inside the `gt_check` block, which is guarded on GLUE reports being present — so a downstream consumer reading `Major_genotype` broke on some runs and silently reported no genotype on others. They were also swapped against the subtype columns on the one sample whose roles are reversed relative to first-mapping (`Major_genotype=2` alongside `Major_subtype=3a`), because `Major_subtype`/`Minor_subtype` are corrected to the role-based assignment but the genotype fields never were. Both columns are now derived from the already-corrected subtypes via `genotype_from_subtype()`, so they are 2k1b-aware, always present, and in lockstep with the subtype columns by construction. **This does change two reported values** on the affected sample; `Major`/`Minor` and `Major_subtype`/`Minor_subtype` are unchanged.
- **Variation-plot grid no longer silently empty after the cand-slot rename:** `bin/summarize.R` now splits variation-plot PNGs by `_cand1.` / `_cand2.` patterns (was `"major"` / `"minor"`, which never matched after the filename migration).
- **`plot_bam_variation.R` cand-slot extraction:** The BAM basename is now parsed by searching for a field matching `^cand[0-9]+$` (with a legacy fallback to position 3), replacing the hard-coded position-3 `str_split` that read `"nodup"` instead of the cand slot after the Phase-9 rename.

### `Changed`

- **Output filename slot migrated `.major.` / `.minor.` → `.cand1.` / `.cand2.`** across all TARGETED_MAPPING outputs (BAM, stats, consensus). `bin/summarize.R` recovers `candidate_rank` via a join against `*.candidates.csv` rather than parsing a hard-coded filename field.
- **Uniform per-candidate TARGETED_MAPPING fan-out:** The asymmetric `MAJOR_MAPPING` / `MINOR_MAPPING` subworkflow aliases are replaced by a single `TARGETED_MAPPING` call fanning out over all rows of `PARSEFIRSTMAPPING.out.candidates` (`splitCsv.flatMap`). Per-candidate `confirmation_status` replaces the former `gate_flag` / `minor_call` plumbing.
- **`Major_role_*` / `Minor_role_*` columns added to `Summary.csv`** alongside the new overall call: `Major_role_reference`, `Major_role_subtype`, `Major_role_dominance_score`, `Minor_role_reference`, `Minor_role_subtype`, `Minor_role_dominance_score`, `overall_sample_call`.
- The legacy `apply_denovo_layer` / `minor_denovo_status` / `coinfection_flag` classification path in `bin/summarize.R` is retired; `bin/denovo_layer.R` is still staged and unit-tested but is no longer called from the main reporting path. `review_flag` is now set from per-sample role roll-up.
- **Continuous assembly-support scoring replaces the binary `own_substantial` threshold:** each candidate's own de novo/BLAST assembly support is now expressed as a continuous `assembly_support_score` (weighted length + identity + k-mer coverage, with a logistic identity term centered at 88%) instead of an ANDed length≥500bp / k-mer≥2 / identity≥90% pass/fail gate. A candidate one point below the old 90% identity cliff with a near-full-length, concordant contig is no longer silently discarded to `background`.
- **Per-candidate evidence state (`confirmed` / `probable` / `weak` / `refuted`)** replaces the binary per-candidate `own_substantial` flag. A candidate's state is computed entirely from its own score and concordance and is never forced by another candidate's dominance. `refuted` now requires the candidate's own assembly to exist AND point to a genuinely different genotype than its own mapping call AND fail quality thresholds — an absent assembly, or a weak-but-concordant one, now yields `weak` rather than `refuted`.
- **`overall_sample_call` is now derived from the count of `confirmed`/`probable` candidates**, independent of dominance ordering — a strong non-dominant candidate is reported as co-infection rather than gated out by which candidate happens to be dominant.
- **The "different-genotype contig" review sentence now carries the measured evidence.** It previously named a subtype and nothing else, then asked for a manual review — while the aligned length, identity and k-mer coverage of that contig existed only in `blastparse/*.assembly_support.csv`, never reaching `Summary.csv` (assembly support is joined at candidate grain, and an off-genotype contig is not a candidate). The sentence now reads e.g. *"…different-genotype contig (6i) — 1620 bp contig, 69 bp aligned (4%), 91.3% identity, k-mer cov 1.0; only 4% of the contig aligns to any reference in the panel, so the subtype assignment is weakly supported — the contig may be largely non-HCV, chimeric, or too divergent to type."*, while a contig aligning over its full length keeps the co-infection wording. The **aligned fraction**, not the contig length, is what separates the two: genuine minors align over 99–100% of their contigs. No flag is suppressed by this — the analyst decides, with the numbers present.
- **The monoinfection "different-genotype contig" review sentence is now gated on contig length and the 2k/1b pair rule.** It previously applied no substantiality floor at all: across five routine runs it fired on 72 of 140 samples (51%), was the only sentence on every one of them, and so accounted for 92% of the cohort's `provisional` calls — median triggering contig 606 bp, shortest 142 bp, i.e. below even `--denovo_min_contig_length`. Two gates now apply: the contig must reach `--review_min_offgenotype_contig_length` (1000 bp), and the pair must not be a 2k/1b recombinant against a genotype 1 or 2 major (the same `is_valid_minor()` rule 2 already applied to candidate promotion, which the review trigger never consulted). Length is the only leg used: the three legacy-typable minors this build demotes to monoinfection sit at k-mer coverage 1.42–1.97, *below* `--denovo_min_kmer_cov`, so a k-mer leg would suppress exactly the samples most worth reviewing. Policy lives in a new pure helper `offgenotype_contig_reviewable()`, replacing a hand-rolled `substr(x, 1, 1)` comparison. **Affects `review_flag` and `call_confidence` only** — no subtype, typability, `overall_sample_call` or resistance value changes.
- **`review_flag` messages now name the specific candidate** they refer to (rank/slot, reference, or subtype) and carry the concrete measured value alongside the floor or expectation it missed (e.g. "candidate 2 contig matched 2c but mapping says 3a"; "identity 89.0 below 90 floor"), rather than only naming the trigger — including flags raised by a non-dominant candidate in an otherwise-monoinfection sample.

### `Deprecated`

- **`Major_*` / `Minor_*` summary columns are aliased for one release** alongside the new `Major_role_*` / `Minor_role_*` equivalents. These legacy columns will be removed in the next release (COMPAT-03 drop).
- **`minor_denovo_status` and `coinfection_flag`** are no longer populated by the main reporting path; they remain in `Summary.csv` as NA-filled stubs for one release.

### `Removed`

### `Dependencies`

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
