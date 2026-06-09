# folkehelseinstituttet/hcvtyper: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### `Added`

### `Fixed`

### `Changed`

### `Removed`

### `Dependencies`

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
  - Implemented resourceLimits instead of deprecated max_* parameters
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
