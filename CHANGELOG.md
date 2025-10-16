# folkehelseinstituttet/hcvtyper: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### `Added`

### `Fixed`

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
