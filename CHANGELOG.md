# folkehelseinstituttet/viralseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### `Added`
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
