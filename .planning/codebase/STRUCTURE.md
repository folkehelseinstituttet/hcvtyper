# Directory Structure

**Analysis Date:** 2026-06-05

## Overview

`hcvtyper` is an nf-core-style Nextflow DSL2 pipeline for HCV (Hepatitis C Virus) genotyping/typing. It follows the standard nf-core project layout: a thin `main.nf` entry point, a primary workflow in `workflows/`, reusable `subworkflows/` and `modules/`, helper scripts in `bin/`, and layered configuration in `conf/`.

## Top-Level Layout

```
hcvtyper/
├── main.nf                              # Pipeline entry point (DSL2)
├── nextflow.config                      # Root config: params, profiles, includes conf/*
├── nextflow_schema.json                 # Parameter schema (nf-core validation)
├── modules.json                         # nf-core module/subworkflow lockfile
├── nf-test.config                       # nf-test configuration
├── contamination_check_standalone.nf    # Standalone entry point for contamination check
├── cloudgene.yaml                       # Cloudgene app descriptor (web deployment)
├── CHANGELOG.md / CITATIONS.md / LICENSE
├── hcvtyper_handoff_denovo_informed_selection.md  # Handoff spec (untracked work — see CONCERNS.md)
├── assets/                              # Static assets (schemas, email templates, MultiQC config)
├── bin/                                 # Executable helper scripts (R, Python, shell)
├── conf/                                # Layered Nextflow config files
├── docs/                                # Documentation (usage, output, images)
├── lib/                                 # Groovy library classes (WorkflowHCVTyper, etc.)
├── modules/                            # Process definitions (local + nf-core)
├── subworkflows/                        # Composed multi-process units (local + nf-core)
├── workflows/                           # Top-level workflow (hcvtyper.nf)
├── tests/                               # nf-test definitions and test data
├── data/                                # Local reference data (blast_db, kraken_db, reads) — gitignored
├── minimal_test/                        # Example/minimal run outputs (per-tool)
├── contamination_runs/                  # Contamination check run outputs
└── .github/                             # CI workflows (ci, linting, fix-linting)
```

## Key Locations

### Entry Points
- `main.nf` — primary pipeline entry; includes and calls the `HCVTYPER` workflow
- `workflows/hcvtyper.nf` — main workflow orchestration (`workflow HCVTYPER { main: ... emit: ... }`)
- `contamination_check_standalone.nf` — standalone contamination-check entry point

### Processes — `modules/local/`
Each subdirectory (or `.nf` file) is one process module:
- `blastparse/` — parse BLAST output, select major/minor references, emit contigs/FASTAs
- `parsefirstmapping/` — parse first-pass mapping results
- `bamvariation/` — BAM-level variation analysis
- `consensus_distance/` — distance between consensus and reference (major + minor)
- `contig_filter_rename/`, `cat_filtered_contigs/` — de novo contig handling
- `glueparse/`, `hcvglue/` — HCV-GLUE genotyping integration
- `contamination_report/`, `plotcoverage/`, `instrumentid/`, `summarize/` — reporting
- `samplesheet_check.nf`, `tanoti.nf` — single-file process modules
- `modules/nf-core/` — vendored nf-core modules (managed via `modules.json`)

### Subworkflows — `subworkflows/local/`
- `input_check.nf` — samplesheet validation → read channels
- `targeted_mapping/` — reference-targeted mapping (aliased as `MAJOR_MAPPING` / `MINOR_MAPPING`)
- `get_mapping_stats/` — mapping statistics (aliased WITHDUP / MARKDUP)
- `contamination_check/` — contamination detection branch
- `utils_nfcore_hcvtyper_pipeline.nf` — nf-core pipeline utility subworkflow
- `subworkflows/nf-core/` — vendored nf-core subworkflows

### Helper Scripts — `bin/`
Called by process `script:` blocks. Mostly R, plus Python and shell:
- Python: `check_samplesheet.py`, `create_samplesheet.R`
- R (parsing/reporting): `blast_parse.R`, `GLUE_json_parser.R`, `parsefirstmapping`-related, `summarize.R`, `summarize_depth.R`, `summarize_mapping_to_all_references.R`, `consensus_distance.R`, `contamination_report.R`, `bam_coverage.R`, `plot_bam_variation.R`, `join_glue_report_with_summary.R`
- Shell: `run_hcvglue.sh`

### Configuration — `conf/`
- `base.config` — default resource allocation, error strategy (retry exit codes)
- `modules.config` / `modules_hcv.config` — per-process `ext.args`, `publishDir`, container overrides
- `test.config` / `test_full.config` — test profiles
- `server.config` — server/HPC profile
- `igenomes.config` — reference genome paths (nf-core boilerplate)

### Groovy Library — `lib/`
- `WorkflowHCVTyper.groovy`, `WorkflowCommons.groovy` — pipeline helper functions (param summaries, validation)

### Tests — `tests/`
- `default.nf.test` — full-pipeline nf-test
- `nextflow.config`, `.nftignore` — test config and snapshot ignore patterns
- `csv/samplesheet_minimal_test.csv` — minimal test input
- `blast_db/` — small test BLAST database

## Naming Conventions

| Artifact | Convention | Example |
|----------|-----------|---------|
| Process module dirs | snake_case | `modules/local/blastparse/` |
| Process names | UPPERCASE_UNDERSCORE | `BLASTPARSE`, `TANOTI_ALIGN` |
| Subworkflow files | lowercase / snake_case | `input_check.nf`, `targeted_mapping/` |
| R/Python scripts | snake_case | `summarize_mapping_to_all_references.R` |
| Groovy lib classes | PascalCase | `WorkflowHCVTyper.groovy` |
| Config files | lowercase | `base.config`, `modules.config` |

## Where to Add New Code

- **New analysis step (single tool):** add a process under `modules/local/<name>/main.nf`, wire it into `workflows/hcvtyper.nf`, add `ext.args`/`publishDir` in `conf/modules.config`.
- **New multi-process unit:** add a subworkflow under `subworkflows/local/<name>/`.
- **New helper logic:** add a script in `bin/` (match language to neighbors — R dominates), call it from the process `script:` block.
- **New parameter:** declare in `nextflow.config` params block and `nextflow_schema.json`.
- **New test expectation:** update `tests/default.nf.test` and run `nf-test test --snapshot-update`.

## Gitignored / Generated (not source)
`work/`, `results/`, `data/`, `.nextflow/`, `.nextflow.log*`, `.nf-test/`, `minimal_test/`, `contamination_runs/`, `testing_*` — runtime/output directories excluded from review and Prettier.

---

*Structure analysis: 2026-06-05*
