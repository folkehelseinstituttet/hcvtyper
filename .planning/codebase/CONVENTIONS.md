# Coding Conventions

**Analysis Date:** 2026-06-05

## Naming Patterns

**Files:**
- Nextflow processes: UPPERCASE with underscores (`SAMPLESHEET_CHECK`, `TANOTI_ALIGN`, `BLASTPARSE`)
- Python scripts: snake_case (`check_samplesheet.py`)
- R scripts: snake_case / descriptive suffix (`blast_parse.R`, `GLUE_json_parser.R`, `summarize.R`, `consensus_distance.R`)
- Groovy utility files: PascalCase (`WorkflowHCVTyper.groovy`, `WorkflowCommons.groovy`)
- Nextflow workflow files: lowercase (`hcvtyper.nf`, `input_check.nf`)

**Functions (Nextflow/Groovy):**
- camelCase: `create_fastq_channel()`, `paramsSummaryMultiqc()`, `validate_and_transform()`

**Variables (Nextflow):**
- Process input/output: lowercase with underscores (`reads`, `blast_out`, `versions`, `meta`, `references`)
- Meta maps: `meta.id`, `meta.single_end` (properties lowercase)
- Task properties: `task.ext.args`, `task.ext.prefix`, `task.cpus`
- Channel variables: prefixed with `ch_` (`ch_versions`, `ch_multiqc_config`, `ch_reads`)
- Local variables in processes: snake_case (`def args`, `def prefix`, `def samtools_command`)

**Variables (Python):**
- Class attributes: snake_case with leading underscore for private (`_sample_col`, `_first_col`, `_seen`)
- Constants: UPPERCASE (`VALID_FORMATS`)
- Module-level logger: `logger = logging.getLogger()`

**Variables (R):**
- Dataframes: snake_case (`trimmed_df`, `kraken_df`, `parsefirstmapping_df`)
- Scalars/vectors: snake_case (`sampleName`, `total_raw_reads`, `pipeline_version`)

**Process Labels (nf-core standard):**
- Resource-based: `process_single`, `process_low`, `process_medium`, `process_high`, `process_long`, `process_high_memory`
- Error handling: `error_ignore`, `error_retry`

**Process Tags:**
- Sample-based processes: `tag "$meta.id"` (per-sample tracking)
- File-based processes: `tag "$samplesheet"` (input file tracking)

## Code Style

**Formatting:**
- Line length: 120 characters (enforced by Prettier for Nextflow, Black for Python)
- Indentation: 4 spaces (Python, R), groovy standard (Nextflow)

**Linting:**
- Python: Black formatter with `line-length = 120`, target versions Python 3.7-3.10
- Nextflow/Groovy: Prettier with `printWidth: 120`
- Nextflow validation: nf-core lint framework
- EditorConfig: checked via editorconfig-checker

**Config Files:**
- `.prettierrc.yml`: `printWidth: 120` (applies to `.nf`, `.groovy`, `.yml`, `.yaml`)
- `pyproject.toml`: Black `line-length = 120`, isort `profile = "black"`, `multi_line_output = 3`
- `.prettierignore`: lists excluded directories (`bin/`, `data/`, `results/`, `work/`, `testing_*`, `contamination_runs/`)

## Import Organization

**Nextflow:**
Standard order with comment-separated sections:
1. Library imports: `import nextflow.Nextflow`, `import groovy.text.SimpleTemplateEngine`
2. Subworkflow includes: `include { WORKFLOW_NAME } from './path/to/workflow'`
3. nf-core module includes: `include { MODULE_NAME } from '../modules/nf-core/path'`
4. Local module includes: `include { MODULE_NAME } from '../modules/local/path'`

Use alias syntax for multiple instances of the same module:
```groovy
include { GET_MAPPING_STATS as GET_MAPPING_STATS_WITHDUP } from '../subworkflows/local/get_mapping_stats'
include { GET_MAPPING_STATS as GET_MAPPING_STATS_MARKDUP } from '../subworkflows/local/get_mapping_stats'
include { TARGETED_MAPPING as MAJOR_MAPPING } from '../subworkflows/local/targeted_mapping'
include { TARGETED_MAPPING as MINOR_MAPPING } from '../subworkflows/local/targeted_mapping'
```

**Python:**
Import order (enforced by isort with Black profile):
1. Standard library: `argparse`, `csv`, `logging`, `sys`
2. Collections: `from collections import Counter`
3. Path utilities: `from pathlib import Path`

**R:**
Load libraries at script start with `library()` or `suppressPackageStartupMessages()`.

## Error Handling

**Python:**
- Raise `AssertionError` with descriptive messages for validation failures
- Catch assertion errors and log with `logger.critical()` before `sys.exit(1)`
- Use specific exit codes: 1 for validation failure, 2 for file not found

Example from `bin/check_samplesheet.py`:
```python
try:
    checker.validate_and_transform(row)
except AssertionError as error:
    logger.critical(f"{str(error)} On line {i + 2}.")
    sys.exit(1)
```

**Nextflow Processes:**
- Process-level error strategy: `errorStrategy "ignore"` for optional tools (e.g., `TANOTI_ALIGN`)
- Label-based retry: `label 'error_retry'` for transient failures
- Default behavior in `conf/base.config`: retries on exit codes 130-145, 104, 175; `maxRetries = 1`
- Workflow functions: use `exit 1, "ERROR: message"` for fatal validation (e.g., `subworkflows/local/input_check.nf`)

**R:**
- Use `stop()` with message to halt execution
- Use `warning()` for non-fatal issues (e.g., `bin/summarize.R`: missing/empty log file)
- Use `try(expr, silent=TRUE)` for graceful handling of expected failures

## Logging

- Python: standard `logging` module with module-level logger; level set via `--log-level` (default WARNING)
- R: base R `print()`, `cat()`, `message()`, `warning()`
- Nextflow: `log` object; logs appear in `.nextflow.log`

## Comments

- Nextflow: `//` single-line; block comments with visual separators for sections
- Python: `#` lines; docstrings (`"""..."""`) for modules/classes/functions with Args sections
- R: `#` with optional dashes (`--------`) for visual section separation

## Function & Module Design

**Nextflow Process Structure** — all processes include:
1. Process name (UPPERCASE_UNDERSCORE)
2. `tag` directive for traceability
3. `label` directive for resource allocation
4. `conda` directive (`conda "${moduleDir}/environment.yml"` or `conda ""`)
5. `container` directive (Docker/Singularity, conditional)
6. `input:` block with explicit type declarations
7. `output:` block with named `emit:` labels
8. `when:` directive: `when: task.ext.when == null || task.ext.when`
9. `script:` or `stub:` block

**Stub Implementations:**
Provide dry-run blocks (invoked with `-stub`) that create deterministic dummy outputs matching declared outputs and include version-capture commands. Example: `modules/local/blastparse/main.nf`.

**Process Parameters (`task.ext`):**
- `task.ext.args`: additional tool-specific arguments
- `task.ext.prefix`: output filename prefix (defaults to `${meta.id}`)
- Access pattern: `def args = task.ext.args ?: ''` (null-coalescing)

**Workflow Structure:**
- `main:` block for orchestration, `emit:` block for results
- Initialize channels early: `ch_versions = Channel.empty()`
- Aggregate versions: `ch_versions = ch_versions.mix(PROCESS.out.versions)`

**Python entry point:**
```python
if __name__ == "__main__":
    sys.exit(main())
```

---

*Convention analysis: 2026-06-05*
