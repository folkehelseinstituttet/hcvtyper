# Testing Patterns

**Analysis Date:** 2026-06-05

## Test Framework

**Runner:**
- nf-test (Nextflow native testing framework)
- Config: `nf-test.config` at project root

**Run Commands:**
```bash
nextflow run hcvtyper -profile test,docker --outdir ./results   # Standard test run
nf-test test tests/default.nf.test --profile test,docker        # Run nf-test
nf-test test --snapshot-update                                   # Update golden snapshots
```

**Assertion Library:**
- nf-test built-in assertions: `assert workflow.success`, `assert snapshot(...).match()`
- nf-test utility plugins for file-format validation (BAM, VCF, FASTQ)

## Test File Organization

**Location:**
- Test definition: `tests/default.nf.test`
- Test data configuration: `tests/nextflow.config`
- Snapshot ignore patterns: `tests/.nftignore`
- Test CSV input: `tests/csv/samplesheet_minimal_test.csv`
- Snapshot data: `.nf-test/tests/` (generated after first run)

**Structure:**
```
tests/
├── default.nf.test              # Test definition with assertions
├── nextflow.config              # Test-specific configuration
├── .nftignore                   # Patterns for snapshot comparison
└── csv/
    └── samplesheet_minimal_test.csv  # Test input data
```

## Test Structure

**Test Definition Format** (`tests/default.nf.test`):
```groovy
nextflow_pipeline {
    name "Test HCV Typer pipeline"
    script "../main.nf"
    tag "pipeline"

    test("-profile test") {
        when {
            params {
                outdir = "$outputDir"
            }
        }
        then {
            assertAll(
                { assert workflow.success },
                { assert snapshot(...).match() }
            )
        }
    }
}
```

**Snapshot Testing:**
```groovy
{ assert snapshot(
    workflow.trace.succeeded().size(),
    removeNextflowVersion("$outputDir/pipeline_info/hcvtyper_software_mqc_versions.yml"),
    getAllFilesFromDir(params.outdir, relative: true, includeDir: true, ignore: ['pipeline_info/*.{html,json,txt}']),
    getAllFilesFromDir(params.outdir, ignoreFile: 'tests/.nftignore')
).match() }
```
- Captures: task count, version info, output file list, file contents
- Ignores: HTML files (timestamps), JSON files (variable content), test config

## Mocking

**Framework:** nf-test stub mode (`-stub` flag) — processes execute stub blocks instead of full scripts.

**What to mock:** external tool execution (BLAST, SPAdes, samtools), heavy file I/O, version extraction (kept in stubs for reproducibility).

**What NOT to mock:** input validation (`INPUT_CHECK`), parameter parsing (Nextflow framework), channel/workflow orchestration logic.

Stub example: `modules/local/blastparse/main.nf` — creates deterministic dummy CSV/PNG/FASTA outputs and emits `versions.yml`.

## Fixtures and Test Data

**Test Data Location:**
- `tests/csv/samplesheet_minimal_test.csv` — minimal samplesheet
- External datasets: `https://raw.githubusercontent.com/nf-core/test-datasets/refs/heads/hcvtyper`

**Test Configuration** (`tests/nextflow.config`):
```groovy
params {
    modules_testdata_base_path = 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/'
    pipelines_testdata_base_path = 'https://raw.githubusercontent.com/nf-core/test-datasets/refs/heads/hcvtyper'
}

aws.client.anonymous = true  // Fixes S3 access on self-hosted runners

process {
    withName: 'HCVGLUE' {
        containerOptions = '-v /var/run/docker.sock:/var/run/docker.sock --privileged'
    }
}
```

## Coverage

- Single integration test (`tests/default.nf.test`) covering the full pipeline end-to-end
- **No unit test framework** for Python/R scripts (validated only via Black/isort formatting, not pytest)
- No explicit coverage target; snapshot-based testing ensures output stability

## Test Types

**Integration / E2E:**
- Full pipeline execution with representative data; all modules and subworkflows together
- nf-test with snapshot matching; runs in `.github/workflows/ci.yml` on push/PR against dev/master

**Unit:**
- Not formally implemented
- Python: Black formatting checks via `nf-core lint` and GitHub Actions
- R: no linting or testing (utility scripts called by Nextflow) — **coverage gap**

## CI/CD Testing

**GitHub Actions Workflows:**
- `ci.yml` — runs `nextflow run ${GITHUB_WORKSPACE} -profile test,docker --outdir ./results`; matrix across Nextflow versions
- `linting.yml` — EditorConfig checker, Prettier, Black, nf-core lint
- `fix-linting.yml` — automatic linting fixes/commits on failing lint checks

## Notable Gaps

- No unit tests for R scripts (high risk of silent failures)
- De novo confirmation / informed-selection logic untested (proposed in handoff)
- No cross-platform testing (Singularity, Conda) — Docker only in CI

---

*Testing analysis: 2026-06-05*
