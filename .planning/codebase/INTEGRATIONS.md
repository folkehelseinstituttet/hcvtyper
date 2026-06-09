# External Integrations

**Analysis Date:** 2026-06-05

## APIs & External Services

**Reference Sequence Databases:**
- NCBI GenBank - HCV reference sequences sourced from NCBI
  - Used for: Initial HCV reference set construction
  - Access: Via Kraken2 databases (pre-built, downloaded as needed)

**GitHub:**
- hcvtyper repository hosting
  - URLs: References at `https://github.com/folkehelseinstituttet/hcvtyper`
  - Test dataset branch: `test-datasets` branch contains BLAST and Kraken2 databases
  - Functions: Pipeline distribution, test data hosting, reference database distribution
  - SDK/Client: Nextflow's built-in Git support for pipeline pulling

**AWS S3:**
- Kraken2 database hosting
  - Service: Amazon S3 (public bucket)
  - Database: `k2_pluspfp_08gb_20250402.tar.gz` (PlusPFP-8 all-species Kraken2 database)
  - URL: `https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_08gb_20250402.tar.gz`
  - Purpose: Complete taxonomic classification of reads
  - Authentication: Public access (no credentials required)

## Data Storage

**Databases:**
- Reference HCV sequences (Local FASTA)
  - Type: FASTA text files
  - Location: `data/blast_db/HCVgenosubtypes_8.5.19_clean.fa` (~2.1MB)
  - Client: BLAST (`blast/makeblastdb` module)
  - Purpose: Sequence homology searching for contig classification
  - Format: Multi-sequence FASTA with modified headers (genotype_accession)

- Kraken2 HCV-focused database (Local/Downloaded TAR)
  - Type: Kraken2 index format (compressed TAR archive)
  - Location: `data/kraken_db/db_hepacivirus/` (unpacked)
  - Client: Kraken2 v2.1.6
  - Purpose: HCV-specific read classification before assembly/mapping
  - Size: ~500MB (unpacked)

- Kraken2 All-species database (Downloaded on-demand)
  - Type: Kraken2 index (PlusPFP-8)
  - Source: AWS S3
  - Client: Kraken2 v2.1.6
  - Purpose: Broad taxonomic classification for contamination assessment
  - Size: ~8GB (requires download on first run)
  - Parameter: `params.kraken_all_db` (default via AWS)

**File Storage:**
- Local filesystem only
  - Input: FASTQ files (paired-end Illumina reads)
  - Output: Results published to `params.outdir` (default structure in `conf/modules_hcv.config`)
  - Intermediate: Temporary files managed by Nextflow (cleaned up after process completion)

**Caching:**
- Nextflow built-in caching for workflow resumption
  - Location: `.nextflow/` directory
  - Purpose: Enable `-resume` flag for restarts from last successful checkpoint
  - Configuration: Automatic (no additional setup required)

## Authentication & Identity

**Auth Provider:**
- None required
- All data sources (GitHub, AWS S3, Docker registries) are public or use standard HTTP
- Docker registry credentials optional (default registry: quay.io - public)

**Environment Variables (Optional):**
- No required authentication environment variables for basic operation
- Optional for Docker registry credentials if using private registries:
  - Docker: Managed through Docker daemon config
  - Singularity: Credential files managed by Singularity
- Tower integration optional (requires Seqera Platform token, not used in default config)

## Monitoring & Observability

**Error Tracking:**
- None (no external error tracking service integrated)
- Local error handling via Nextflow's native error strategy
- Process-level error strategy: Retry logic (default: 1 retry on specific exit codes)

**Logs:**
- Nextflow execution logs: Standard output to console and `${params.outdir}/pipeline_info/`
- Timeline report: `execution_timeline_[timestamp].html`
- Execution report: `execution_report_[timestamp].html`
- Trace log: `execution_trace_[timestamp].txt`
- DAG visualization: `pipeline_dag_[timestamp].html`
- MultiQC aggregation: `${params.outdir}/multiqc/multiqc_report.html` (aggregates all QC metrics)

**Software Versions:**
- Tracked in `versions.yml` files for each process (`modules_hcv.config` line 21-22)
- Aggregated into final pipeline report
- Located: Output directory under each process subdirectory

## CI/CD & Deployment

**Hosting:**
- GitHub repository (`folkehelseinstituttet/hcvtyper`)
- Docker Hub/Quay.io for container images
- Seqera Wave for dynamic container builds (used by nf-core modules)

**CI Pipeline:**
- GitHub Actions workflows in `.github/workflows/`:
  - `ci.yml` - Continuous integration testing
  - `linting.yml` - Code quality checks
  - `fix-linting.yml` - Auto-fix linting issues
  - `release.yml` - Release automation
  - `auto-release.yml` - Automatic release creation
  - `release-announcements.yml` - Release notifications
  - `branch.yml` - Branch management
  - `linting_comment.yml` - PR comments for lint failures
  - `clean-up.yml` - Cleanup workflows

**Deployment Target:**
- Nextflow can execute on:
  - Local workstations (Docker/Conda profile)
  - HPC clusters (SLURM, SGE, LSF via Nextflow config)
  - Cloud platforms (AWS, Google Cloud, Azure via Nextflow)
  - Seqera Platform for centralized execution and monitoring

**Dockstore Integration:**
- Dockstore configuration: `.github/.dockstore.yml`
- Provides containerized pipeline distribution via Dockstore registry

**Cloudgene Integration:**
- Cloudgene configuration: `cloudgene.yaml`
- Enables pipeline execution on Cloudgene cloud platform

## Environment Configuration

**Required env vars:**
- None (all critical parameters set via `nextflow.config` or `-params-file`)

**Optional env vars:**
- `PYTHONNOUSERSITE=1` - Prevents Python user site packages (set automatically)
- `R_PROFILE_USER` - Custom R environment (set automatically)
- `JULIA_DEPOT_PATH` - Julia package path (set for future compatibility)
- Nextflow-specific: `NF_HOME`, `NXF_VER`, `NXF_OPTS` (for Nextflow configuration)

**Secrets location:**
- No secrets management integrated
- All configuration via public repositories and environment-agnostic defaults
- Credentials (if needed for private registries) managed by:
  - Docker: ~/.docker/config.json
  - Singularity: Singularity credential files
  - Tower: TOWER_ACCESS_TOKEN environment variable (optional)

## Webhooks & Callbacks

**Incoming:**
- Pipeline input: CSV samplesheet validated via `INPUT_CHECK` subworkflow
  - Validator: `bin/check_samplesheet.py`
  - Purpose: Validate fastq file paths and sample metadata

**Outgoing:**
- Email notifications (optional):
  - `params.email` - Success email recipient
  - `params.email_on_fail` - Failure notification recipient
  - Configured via `subworkflows/local/utils_nfcore_hcvtyper_pipeline`
  
- Webhook callbacks (optional):
  - `params.hook_url` - Generic webhook endpoint for pipeline status
  - Used for integration with external notification systems

**Data Outputs:**
- Summary CSV: `${params.outdir}/summary/Summary.csv` (main results)
- MultiQC HTML: `${params.outdir}/multiqc/multiqc_report.html` (aggregated QC)
- BLAST results: `${params.outdir}/blastparse/*.csv` (contig classification)
- BAM files: `${params.outdir}/samtools/*.bam` (aligned reads)
- Kraken2 reports: `${params.outdir}/kraken2/*report.txt` (taxonomic classification)
- HCV-GLUE analysis: `${params.outdir}/hcvglue/*.{json,html,tsv}` (resistance mutations)
- Consensus sequences: Across multiple output directories
- Coverage plots: `${params.outdir}/QC/coverage_plots/*.png`
- Contamination reports: `${params.outdir}/contamination_check/*.{tsv,png,json}` (optional)

## External Tool Dependencies

**HCV-GLUE:**
- Purpose: Detect drug resistance mutations in HCV sequences
- Implementation: Docker containers (`cvrbioinformatics/gluetools-mysql`, `cvrbioinformatics/gluetools`)
- Orchestration: `bin/run_hcvglue.sh` - Manages HCV-GLUE container lifecycle
- Database: Pre-built NCBI HCV project installed in gluetools-mysql container
- Configuration: `params.hcvglue_threshold` (default: 15)
- Status: Optional, can be skipped with `--skip_hcvglue` parameter
- Limitation: Currently only available with Docker profile (not Singularity)

**TANOTI:**
- Purpose: Read mapping with support for viral haplotypes
- Implementation: Custom Docker image (`docker.io/jonbra/viral_haplo:1.3`)
- Location: `modules/local/tanoti.nf`
- Parameters: Stringency levels via `params.tanoti_stringency_1`, `params.tanoti_stringency_2`
- Alternative: Can use Bowtie2 mapping (selected via `params.mapper`)

**BLAST:**
- Purpose: Sequence homology searching for contig classification
- Implementation: NCBI BLAST v2.17.0
- Databases: Local HCV reference FASTA converted to BLAST databases at runtime
- Module: `modules/nf-core/blast/makeblastdb/` and `modules/nf-core/blast/blastn/`

**Kraken2:**
- Purpose: Read classification against HCV-specific and all-species databases
- Implementation: Kraken2 v2.1.6
- Databases: 
  - HCV-focused (local): `data/kraken_db/db_hepacivirus/`
  - All-species (AWS S3): PlusPFP-8 database auto-downloaded
- Parameters: `params.kraken_all`, `params.kraken_focused_db`, `params.kraken_all_db`

## External Analysis Platforms

**Seqera Tower (Optional):**
- Configuration: `tower.yml`
- Purpose: Centralized pipeline monitoring and execution
- Integration: Optional (not required for basic operation)

**Seqera Wave (Automatic):**
- Provides on-demand Singularity container image builds
- Used for nf-core modules and custom R-based processes
- Reduces local build time and storage requirements
- Registry: community.wave.seqera.io and community-cr-prod.seqera.io

---

*Integration audit: 2026-06-05*
