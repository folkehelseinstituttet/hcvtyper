# Technology Stack

**Analysis Date:** 2026-06-05

## Languages

**Primary:**
- Nextflow (DSL2) - Workflow orchestration and pipeline definition (`main.nf`, `nextflow.config`, `workflows/hcvtyper.nf`)
- R - Statistical analysis and visualization (`bin/blast_parse.R`, `bin/contamination_report.R`, `bin/summarize.R`, etc.)
- Python 3 - Utility scripts and validation (`bin/check_samplesheet.py`)
- Bash/Shell - Container orchestration and environment management (`bin/run_hcvglue.sh`)

**Secondary:**
- Groovy - Nextflow configuration syntax (embedded in `.nf` and `.config` files)

## Runtime

**Environment:**
- Nextflow 23.04.0+ (`nextflowVersion = '!>=23.04.0'` in `nextflow.config`)

**Execution Engines:**
- Docker (primary - see `/home/jon.brate/hcvtyper/nextflow.config` lines 121-133)
- Singularity (line 137-146)
- Podman (line 147-155)
- Apptainer (line 174-183)
- Conda/Mamba (line 101-119)
- Local (fallback mode)

**Package Managers:**
- Conda/Mamba - Primary dependency management across all environments
- Bioconda - Specialized channel for bioinformatics tools
- Conda-forge - Standard conda packages

## Frameworks

**Core:**
- Nextflow DSL2 - Modular workflow framework with plugins
- nf-schema plugin v2.1.0 - Parameter validation and sample sheet handling

**Data Processing:**
- Tidyverse (R package suite) - Data manipulation and visualization
- Seqinr (R) - Biological sequence analysis
- Jsonlite (R) - JSON parsing

**Quality Control:**
- FastQC v0.12.1 - Read quality assessment
- MultiQC v1.31 - Aggregate quality reports
- Picard v3.4.0 - BAM file metrics and validation

**Build/Dev:**
- Black - Python code formatting (configured in `pyproject.toml`)
- isort - Python import sorting

## Key Dependencies

**Bioinformatics - Sequence Processing:**
- fastp v1.0.1 - Read trimming and quality filtering
- cutadapt v5.0 - Adapter trimming (alternative to fastp)
- PRINSEQ++ v1.2.3 - Quality filtering and complexity analysis
- seqkit - Sequence manipulation (installed from conda-forge)

**Bioinformatics - Assembly:**
- SPAdes v4.1.0 - de novo assembly with RNA-viral mode

**Bioinformatics - Alignment/Mapping:**
- Bowtie2 v2.5.4 - Read alignment to references (`modules/nf-core/bowtie2/`)
- TANOTI - Read mapping tool (custom docker image: `docker.io/jonbra/viral_haplo:1.3`)
- Samtools v1.22.1 - BAM/SAM file processing

**Bioinformatics - Taxonomy/Classification:**
- Kraken2 v2.1.6 - Taxonomic read classification against NCBI and HCV-specific databases

**Bioinformatics - Sequence Alignment:**
- BLAST v2.17.0 - Sequence homology search for contig classification

**Bioinformatics - Consensus/Variant:**
- iVar v1.4.4 - Consensus sequence generation from BAM files

**Bioinformatics - Utilities:**
- csvtk v0.31.0 - CSV/TSV manipulation and sorting
- pigz v2.8 - Parallel gzip compression
- coreutils v9.5 - Standard Unix utilities

**Bioinformatics - Drug Resistance Analysis:**
- HCV-GLUE - Docker-based resistance mutation detection (container: `cvrbioinformatics/gluetools-mysql:latest`, `cvrbioinformatics/gluetools:latest`)

**Container/Virtualization:**
- Docker - Primary container runtime
- Singularity/Apptainer - HPC-friendly container runtime
- Podman - Open-source container management

**Wave (Seqera):**
- Community Wave registry for on-demand Singularity/Docker image builds
- Referenced in container specifications for R-based modules

## Configuration

**Environment:**
- Environment variables set in `nextflow.config` lines 210-215:
  - `PYTHONNOUSERSITE=1` - Isolate Python from user site packages
  - `R_PROFILE_USER` and `R_ENVIRON_USER` - Custom R environment
  - `JULIA_DEPOT_PATH` - Julia package path (for future use)

**Build:**
- Main config: `nextflow.config` (lines 10-74) - Default parameters and profiles
- Base config: `conf/base.config` - Resource limits (CPU, memory, time allocations)
- Module config: `conf/modules.config` - Per-process settings
- HCV-specific config: `conf/modules_hcv.config` - HCV pipeline customizations
- Test configs: `conf/test.config`, `conf/test_full.config` - Minimal and full test profiles
- Server config: `conf/server.config` - Server-specific settings

**Process Resource Defaults:**
- Standard: 1 CPU, 6GB memory, 4 hours
- Low: 2 CPU, 12GB memory, 4 hours
- Medium: 6 CPU, 36GB memory, 8 hours
- High: 12 CPU, 72GB memory, 16 hours
- Global limits: 16 CPUs, 250GB memory, 240 hours

**Execution Profiles:**
- `docker` - Docker container execution (default)
- `singularity` - Singularity/Apptainer container execution
- `conda` / `mamba` - Conda environment management
- `debug` - Debugging mode with hash dumps
- `test` - Minimal test dataset
- `test_full` - Full test dataset from SRA
- `gitpod` - Cloud IDE with resource constraints
- `arm` - ARM architecture support with AMD emulation

## Platform Requirements

**Development:**
- Nextflow 23.04.0 or later
- Docker and/or Singularity (Docker strongly recommended for HCV-GLUE)
- At least 6GB RAM per process by default
- Disk space: 5GB+ for Kraken2 databases, additional space for outputs

**Production:**
- Deployment: Linux-based HPC clusters, cloud platforms (AWS, GCP, Azure via Nextflow Tower)
- Seqera Platform (formerly Nextflow Tower) optional for centralized pipeline management
- Git repository access for pipeline installation

## Database and Reference Configuration

**Default References:**
- HCV reference sequences: GitHub-hosted FASTA file (~2.1MB)
  - URL: `https://raw.githubusercontent.com/folkehelseinstituttet/hcvtyper/test-datasets/blast_db/HCVgenosubtypes_8.5.19_clean.fa`
  - Local copy: `data/blast_db/HCVgenosubtypes_8.5.19_clean.fa`
  - Contains ~200 HCV reference genotypes/subtypes (headers formatted as genotype_accession)

**Kraken2 Databases:**
- HCV-focused database: `db_hepacivirus.tar.gz` (~200 HCV reference sequences)
- All-species database: PlusPFP-8 from AWS S3 (~8GB)
  - URL: `https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_08gb_20250402.tar.gz`
  - Auto-downloads and unpacks on first use

**Database Parameters:**
- `params.references` - Custom HCV reference FASTA file path
- `params.kraken_focused_db` - HCV-specific Kraken2 database
- `params.kraken_all_db` - Full-species Kraken2 database
- `params.kraken_all` - Boolean to run full classification (default: true)

## Container Registry

**Default Registry:** quay.io (configurable in `nextflow.config` lines 196-199)
- Docker: `quay.io`
- Singularity: `quay.io`
- Podman: `quay.io`
- Apptainer: `quay.io`

**Container Images Used:**
- nf-core modules: Biocontainers and community Wave registry
- Custom modules: Seqera Wave-built images and community.wave.seqera.io registry
- HCV-GLUE: CVR Bioinformatics registry (`cvrbioinformatics/gluetools-mysql`, `cvrbioinformatics/gluetools`)
- TANOTI: Custom community image (`docker.io/jonbra/viral_haplo:1.3`)

---

*Stack analysis: 2026-06-05*
