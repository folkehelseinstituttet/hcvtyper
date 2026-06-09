<!-- refreshed: 2026-06-05 -->
# Architecture

**Analysis Date:** 2026-06-05

## System Overview

```text
┌─────────────────────────────────────────────────────────────────────────────┐
│                         Main Pipeline Entry                                 │
│                  `main.nf` (DSL2 Nextflow)                                  │
│         Orchestrates PIPELINE_INITIALISATION → HCVTYPER → COMPLETION       │
└────────────────────────┬────────────────────────────────────────────────────┘
                         │
          ┌──────────────┴──────────────┐
          ▼                             ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                        HCVTYPER Workflow                                    │
│              `workflows/hcvtyper.nf` (Primary analysis engine)              │
│  Coordinates input validation, QC, assembly, mapping, genotyping pipeline   │
└───────────────────┬──────────┬──────────────┬───────────────┬──────────────┘
                    │          │              │               │
        ┌───────────┴──┐   ┌──┴──────┐   ┌──┴──────┐    ┌────┴────┐
        │              │   │         │   │         │    │         │
        ▼              ▼   ▼         ▼   ▼         ▼    ▼         ▼
  ┌──────────┐  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐
  │  Input   │  │ Read QC &    │  │ De novo      │  │ Reference    │
  │ Validation  Trimming       │  │ Assembly     │  │ Assembly &   │
  │ & Staging  & Cleanup       │  │ (optional)   │  │ Genotyping   │
  │            - Kraken2       │  │              │  │              │
  │            - FastQC        │  │ - SPAdes     │  │ - BOWTIE2/   │
  │            - FASTP/        │  │ - BLAST      │  │   TANOTI     │
  │              CUTADAPT      │  │             │  │ - GLUE        │
  └──────────┘  └──────────────┘  └──────────────┘  └──────────────┘
        │              │                  │                │
        └──────────────┴──────────────────┴────────────────┘
                       │
        ┌──────────────┴──────────────┐
        ▼                             ▼
  ┌────────────────────┐      ┌──────────────────┐
  │ Contamination      │      │ Results          │
  │ Detection          │      │ Aggregation      │
  │ (Optional)         │      │ & MultiQC        │
  │                    │      │                  │
  │ - CONTAMINATION    │      │ - SUMMARIZE      │
  │   _CHECK           │      │ - MULTIQC        │
  └────────────────────┘      └──────────────────┘
```

## Component Responsibilities

| Component | Responsibility | File |
|-----------|----------------|------|
| INPUT_CHECK | Validate and parse samplesheet CSV; stage input FASTQ files into Nextflow channels with metadata | `subworkflows/local/input_check.nf` |
| SAMPLESHEET_CHECK | Verify samplesheet format and read file existence | `modules/local/samplesheet_check.nf` |
| INSTRUMENTID | Identify sequencing instrument from FASTQ headers | `modules/local/instrumentid/main.nf` |
| FASTQC_RAW, FASTQC_TRIM | Quality assessment of raw and trimmed reads | `modules/nf-core/fastqc/main.nf` |
| FASTP / CUTADAPT | Trim adapters and filter reads (configurable) | `modules/nf-core/fastp/main.nf`, `modules/nf-core/cutadapt/main.nf` |
| PRINSEQPLUSPLUS | Optional: Remove low-complexity reads | `modules/nf-core/prinseqplusplus/main.nf` |
| KRAKEN2_KRAKEN2 | Classify reads against all domains (general contamination check) | `modules/nf-core/kraken2/kraken2/main.nf` |
| KRAKEN2_FOCUSED | Classify reads as HCV vs. non-HCV; extract HCV-classified reads | `modules/nf-core/kraken2/kraken2/main.nf` (aliased) |
| BOWTIE2_BUILD | Build BOWTIE2 index from reference FASTA | `modules/nf-core/bowtie2/build/main.nf` |
| BLAST_MAKEBLASTDB | Build BLAST nucleotide database from reference FASTA | `modules/nf-core/blast/makeblastdb/main.nf` |
| SPADES | De novo assembly of classified HCV reads | `modules/nf-core/spades/main.nf` |
| BLAST_BLASTN | BLAST assembled contigs against reference database for species/genotype inference | `modules/nf-core/blast/blastn/main.nf` |
| BLASTPARSE | Parse BLAST output; select major and minor reference sequences for detailed mapping | `modules/local/blastparse/main.nf` |
| TARGETED_MAPPING | Map reads to major/minor reference; compute consensus; generate plots (subworkflow) | `subworkflows/local/targeted_mapping/main.nf` |
| GET_MAPPING_STATS | Compute BAM indexing, idxstats, depth, samtools stats | `subworkflows/local/get_mapping_stats/main.nf` |
| PARSEFIRSTMAPPING | Identify top 2 references by read count; filter by coverage thresholds | `modules/local/parsefirstmapping/main.nf` |
| MAJOR_MAPPING | Map classified reads to highest-coverage reference; generate stats and consensus | Calls TARGETED_MAPPING |
| MINOR_MAPPING | Map classified reads to second-highest reference (if sufficient coverage); generate stats | Calls TARGETED_MAPPING |
| HCVGLUE | Run HCV-GLUE genotyping engine on consensus sequences | `modules/local/hcvglue/main.nf` |
| GLUE_PARSER | Parse HCVGLUE JSON output; extract genotype and resistance calls | `modules/local/glueparse/main.nf` |
| CONTAMINATION_CHECK | Detect cross-sample contamination via all-vs-all BLAST on assembled contigs | `subworkflows/local/contamination_check/main.nf` |
| SUMMARIZE | Aggregate results from all samples into single TSV and JSON summary files | `modules/local/summarize/main.nf` |
| MULTIQC | Generate interactive HTML QC report aggregating logs from all tools | `modules/nf-core/multiqc/main.nf` |

## Pattern Overview

**Overall:** Nextflow DSL2 modular pipeline orchestrating bioinformatics workflow for HCV genotyping and assembly from NGS data

**Key Characteristics:**
- **Multi-stage analysis:** Input validation → quality control → read classification → assembly (optional) → reference mapping → genotyping
- **Flexible mapping strategy:** Two modes via `params.strategy`: `mapping` (reference-based) or `denovo` (assembly-informed)
- **Configurable tools:** Trimmer (fastp/cutadapt), mapper (bowtie2/tanoti), optional assembly filtering and contamination detection
- **nf-core integration:** Heavy use of nf-core modules and patterns; local modules for HCV-specific logic (BLASTPARSE, GLUE, contamination detection)
- **Channel-based data flow:** Nextflow channels propagate tuples of (metadata, file paths) through subworkflows and modules

## Layers

**Input Layer:**
- Purpose: Validate samplesheet, read FASTQ files, construct metadata
- Location: `subworkflows/local/input_check.nf`, `modules/local/samplesheet_check.nf`
- Contains: CSV parsing, file existence validation, metadata enrichment
- Depends on: User samplesheet CSV
- Used by: All downstream steps

**Quality Control Layer:**
- Purpose: Assess and filter read quality; classify reads by organism
- Location: `modules/nf-core/fastqc/`, `modules/nf-core/fastp/`, `modules/nf-core/cutadapt/`, `modules/nf-core/kraken2/`
- Contains: FastQC reports, trimmed FASTQ files, Kraken2 classifications
- Depends on: Raw/trimmed reads from input layer
- Used by: Assembly layer, mapping layer

**Assembly Layer (Optional):**
- Purpose: De novo assemble classified reads; identify major/minor genotypes via BLAST
- Location: `modules/nf-core/spades/`, `modules/local/blastparse/`, `modules/nf-core/blast/blastn/`
- Contains: Contig FASTA, BLAST results, parsed genotype calls
- Depends on: Classified reads from QC layer
- Used by: Mapping strategy selection (denovo mode), contamination detection

**Reference Mapping Layer:**
- Purpose: Map classified reads to reference sequences; compute consensus; gather depth/variation statistics
- Location: `subworkflows/local/targeted_mapping/`, `modules/nf-core/bowtie2/`, `modules/nf-core/ivar/`, `modules/nf-core/samtools/`
- Contains: BAM files, depth TSV, consensus FASTA, coverage plots
- Depends on: Classified reads and reference selection from prior layers
- Used by: Genotyping layer, summary layer

**Genotyping Layer:**
- Purpose: Run HCV-GLUE on consensus sequences; extract genotype/subtype and resistance predictions
- Location: `modules/local/hcvglue/`, `modules/local/glueparse/`
- Contains: GLUE JSON/HTML reports, parsed genotype calls
- Depends on: Consensus sequences from mapping layer
- Used by: Summary layer

**Contamination Detection Layer (Optional):**
- Purpose: Identify cross-sample contig sharing via all-vs-all BLAST
- Location: `subworkflows/local/contamination_check/`, `modules/local/contig_filter_rename/`, `modules/local/contamination_report/`
- Contains: TSV pairs of potentially contaminated samples, heatmaps
- Depends on: Assembled contigs from assembly layer
- Used by: Summary/reporting

**Output Layer:**
- Purpose: Aggregate all results into single report; generate MultiQC
- Location: `modules/local/summarize/`, `modules/nf-core/multiqc/`
- Contains: Summary TSV (one row per sample), MultiQC HTML report
- Depends on: Outputs from all prior layers
- Used by: User inspection

## Data Flow

### Primary Request Path (Reference-Based Mapping)

1. **Input parsing** — `INPUT_CHECK` loads samplesheet, creates channel `[val(meta), [fastq_1, fastq_2]]` (`subworkflows/local/input_check.nf:14`)
2. **Sequencer ID detection** — `INSTRUMENTID` identifies instrument from FASTQ headers (`modules/local/instrumentid/main.nf`)
3. **Index preparation** — `BOWTIE2_BUILD` and `BLAST_MAKEBLASTDB` create alignment/search indices from reference FASTA (`workflows/hcvtyper.nf:142-153`)
4. **Raw QC** — `FASTQC_RAW` generates quality metrics on raw reads (`workflows/hcvtyper.nf:158-161`)
5. **Trimming** — `FASTP` or `CUTADAPT` removes adapters; filters if specified (`workflows/hcvtyper.nf:166-187`)
6. **Trimmed QC** — `FASTQC_TRIM` assesses trimmed read quality (`workflows/hcvtyper.nf:192-195`)
7. **Low-complexity filtering** (optional) — `PRINSEQPLUSPLUS` removes low-complexity sequences if enabled (`workflows/hcvtyper.nf:203-227`)
8. **General classification** — `KRAKEN2_KRAKEN2` classifies against all domains (contamination check) (`workflows/hcvtyper.nf:232-238`)
9. **HCV-specific classification** — `KRAKEN2_FOCUSED` extracts HCV-classified reads from input (`workflows/hcvtyper.nf:243-248`)
10. **First-pass mapping (all references)** — `BOWTIE2_ALIGN` or `TANOTI` maps classified reads to all reference sequences; `GET_MAPPING_STATS_WITHDUP` computes coverage (`workflows/hcvtyper.nf:321-350`)
11. **Reference selection** — `PARSEFIRSTMAPPING` identifies top 2 references by mapped read count and coverage (`workflows/hcvtyper.nf:384-388`)
12. **Major mapping** — `MAJOR_MAPPING` (→ `TARGETED_MAPPING`) maps reads to highest-coverage reference; generates consensus via `IVAR_CONSENSUS`; computes stats and depth (`workflows/hcvtyper.nf:421-424`)
13. **Minor mapping** — `MINOR_MAPPING` (→ `TARGETED_MAPPING`) repeats for second reference if sufficient coverage (`workflows/hcvtyper.nf:474-476`)
14. **Genotyping** — `HCVGLUE` runs against collected consensus files; `GLUE_PARSER` extracts results (`workflows/hcvtyper.nf:481-493`)
15. **Result aggregation** — `SUMMARIZE` collects all metrics into single summary file (`workflows/hcvtyper.nf:524-542`)
16. **QC report** — `MULTIQC` aggregates all logs and reports into HTML (`workflows/hcvtyper.nf:586-593`)

### De Novo Assembly Path (Alternative Strategy)

1. **Steps 1–9** same as reference-based path
2. **Assembly** — `SPADES` performs de novo assembly on classified reads (`workflows/hcvtyper.nf:266-271`)
3. **Contig search** — `BLAST_BLASTN` searches contigs against reference; `BLASTPARSE` identifies major and minor contigs by BLAST hit quality (`workflows/hcvtyper.nf:285-303`)
4. **Major mapping** — `MAJOR_MAPPING` maps classified reads to best-hit contig (from BLASTPARSE) (`workflows/hcvtyper.nf:421-424`)
5. **Minor mapping** — `MINOR_MAPPING` maps reads to secondary contig if contig length > `params.minDenovoLength` (`workflows/hcvtyper.nf:468-472`)
6. **Steps 14–16** same as reference-based

### Contamination Detection Branch (Optional)

- Runs in parallel if `!params.skip_contamination_check` and assembly is enabled
- `CONTAMINATION_CHECK` filters contigs ≥ min length, concatenates all samples' FASTAs, runs BLAST all-vs-all (`subworkflows/local/contamination_check/main.nf:1-87`)
- Produces TSV pairs and heatmap visualization (`subworkflows/local/contamination_check/main.nf:79-85`)

**State Management:**
- Metadata (`val(meta)`) carries sample ID, single-end flag, and optional fields added during processing (e.g., `reference` name for mapping)
- Channels use Nextflow tuple structures: `[val(meta), path(...)]` for file pairs
- Filter operations reduce channels on coverage/read-count thresholds (e.g., `params.minRead`, `params.minCov`)
- `.join()` operations synchronize multi-input modules by metadata key (critical in `TARGETED_MAPPING` to avoid index misalignment)

## Key Abstractions

**Metadata Map (val(meta)):**
- Purpose: Carries immutable sample context and enriched parameters through pipeline
- Examples: `[id: "sample_name", single_end: false, reference: "HCV_1a", major_reads: 1000, major_cov: 50]`
- Pattern: Tuples always pair metadata with file paths; filter/join operations use metadata keys

**Reference Selection Channel:**
- Purpose: Routes classified reads and reference FASTA pairs to mapping processes
- Examples: Outputs from `PARSEFIRSTMAPPING` (mapping mode) or `BLASTPARSE` (denovo mode)
- Pattern: Lazy evaluation; channels split into major and minor branches based on coverage criteria

**Mapped BAM Channel:**
- Purpose: Carries aligned read files through statistics and consensus generation
- Examples: Output from `BOWTIE2_ALIGN`, `TANOTI`, forwarded to samtools processes
- Pattern: Indexed with `.bai` files for fast random access

**Contamination Detection Subworkflow:**
- Purpose: Encapsulates all-vs-all BLAST logic independent of main pipeline
- Examples: Reusable via `contamination_check_standalone.nf`; can be invoked with custom contig/fastp/GLUE directories
- Pattern: Takes optional channels; emits TSV, heatmap, and MultiQC JSON

## Entry Points

**Main Pipeline (Standard):**
- Location: `main.nf`
- Triggers: `nextflow run main.nf --input samplesheet.csv --outdir results/`
- Responsibilities: Load main workflow; run pipeline initialisation (logging, version info) → run HCVTYPER workflow → run completion tasks (email, cleanup)

**Standalone Contamination Check:**
- Location: `contamination_check_standalone.nf`
- Triggers: `nextflow run contamination_check_standalone.nf --contigs_dir spades_output/ --outdir results/`
- Responsibilities: Bypass main pipeline; run only CONTAMINATION_CHECK subworkflow on user-provided contig directory (allows rerunning contamination detection with different thresholds)

**Workflow Executor (HCVTYPER):**
- Location: `workflows/hcvtyper.nf:77`
- Triggers: Implicitly called from `main.nf`
- Responsibilities: Orchestrate all analysis steps (QC → assembly → mapping → genotyping); emit MultiQC report

## Architectural Constraints

- **Threading:** Nextflow uses multi-threaded event loop; individual processes may use multiple CPUs per task label (e.g., `process_low` = 2 CPUs, `process_high` = 12 CPUs)
- **Global state:** Nextflow parameters defined in `nextflow.config` and `conf/*.config` files; no mutable global state in Groovy code. Reference files (`params.references`, `params.kraken_*_db`) are immutable channels
- **Circular imports:** None detected; module dependencies form a DAG (input validation → QC → assembly/mapping → genotyping → summarization)
- **Conditional branching:** Strategy selection (`params.strategy: mapping|denovo`), mapper selection (`params.mapper: bowtie2|tanoti`), trimmer selection (`params.trimmer: fastp|cutadapt`), and optional steps (assembly, contamination check, GLUE) all use Nextflow `if/else` blocks to enable/disable entire subworkflows
- **Channel synchronization:** Critical joins in `TARGETED_MAPPING.nf` (lines 51–53) explicitly join by metadata key to avoid mis-pairing indices with samples during parallel execution
- **Reference database handling:** Kraken2 databases can be `.tar.gz` (auto-extracted via `UNTAR`) or pre-extracted directories. Extraction is conditional and version-tracked

## Anti-Patterns

### Implicit Channel Ordering Assumption in First Implementation (Fixed in TARGETED_MAPPING)

**What happens:** Earlier code might have assumed positional pairing: `ch_reads.join(ch_index)` where indices from `BOWTIE2_BUILD` are emitted in completion order, not submission order
**Why it's wrong:** When multiple samples run in parallel, downstream join may pair wrong index with wrong sample if execution times vary
**Do this instead:** Use explicit `.join()` by metadata key as in `subworkflows/local/targeted_mapping/main.nf:51-53`, ensuring each sample's reads pair with its own index

### Dummy File Fallback Instead of Optional Outputs

**What happens:** When assembly is skipped, pipeline passes `file("dummy_file")` to MULTIQC and SUMMARIZE for BLAST and GLUE outputs (`workflows/hcvtyper.nf:514, 519`)
**Why it's wrong:** Creates misleading artifact files in output; complicates post-processing
**Do this instead:** Use optional output channels with `.ifEmpty([])` and handle in downstream modules (e.g., check file size before processing)

## Error Handling

**Strategy:** Conservative error reporting with selective retry

**Patterns:**
- `errorStrategy = "finish"` (default): Continue pipeline on process failure; fail overall job if non-retryable
- `errorStrategy = "retry"` with `maxRetries = 1` (base config): Retry once on resource/timeout errors (exit codes 130–145, 104, 175)
- `errorStrategy = "ignore"` (TANOTI only): Ignore failures; permits downstream to decide (allows some samples to proceed if one fails alignment)
- Samplesheet validation early (`INPUT_CHECK` → `SAMPLESHEET_CHECK`) catches missing files before compute
- Module `when` clauses allow runtime skipping (e.g., `!params.skip_assembly`)

## Cross-Cutting Concerns

**Logging:** Nextflow built-in logging via `log` object in Groovy; process logs captured to `pipeline_info/execution_*.txt`

**Validation:** Parameter validation via nf-schema plugin (`plugins { id 'nf-schema@2.1.0' }`); samplesheet validation via custom Python script `bin/check_samplesheet.py`

**Authentication:** None required; all reference databases are public URLs or file paths

---

*Architecture analysis: 2026-06-05*
