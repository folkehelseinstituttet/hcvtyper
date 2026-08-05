# ![folkehelseinstituttet/hcvtyper](docs/images/logo-engelsk-hele-navnet.jpg#gh-light-mode-only) ![folkehelseinstituttet/hcvtyper](docs/images/logo-engelsk-hele-navnet-hvit.png#gh-dark-mode-only)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## Table of Contents

- [About HCVTyper](#about-hcvtyper)
- [How the pipeline calls strains](#how-the-pipeline-calls-strains)
- [Requirements](#requirements)
- [Run the pipeline](#run-the-pipeline)
- [Test the pipeline](#test-the-pipeline)
- [Required parameters](#required-parameters)
  - [Samplesheet input](#samplesheet-input)
  - [Output directory](#output-directory)
  - [Profiles](#profiles)
  - [Provide parameters in a file](#provide-parameters-in-a-file)
- [Optional parameters](#optional-parameters)
  - [Kraken2 databases](#kraken2-databases)
  - [HCV reference sequences](#hcv-reference-sequences)
  - [Co-infections (major and minor strains)](#co-infections-major-and-minor-strains)
  - [Candidate selection](#candidate-selection)
  - [De novo corroboration thresholds](#de-novo-corroboration-thresholds)
  - [De novo subtype rescue](#de-novo-subtype-rescue)
  - [Dominance-score weights](#dominance-score-weights)
  - [Review-flag thresholds](#review-flag-thresholds)
  - [Cross-sample contamination check](#cross-sample-contamination-check)
- [Starting and stopping the pipeline](#starting-and-stopping-the-pipeline)
- [Customizing the pipeline](#customizing-the-pipeline)
- [Output files](#output-files)
- [Citations](#citations)

## About HCVTyper

**folkehelseinstituttet/hcvtyper** is a bioinformatics pipeline used at the [Norwegian Institute of Public Health](https://www.fhi.no/en/) that is designed for highly variable viruses, and viruses that are likely to appear as co-infections between multiple strains, such as Hepatitis C Virus. The pipeline identifies the strains present in a sample sequenced on the Illumina platform, maps the reads to the corresponding references using [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml), and creates a consensus sequence per strain. For Hepatitis C Virus the pipeline can also run a [GLUE-analysis](http://hcv-glue.cvr.gla.ac.uk/#/home) to identify drug resistance mutations.

Every strain call is cross-checked against an independent line of evidence — _de novo_ assembly ([SPAdes](https://github.com/ablab/spades)) followed by BLAST against the reference panel — and the pipeline reports both **what the data indicate** and **how much to trust it**. Samples that need a human eye are labelled as such, so a marginal call never reads like a clean one.

> [!TIP]
> For how to read the results and how to confirm a call against the published evidence, see the
> [results interpretation guide](docs/output_interpretation.md).

## How the pipeline calls strains

Understanding the output is much easier with the shape of the analysis in mind:

1. **Read QC and classification.** Reads are trimmed ([fastp](https://github.com/OpenGene/fastp) or [Cutadapt](https://cutadapt.readthedocs.io/)), classified broadly with [Kraken2](https://github.com/DerrickWood/kraken2) for an overview, then classified against a small HCV-specific database. Only the HCV-classified reads go on to mapping and assembly.
2. **First mapping — neutral candidate selection.** All HCV-classified reads are mapped against the full reference panel. The top reference per distinct subtype is selected, up to `--n_candidates` (default 2). This step is deliberately **neutral**: candidates are ranked by read recruitment only, with no major/minor semantics yet.
3. **_De novo_ assembly in parallel.** SPAdes assembles the same reads, the contigs are BLASTed against the reference panel, and the results are rolled up per subtype into an assembly-support table carrying contig length, BLAST identity, alignment length and k-mer coverage.
4. **_De novo_ subtype rescue.** Where a high-quality contig contradicts a mapped candidate's subtype, the candidate's reference can be **reassigned** before the targeted mapping. Every rescue, nomination and block decision is written to a per-sample audit file.
5. **Competitive joint mapping.** Rather than mapping to each candidate independently, the reads are mapped **once** against a combined index of all candidate references, deduplicated, and then split per candidate. Reads are therefore assigned to the reference they fit best, instead of being counted several times over — which is what makes cross-mapping artefacts visible.
6. **Dominance scoring and role classification.** Each candidate gets a dominance score combining mapped reads, coverage breadth, depth evenness and k-mer coverage — with **evenness weighted above raw read count**, so a spiky high-read cross-mapping artefact loses to a genuine even minor. Each candidate is then assigned a `role` (`dominant` / `co-infection` / `background` / `indeterminate`) with a coded reason, and an evidence state (`confirmed` / `probable` / `weak`) computed from its own assembly support.
7. **Sample call and confidence.** The per-candidate roles roll up into one `overall_sample_call` and one `call_confidence` tier per sample, plus a human-readable `review_flag` naming any specific conflict.
8. **Consensus, GLUE and reporting.** Consensus sequences are generated per candidate, submitted to HCV-GLUE for genotype/subtype confirmation and drug resistance, and everything is collected into `Summary.csv` and a MultiQC report.

Candidates that are demoted are **never silently dropped** — they appear in the candidate files with the reason they were demoted.

## Requirements

The pipeline only requires [Nextflow](https://nextflow.io/) and [Docker](https://www.docker.com/) in order to run. Note that you must be able to run Docker as a non-root user as described [here](https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user).

> [!IMPORTANT]
> HCV-GLUE is currently only available with the Docker profile. We recommend that you always run the pipeline with Docker.

## Run the pipeline

The pipeline does not require any installation, only an internet connection. The pipeline is typically run with the following command:

```
nextflow run folkehelseinstituttet/hcvtyper -r v1.1.7 \
    --input samplesheet.csv \
    --outdir <OUTDIR> \
    -profile docker
```

Nextflow will pull the pipeline from the GitHub repo automatically when it is launched. Here, the version of the 1.1.7 release is downloaded and run. You can omit `-r` and the code from the master branch will be used. But we always recommend that you specify either branch or release using `-r`.

> [!NOTE]
> The evidence-based classification described in this README (neutral candidate selection, dominance
> scoring, per-candidate roles, `call_confidence`, _de novo_ rescue and competitive joint mapping) is on
> the `dev` branch and ships in the next release. To try it before then, use `-r dev`. See the
> [`CHANGELOG.md`](CHANGELOG.md) `[Unreleased]` section for the full list of changes.

If you want to download a local copy of the pipeline you can run:

```
nextflow pull folkehelseinstituttet/hcvtyper -r v1.1.7
```

Again, `-r` is optional.

## Test the pipeline

To run a minimal test:

```
nextflow run folkehelseinstituttet/hcvtyper -profile docker,test
```

This is only to see if you can get the pipeline up and running and will not run the entire pipeline such as HCV-GLUE. The results will be in a directory called `minimal_test`.

To run a full test on a real dataset type:

```
# First download the test dataset using nf-core/fetchngs
nextflow run nf-core/fetchngs -profile docker --input 'https://raw.githubusercontent.com/folkehelseinstituttet/hcvtyper/refs/heads/dev/assets/test_ids.csv' --outdir full_test

# Then run the pipeline on the downloaded dataset
nextflow run folkehelseinstituttet/hcvtyper -profile docker,test_full
```

This will download a HCV Illumina dataset from SRA and run the entire pipeline. The results will be in a directory called `full_test`.
Note that the pipeline will by default download and use the Kraken 2 [PlusPFP-8](https://benlangmead.github.io/aws-indexes/k2) database. This reqires at least 5 GB of free disk space and will take a few minutes to download and unpack. In addition, the default memory and cpu requirements of 12 cpus and 72 GB have been overridden to `50.GB` and `8`.

## Required parameters

### Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown below. The sample names can contain numbers and underscores (\_), but not spaces, dots (.) or other symbols. The fastq_1 and fastq_2 columns must contain the full path to the gzipped paired fastq files corresponding to the same sample.

```
sample,fastq_1,fastq_2
Sample_1,/path/to/sample1_fastq_R1.fastq.gz,/path/to/sample1_fastq_R2.fastq.gz
Sample_2,/path/to/sample2_fastq_R1.fastq.gz,/path/to/sample2_fastq_R2.fastq.gz
```

The samplesheet is input to the pipeline using the `--input` parameter, e.g.:
`--input assets/samplesheet_illumina.csv`

An [example samplesheet](assets/samplesheet_illumina.csv) has been provided with the pipeline in the assets directory.

**File naming requirements:**

- FASTQ files should be gzipped and paired-end
- Files should follow the naming pattern: `*_R1.fastq.gz` and `*_R2.fastq.gz` (or similar R1/R2 designation)
- All FASTQ files for a project should be organized in a single directory or subdirectories

**Creating a samplesheet automatically:**
If you have many samples, you can use the provided Docker container to automatically generate a samplesheet from a directory containing FASTQ files. The paired fastq files can be in subdirectories, and you need to point to the directory above the sub-directories. You must also point to an existing directory where you want to write the samplesheet. The container can be run like this:

```bash
# Generate samplesheet from a directory containing FASTQ files
docker run --rm \
    -v /path/to/fastq/directory:/data \
    -v /path/to/output/directory:/out \
    ghcr.io/jonbra/viralseq_utils:latest \
    /data /out/samplesheet.csv
```

### Output directory

The output directory is specified using the `--outdir` parameter, e.g.:
`--outdir results`

### Profiles

The pipeline can be run using different profiles, which will determine how the pipeline is executed. The default profile is `docker`, which uses Docker containers to run the pipeline. You can also use `singularity` or `conda` profiles if you prefer those environments. To set the profile use the `-profile` parameter, e.g.: `-profile docker/singularity/conda`.

### Provide parameters in a file

The different parameters can be provided in a file using the argument `-params-file path/to/params-file.yml`. The file can be either YAML-formatted:

```yml
input: "samplesheet.csv"
outdir: "results"
```

or JSON-formatted:

```json
{
  "input": "samplesheet.csv",
  "outdir": "results"
}
```

## Optional parameters

### Kraken2 databases

The pipeline uses [Kraken2](https://github.com/DerrickWood/kraken2) for two purposes. One is to classify the reads against a general database to get a broad overview of the taxonomic diversity within the sample (e.g., are there a lot of human reads?). The second is to classify the reads against a specific HCV-database and then use only the classified reads for the rest of the pipeline. This is done to reduce the computational load and time needed to run mapping and _de novo_ assembly.

By default, the pipeline will download and use the [PlusPFP-8 database](https://benlangmead.github.io/aws-indexes/k2) compiled by Ben Langmead for the broad classification. This requires the download and upacking of a fairly large file (>5 GB) and we recommend that you download and unpack this yourself and specify the path to the database using the `--kraken_all_db` parameter.

For the HCV-specific classification, the pipeline will use a very small and provided database which consists of around 200 different HCV strains. You can specify a custom HCV-datavase using the `--kraken_focused_db` paramter.

### HCV reference sequences

The database comes with a provided set of about 200 HCV reference sequences downloaded from NCBI. See the file [data/blast_db/HCVgenosubtypes_8.5.19_clean.fa](data/blast_db/HCVgenosubtypes_8.5.19_clean.fa). The fasta headers have been modified to begin with the genotype and subtype information (e.g., `1a`, `3b`, etc.) followed by an underscore and the NCBI accession number (e.g, `1a_AF009606`). You can for example add or remove HCV strains by modifying this file. Remember to format the fasta headers accordingly. This file will then be used in the mapping and analysis of the de novo assembled contigs to identify genotype and subtype. You need to provide the path to this file like this: `--references /path/to/HCV-sequences.fasta`.

### Co-infections (major and minor strains)

All HCV-classified reads are first mapped against the full reference panel. The best reference per distinct subtype is then selected as a **candidate** (see [Candidate selection](#candidate-selection)), and each candidate is scored and classified independently — there is no fixed "major slot" and "minor slot" during the analysis. The `Major_*` / `Minor_*` columns in `Summary.csv` are the **reported view** of that model: the dominant candidate and the highest-ranked co-infection candidate respectively.

A candidate is reported as a genuine co-infection only when it both:

1. **Clears the abundance floor** — more than `--minRead` mapped reads (default 499) and more than `--minCov` percent genome coverage at ≥5× depth (default 29). Below that it is classified `background` with reason `below_floor`.
2. **Has _de novo_ support** — a contig assembling to the same genotype (or subtype, see `--denovo_match_level`) that clears the corroboration thresholds. Without it the candidate is demoted to `background`, with reason `no_own_assembly` when no such contig exists at all, or `weak_own_assembly_below_floor` when one exists but scores below the band cut. Either way a review note names the candidate and quotes the measured contig values.

Two long-standing exceptions still apply and are unchanged:

- **Genotypes 1a and 1b** are considered distinguishable enough to be reported as a co-infection with each other. All other same-genotype pairs are collapsed — the weaker candidate is demoted with reason `same_genotype_as_dominant`.
- **The 2k/1b recombinant** present in the reference panel does not permit a co-infection with either genotype 1 or genotype 2 (reason `recombinant_2k1b`).

Demoted candidates are always written to the candidate files with their `role_reason`, so nothing disappears silently. See the [results interpretation guide](docs/output_interpretation.md) for the full role and reason glossary.

### Candidate selection

| Parameter | Default | Description |
| --- | --- | --- |
| `--n_candidates` | `2` | Number of neutrally-ranked candidate references to select and map (`cand1`..`candN`). The top reference per distinct subtype is taken. The default of 2 reproduces the classic two-slot major/minor topology. |
| `--minRead` | `499` | Minimum mapped reads for a candidate to clear the abundance floor. |
| `--minCov` | `29` | Minimum percent of the reference covered at ≥5× for a candidate to clear the abundance floor. |

### De novo corroboration thresholds

These control when an assembled contig counts as evidence for a candidate.

| Parameter | Default | Description |
| --- | --- | --- |
| `--denovo_confirm_minor` | `true` | Require orthogonal _de novo_/BLAST evidence before reporting a minor strain. |
| `--denovo_min_contig_length` | `500` | Minimum contig length (bp) to count as evidence. Validated at 500 to keep genuine short minor contigs. |
| `--denovo_min_kmer_cov` | `2.0` | Minimum SPAdes k-mer coverage of the contig. |
| `--denovo_min_blast_identity` | `90` | Minimum BLAST % identity of the contig hit. |
| `--denovo_match_level` | `genotype` | Whether _de novo_ evidence must match a candidate at `genotype` or `subtype` level. |

> [!NOTE]
> Assembly support is no longer a pass/fail gate. Each candidate receives a continuous
> `assembly_support_score` (weighted length, identity and k-mer coverage), so a contig one point below a
> threshold is no longer discarded outright. The thresholds above set the shape of that score and the
> floors used for refutation.

### De novo subtype rescue

When a high-quality contig contradicts a mapped candidate's subtype, the candidate's reference can be reassigned before targeted mapping. All four quality floors must be cleared.

| Parameter | Default | Description |
| --- | --- | --- |
| `--rescue_min_length` | `3000` | Minimum _de novo_ contig length (bp) to trigger a rescue. |
| `--rescue_min_pident` | `85` | Minimum BLAST % identity of the contig's top hit. |
| `--rescue_min_aln_length` | `3000` | Minimum BLAST alignment length (bp) of the top hit. |
| `--rescue_min_kmer_cov` | `2.0` | Minimum SPAdes k-mer coverage of the contig. |
| `--rescue_1a1b_length` | `5000` | Stricter contig-length floor when both the candidate and the _de novo_ subtype are in {1a, 1b}. |
| `--rescue_kmer_cov_ratio` | `10` | Block a replacement when the candidate's own-subtype contig k-mer coverage exceeds the rescue target's by at least this factor (cross-mapping-noise guard). `0` disables. |
| `--rescue_dominant_protect_cov` | `90` | A first-mapping candidate at or above this coverage whose own subtype _is_ assembled is self-confirming and is never replaced. Blank disables. |

A rescue that reassigns the **Major** reference always forces `call_confidence = review` and a review note — it is an automated override of the primary call. The full from/to and trigger for every rescue decision is written to `blastparse/{sample}.rescue_audit.csv`.

### Dominance-score weights

Each candidate's dominance score is a weighted sum over log10(mapped reads), coverage breadth, a CV-of-depth evenness factor, and a bonus-only log(k-mer coverage) term. Evenness deliberately outweighs raw read count so that a high-read but spiky cross-mapping artefact loses to a genuine, evenly covered minor.

| Parameter | Default | Description |
| --- | --- | --- |
| `--score_weight_evenness` | `3.0` | Weight on the CV-evenness factor — the dominant term. |
| `--score_weight_reads` | `1.0` | Weight on log10(candidate reads). |
| `--score_weight_kmercov` | `0.5` | Weight on the bonus-only log(k-mer coverage) term (never penalises). |
| `--score_evenness_k` | `1.0` | Constant in the evenness transform `1 / (1 + k × CV)`. |

### Review-flag thresholds

| Parameter | Default | Description |
| --- | --- | --- |
| `--review_min_offgenotype_contig_length` | `1000` | Minimum contig length (bp) for the monoinfection "_de novo_ assembly found a different-genotype contig" review sentence to fire. Deliberately separate from — and higher than — `--denovo_min_contig_length`: that floor confirms a minor that mapping already supports, whereas here the contig is the _only_ evidence. Validated across 140 samples: with no floor the flag fired on 51% of samples (median contig 606 bp); at 1000 bp it fires on 12% while retaining every genuine case. |

### Cross-sample contamination check

An all-vs-all BLAST of the assembled contigs across samples in the run, used to detect index hopping and cross-contamination. Results are published to `contamination_check/`.

| Parameter | Default | Description |
| --- | --- | --- |
| `--skip_contamination_check` | `false` | Skip the check entirely. |
| `--contamination_min_length` | `1000` | Minimum contig length (bp) to include. |
| `--contamination_min_id` | `95` | Minimum BLAST % identity to flag a hit. |
| `--contamination_min_aln_cov` | `0.9` | Minimum fraction of the shorter contig that must be covered by the alignment. |
| `--contamination_hop_rate` | `0.001` | Expected index-hopping rate. |
| `--contamination_min_dir_ratio` | `10.0` | Minimum source/recipient coverage ratio for a direction to be reported. |
| `--contamination_genome_size` | `9500` | Expected viral genome size (bp), used for the index-hop threshold. |
| `--contamination_kmer_size` | `127` | SPAdes final k-mer size, used to convert k-mer coverage to read depth. |

## Starting and stopping the pipeline

If the pipeline crashes, or stopped deliberately, it can be restarted from the last completed step by running the same command but with the `-resume` option. Read more about resuming a Nextflow pipeline [here](https://www.nextflow.io/docs/latest/cache-and-resume.html).

## Customizing the pipeline

Changing the arguments given to the various sub-tools can be done in several ways, perhaps the easiest is to create a custom config file. Described in more detail [here](https://nf-co.re/docs/usage/configuration#custom-configuration-files).

## Output files

The pipeline generates a comprehensive set of output files from various processes to facilitate result interpretation and quality control. By default, many intermediate files are published to help you understand the analysis. You can customize which files are published by modifying the `publishDir` settings in the configuration files. For example, to disable publishing for a specific process:

```groovy
withName: 'PROCESS_NAME' {
    publishDir = [enabled: false]
}
```

### Main output files

> [!TIP]
> This section is a reference for _what_ each file and column is. For _how to read a call_ — the
> two-axis call/confidence model, the review-flag triggers, the role glossary and a step-by-step
> procedure for confirming a call against the evidence — see the
> [results interpretation guide](docs/output_interpretation.md).

#### Summary.csv

The primary output file (`summary/Summary.csv`) containing one row per sample with genotyping and quality metrics. Key columns include:

**Read statistics:**

- `sampleName` - Sample identifier
- `total_raw_reads` - Total number of raw reads
- `total_trimmed_reads` - Reads after quality trimming
- `total_classified_reads` - Reads classified by Kraken2 as target organism
- `total_mapped_reads` - Reads mapped to all reference genomes
- `fraction_mapped_reads_vs_median` - Fraction of mapped reads relative to median across all samples. Useful for identifying outliers in a sequencing batch.

**The call and how much to trust it — read these two first:**

- `overall_sample_call` - What the data indicate: `monoinfection`, `co-infection`, `co-infection (indeterminate dominance)`, `indeterminate`, or `untypable`. Derived from the count of candidates in a `confirmed`/`probable` evidence state, independently of which candidate happens to be dominant.
- `call_confidence` - How much to trust it: `high` (no conflicting signals — accept), `provisional` (a soft caveat worth noting), `review` (a hard conflict that should block automatic reporting), or `indeterminate` (no actionable call). A sample with a non-empty `review_flag` is never `high`.

**Genotyping results:**

- `Major_genotype_mapping` / `Minor_genotype_mapping` - Identified genotypes (major/minor strains) from the reference mapping
- `Major_reference` / `Minor_reference` - Closest references identified in the mapping against all references. These were used for genotyping and re-mapping
- `Major_genotype` / `Major_subtype` / `Minor_genotype` / `Minor_subtype` - The final reported genotype and subtype per strain, derived from the role-corrected assignment. `*_genotype` is always in lockstep with `*_subtype` and is 2k/1b-aware.
- `major_typable` / `minor_typable` - Whether the strain meets quality thresholds for reliable genotyping (YES/NO)

**Per-candidate roles and evidence:**

These columns expose the role model behind the reported major/minor view.

- `Major_role_reference` / `Minor_role_reference` - The reference assigned to each role after classification (which may differ from the first-mapping reference if a _de novo_ rescue fired).
- `Major_role_subtype` / `Minor_role_subtype` - The subtype of that reference.
- `Major_dominance_score` / `Minor_dominance_score` - The combined dominance score (reads, breadth, evenness, k-mer coverage).
- `Major_role_reason` / `Minor_role_reason` - Why the candidate landed in its role: `dominant`, `corroborated`, `no_own_assembly`, `weak_own_assembly_below_floor`, `same_genotype_as_dominant`, `recombinant_2k1b`, `discordant_identity`, `indeterminate_dominance_conflict`.
- `Major_evidence_state` / `Minor_evidence_state` - The evidence band behind that reason: `confirmed`, `probable` or `weak`. Computed from the candidate's own assembly support and concordance — never forced by another candidate's dominance.
- `Major_evidence` / `Minor_evidence` - A single compact string summarising the basis for the strain, e.g. `1a_HQ850279 | 97,242 reads (nodup) | 91% breadth@10x | de novo 1a | 99.2% consensus id | ref from first-mapping`. The final token states whether the reference came from first-mapping or was reassigned by the rescue. Empty for monoinfections in the Minor slot.
- `gate_flag` - First-mapping quality gate status for the dominant strain. `ok` = passed; any other value means the major failed coverage or depth requirements and the genotype call is uncertain.

**De novo rescue:**

- `rescue_flag` - Whether any rescue survived into the final call.
- `rescue_effect` - Where a surviving rescue landed: `none`, `minor_ref_changed`, or `major_ref_changed`. `major_ref_changed` always forces `call_confidence = review` — it reassigned the dominant strain. Open `blastparse/{sample}.rescue_audit.csv` for the from/to and trigger.
- `cand_N_rescued_from` / `cand_N_rescue_trigger` - Per candidate slot, the original reference and the trigger, where a rescue applied.

**Mapping statistics (major/minor):**

- `Reads_withdup_mapped_major/minor` - Mapped reads including duplicates
- `Reads_nodup_mapped_major/minor` - Mapped reads after duplicate removal
- `Percent_reads_mapped_of_trimmed_with_dups_major/minor` - Percentage of trimmed reads that mapped, duplicates included
- `Major/Minor_cov_breadth_min_5/10` - Percentage of reference covered at ≥5× or ≥10× depth
- `Major/Minor_avg_depth` - Average sequencing depth across the reference

**HCV-specific outputs (if applicable):**

- `GLUE_genotype` / `GLUE_subtype` - Genotype and subtype determined by HCV-GLUE. "Typable" only if this matches the mapping genotype.
- `Reference` - GLUE reference sequence
- Drug resistance markers for NS3/4A inhibitors (glecaprevir, grazoprevir, paritaprevir, voxilaprevir)
- Drug resistance markers for NS5A inhibitors (daclatasvir, elbasvir, ledipasvir, ombitasvir, pibrentasvir, velpatasvir)
- Drug resistance markers for NS5B inhibitors (dasabuvir, sofosbuvir)
- `*_mut` columns - Detailed mutation information
- `*_mut_short` columns - Abbreviated mutation notation

**Technical metadata:**

- `sequencer_id` - Sequencing instrument identifier
- `pipeline_version` - Version of the folkehelseinstituttet/hcvtyper pipeline
- `HCV_project_version` - HCV-GLUE version
- `GLUE_engine_version` - GLUE engine version
- `PHE_drug_resistance_extension_version` - Version of the Public Health England (PHE) drug resistance extension applied in HCV-GLUE

**De novo confirmation columns:**

The pipeline runs _de novo_ assembly (SPAdes) and BLAST in parallel with reference mapping. These columns cross-check the two approaches and flag discrepancies for review.

- `denovo_major_ref` / `denovo_minor_ref` - Best-matching reference from the _de novo_ BLAST for the major and minor strain respectively.
- `denovo_minor_contig` - Name of the contig that produced `denovo_minor_ref`. This, `denovo_minor_ref` and `denovo_minor_contig_length` all describe **the same contig**.
- `denovo_major_contig_length` / `denovo_minor_contig_length` - Length (bp) of the supporting contig.
- `denovo_major_subtype` / `denovo_minor_subtype` - Subtype extracted from the _de novo_ BLAST hit (first field before `_` in the reference name).
- `denovo_major_subtype_match` / `denovo_minor_subtype_match` - Whether the _de novo_ subtype agrees with the reference-mapping subtype (`YES` / `NO` / `NA` if one method had no result).
- `cand_N_assembly_support` - Whether candidate slot `N` has matching assembly evidence (`none` when no contig matched).
- `cand_N_assembly_support_subtype` - Subtype of the supporting contig.
- `cand_N_assembly_support_best_contig_length` / `_pident` / `_aln_length` / `_kmer_cov` - The four measured metrics of the best supporting contig for that candidate: length (bp), BLAST % identity, BLAST alignment length (bp), and SPAdes k-mer coverage.

**Review flag:**

- `review_flag` - Human-readable summary of any issue worth manual inspection, empty when all checks pass. Multiple issues are joined with ` | `. Each sentence names **the specific candidate** it refers to and carries the **measured value** alongside the floor it missed — e.g. _"…different-genotype contig (6i) — 1620 bp contig, 69 bp aligned (4%), 91.3% identity, k-mer cov 1.0; only 4% of the contig aligns to any reference in the panel, so the subtype assignment is weakly supported"_ — so the flag can be adjudicated without re-running the pipeline. The ten triggers are listed in full in the [results interpretation guide](docs/output_interpretation.md#review_flag--the-human-readable-why).

Samples with a non-empty `review_flag` are highlighted in orange in the MultiQC Results summary table.

**Deprecated columns:**

- `minor_denovo_status` and `coinfection_flag` are no longer populated by the reporting path and remain as NA-filled stubs for one release. The information they carried is now in `Major_evidence_state` / `Minor_evidence_state`, `Minor_role_reason` and `overall_sample_call`.
- The `Major_*` / `Minor_*` columns are aliased for one release alongside their `Major_role_*` / `Minor_role_*` equivalents, and will be removed in the next release.

#### candidates.csv

`summary/candidates.csv` is the long-format companion to `Summary.csv`: **one row per candidate per sample**, including candidates that were demoted to `background`. This is where you look when you want to know why a strain was or was not reported. Columns include `candidate_rank`, `candidate_ref`, `candidate_subtype`, `candidate_reads`, `candidate_cov`, `dominance_score`, `role`, `role_reason`, `evidence_state`, all the `assembly_support_*` metrics, and:

- `evidence_summary` - A plain-language sentence describing the contig corroboration for that candidate.
- `contig_identity_contribution` / `contig_length_contribution` / `contig_kmer_contribution` - The weighted per-metric contributions to the candidate's `assembly_support_score`, so the score can be reconstructed from the file.

#### Per-sample evidence files

| File | Contents |
| --- | --- |
| `parsefirstmapping/{sample}.candidates.csv` | The neutrally-ranked candidates as selected from the first mapping, before rescue and classification. |
| `parsefirstmapping/{sample}_cand*.fa` | The selected candidate reference sequences. |
| `blastparse/{sample}.assembly_support.csv` | Per-subtype _de novo_ assembly-support roll-up: best contig length, BLAST % identity, alignment length and k-mer coverage. |
| `blastparse/{sample}.blastparse.csv` | Parsed BLAST hits of the assembled contigs against the reference panel. |
| `blastparse/{sample}_blast_out.csv` | The raw BLAST output (every HSP). |
| `blastparse/{sample}.rescued.candidates.csv` | The candidate set after _de novo_ subtype rescue. |
| `blastparse/{sample}.rescue_audit.csv` | The authoritative ledger of **every** rescue, nomination and block decision, with from/to and trigger — including rescues that fired and were then dropped by genotype collapse or the candidate cap, which leave no trace anywhere else. |

#### MultiQC Report

A comprehensive HTML report (`multiqc_report.html`) that summarizes:

- Run information and pipeline parameters
- Command line and configuration used
- Pipeline version and software versions
- Quality control metrics (FastQC, trimming statistics)
- Read classification and mapping statistics
- Genotyping results and drug resistance summaries (for HCV)
- Visualization of coverage and variant distributions

The MultiQC report provides an interactive overview of all samples and is the recommended starting point for result interpretation.

### Additional output directories

- `summary/` - `Summary.csv` and `candidates.csv` — the primary results
- `fastqc/` - Raw and trimmed read quality reports
- `fastp/` or `cutadapt/` - Read trimming logs and statistics
- `kraken2/` - Taxonomic classification reports
- `parsefirstmapping/` - Neutral candidate selection from the first mapping, and the selected reference FASTAs
- `blastparse/` - Parsed contig BLAST results, assembly support, rescued candidates and the rescue audit
- `samtools/` - BAM file statistics and mapping metrics
- `bowtie2/` - Alignment files and indices
- `ivar/` - Per-candidate consensus sequences
- `consensus/` - Consensus-vs-reference distance metrics (`Major/Minor_consensus_similarity_pct`)
- `spades/` - De novo assembly results
- `blast/` - BLAST results against reference database
- `hcvglue/` - HCV-GLUE genotyping and resistance reports (for HCV samples)
- `contamination_check/` - Cross-sample contamination and index-hopping report
- `QC/coverage_plots/`, `QC/bam_variation/`, `QC/denovo/`, `QC/insert_size/` - Per-strain coverage and evenness plots, variation grids, assembly plots and insert-size metrics
- `pipeline_info/` - Execution reports, timeline, and software versions

## Citations

If you use folkehelseinstituttet/hcvtyper for your analysis, please cite it using the following doi: [https://doi.org/10.1099/acmi.0.001193.v1](https://doi.org/10.1099/acmi.0.001193.v1)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
