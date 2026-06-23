# HCVtyper Pipeline — Flow Diagram Specification

Use this file to generate a flow diagram. The pipeline has 6 main sections arranged vertically,
with one section splitting into two parallel paths that merge again.

---

## Diagram layout

```
[1. Input & Read Processing]
          |
    ______|______
   |             |
[2a. First      [2b. De novo
 Mapping]        Assembly & BLAST]
   |             |
   |_____________|
          |
[3. Reference Evaluation & Candidate Selection]
          |
[4. Joint Competitive Mapping]
          |
[5. Resistance Identification (HCV-GLUE)]
          |
[6. Final Summary & Report]
```

---

## Box descriptions

### Box 1 — Input & Read Processing

**Label:** Input & Read Processing

**Contents / processes:**
- Raw paired-end FASTQ input (Illumina short reads)
- Quality control on raw reads (FastQC)
- Adapter trimming and quality filtering (Fastp / PRINSEQ++)
- Quality control on trimmed reads (FastQC)
- Taxonomic classification against full database (Kraken2)
- Focused re-classification to extract HCV-classified reads (Kraken2 focused)
- Contamination check across samples (all-vs-all BLAST)
- Instrument ID detection

**Output to next step:** HCV-classified trimmed reads (FASTQ)

---

### Box 2a — First Read Mapping (left parallel path)

**Label:** First Read Mapping

**Contents / processes:**
- Build Bowtie2 index from comprehensive HCV reference database
- Map HCV-classified reads to all references (Bowtie2 align)
- Mark and remove duplicates (Samtools / Picard SortMarkDup)
- Parse mapping results: count reads per reference, compute coverage (PARSEFIRSTMAPPING)
- Rank candidate references by neutral read recruitment — top reference per distinct subtype, top N candidates by read count
- Output: candidates CSV with per-reference read counts, coverage, and confirmation status

**Output to merge step:** Ranked candidate list with read-based evidence

---

### Box 2b — De Novo Assembly & BLAST (right parallel path)

**Label:** De Novo Assembly & BLAST

**Contents / processes:**
- De novo genome assembly from HCV-classified reads (SPAdes)
- Build BLAST database from HCV reference sequences
- BLAST assembled contigs against HCV reference database
- Parse BLAST results: extract best contig per subtype/genotype (length, % identity, alignment length, k-mer coverage)
- Compute assembly support score per genotype

**Output to merge step:** Assembly support evidence per genotype (contig quality metrics)

---

### Box 3 — Reference Evaluation & Candidate Selection

**Label:** Reference Evaluation & Candidate Selection

**Contents / processes:**
- Join first-mapping candidates with de novo/BLAST assembly evidence at genotype level
- Concordance check: compare mapping-derived subtype with de novo-derived subtype per candidate
- Rescue evaluation: if a candidate's reference does not match de novo evidence, replace with the de novo-supported reference (RESCUE_EVALUATION)
- Dominance scoring per candidate:
  - Score = (weight_reads × log10(targeted_reads_nodup)) + (weight_evenness × breadth_fraction) + (weight_evenness × cv_evenness) + (weight_kmercov × log10(1 + k-mer_coverage))
  - Default weights: reads=1.0, evenness=3.0, k-mer=0.5
- Role classification per candidate:
  - Dominant (major strain): highest dominance score, passes read/coverage floor
  - Co-infection (minor strain): clears abundance floor AND has genotype-level assembly support
  - Background/artefact: below floor or no assembly support (reported but not called)
- Overall sample call: monoinfection / co-infection / indeterminate
- Output: finalized candidate list with references, roles, dominance scores, rescue flags

**Output to next step:** Final candidate FASTA references + candidates CSV

---

### Box 4 — Joint Competitive Mapping

**Label:** Joint Competitive Mapping

**Contents / processes:**
- Concatenate all candidate reference FASTAs into one combined multi-reference FASTA
- Build single Bowtie2 index from combined reference
- Map HCV-classified reads competitively: each read assigned to its single best-matching candidate only (no read counted in multiple candidates)
- Mark and remove duplicates on combined BAM (Samtools SortMarkDup)
- Collect read-count statistics: with-duplicate and deduplicated counts per candidate reference (Samtools idxstats)
- Split combined BAM into per-candidate BAMs by reference region
- Reheader per-candidate BAMs to single-reference headers
- Per-candidate downstream processing:
  - Coverage depth profiling (Samtools depth)
  - Coverage breadth and evenness statistics
  - Coverage plot per candidate
  - BAM variation plot
  - Consensus sequence generation (iVar consensus)
  - Insert size metrics (Picard CollectInsertSizeMetrics)
  - Consensus distance calculation between candidates

**Output to next step:** Per-candidate consensus sequences + BAMs + coverage statistics

---

### Box 5 — Resistance Identification (HCV-GLUE)

**Label:** Resistance Identification (HCV-GLUE)

**Contents / processes:**
- Run HCV-GLUE per candidate BAM: genotype confirmation and drug resistance mutation analysis
- Parse GLUE JSON output (HCV_GLUE_PARSER)
- Collect and aggregate per-sample GLUE reports (major and minor candidate)
- Output: per-candidate subtype confirmation, resistance mutation calls, GLUE HTML report

**Output to next step:** GLUE TSV reports (subtype, resistance mutations) per candidate

---

### Box 6 — Final Summary & Report

**Label:** Final Summary & Report

**Contents / processes:**
- Aggregate all per-sample results: read counts, mapping stats, assembly support, dominance scores, role classifications, rescue flags, GLUE resistance calls (SUMMARIZE / bin/summarize.R)
- Compute final per-sample call: major subtype, minor subtype (co-infection), overall_sample_call
- Generate review flag: highlights samples needing manual inspection (conflicting mapping vs de novo, rescue triggered, low coverage, low mapped reads, co-infection, etc.)
- Produce Summary.csv: one row per sample, all metrics
- MultiQC report: aggregated QC metrics across all samples (Fastp, Bowtie2, Samtools, custom HCVtyper sections)

---

## Styling notes for the diagram

- Use a clean left-to-right or top-to-bottom flow
- Box 1, 3, 4, 5, 6: single-column full-width boxes
- Boxes 2a and 2b: side by side, same width, connected by a fork from Box 1 and a merge into Box 3
- Suggested color scheme:
  - Box 1 (read processing): blue/teal
  - Box 2a (mapping): orange
  - Box 2b (de novo/BLAST): green
  - Box 3 (candidate selection): purple — this is the intellectual core of the pipeline
  - Box 4 (joint mapping): orange
  - Box 5 (GLUE): red/salmon
  - Box 6 (summary): grey/neutral
- Key data flowing between boxes can be shown as labeled arrows:
  - Box 1 → split: "HCV-classified reads"
  - Box 2a → Box 3: "Candidate references + neutral read counts"
  - Box 2b → Box 3: "Assembly support per genotype"
  - Box 3 → Box 4: "Final candidate FASTAs"
  - Box 4 → Box 5: "Per-candidate consensus + BAM"
  - Box 5 → Box 6: "GLUE resistance reports"
  - Box 4 → Box 6: "Coverage stats, read counts, mapping metrics"
