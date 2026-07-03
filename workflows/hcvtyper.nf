/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { softwareVersionsToYAML                         } from '../subworkflows/nf-core/utils_nfcore_pipeline/main.nf'
include { GET_MAPPING_STATS as GET_MAPPING_STATS_WITHDUP } from '../subworkflows/local/get_mapping_stats'
include { GET_MAPPING_STATS as GET_MAPPING_STATS_MARKDUP } from '../subworkflows/local/get_mapping_stats'
include { JOINT_MAPPING                                  } from '../subworkflows/local/joint_mapping/main'
include { CONTAMINATION_CHECK                            } from '../subworkflows/local/contamination_check/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
//include { fromSamplesheet                    } from 'plugin/nf-validation'
include { INPUT_CHECK                        } from '../subworkflows/local/input_check'
include { BOWTIE2_BUILD                      } from '../modules/nf-core/bowtie2/build/main'
include { BLAST_MAKEBLASTDB                  } from '../modules/nf-core/blast/makeblastdb/main'
include { FASTP                              } from '../modules/nf-core/fastp/main'
include { FASTQC as FASTQC_RAW               } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM              } from '../modules/nf-core/fastqc/main'
include { CUTADAPT                           } from '../modules/nf-core/cutadapt/main'
include { PRINSEQPLUSPLUS                    } from '../modules/nf-core/prinseqplusplus/main'
include { MULTIQC                            } from '../modules/nf-core/multiqc/main'
include { KRAKEN2_KRAKEN2                    } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_FOCUSED } from '../modules/nf-core/kraken2/kraken2/main'
include { SPADES                             } from '../modules/nf-core/spades/main'
include { BLAST_BLASTN                       } from '../modules/nf-core/blast/blastn/main'
include { BOWTIE2_ALIGN                      } from '../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_SORMADUP                  } from '../modules/nf-core/samtools/sormadup/main'
include { PICARD_COLLECTINSERTSIZEMETRICS    } from '../modules/nf-core/picard/collectinsertsizemetrics/main'
include { UNTAR as UNTAR_KRAKEN_ALL          } from '../modules/nf-core/untar/main'
include { UNTAR as UNTAR_KRAKEN_FOCUSED      } from '../modules/nf-core/untar/main'

//
// Local modules
//
include { INSTRUMENTID                       } from '../modules/local/instrumentid/main'
include { BLASTPARSE                         } from '../modules/local/blastparse/main'
include { PARSEFIRSTMAPPING                  } from '../modules/local/parsefirstmapping/main'
include { RESCUE_EVALUATION                  } from '../modules/local/rescueevaluation/main'
include { GLUEPARSE as HCV_GLUE_PARSER       } from '../modules/local/glueparse/main'
include { HCVGLUE                            } from '../modules/local/hcvglue/main'
include { SUMMARIZE                          } from '../modules/local/summarize/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HCVTYPER {

    main:

    ch_versions = Channel.empty()

    //
    // CONFIG FILES
    //
    ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    reads = INPUT_CHECK.out.reads
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // Prepare Kraken2 database for all domains of life
    //
    ch_kraken_all_db = Channel.empty()
    if (params.kraken_all_db) {
        if (params.kraken_all_db.endsWith('.tar.gz')) {
            UNTAR_KRAKEN_ALL (
                [ [:], params.kraken_all_db ] // Add empty meta map
            )
            ch_kraken_all_db = UNTAR_KRAKEN_ALL.out.untar.map { it[1] } // Do not extract the meta map which is emitted by default
            ch_versions = ch_versions.mix(UNTAR_KRAKEN_ALL.out.versions.first())
        } else {
            ch_kraken_all_db = Channel.value(file(params.kraken_all_db))
        }
    }

    //
    // Prepare Kraken2 database for HCV only
    //
    ch_kraken_focused_db = Channel.empty()
        if (params.kraken_focused_db.endsWith('.tar.gz')) {
            UNTAR_KRAKEN_FOCUSED (
                [ [:], params.kraken_focused_db ] // Add empty meta map
            )
            ch_kraken_focused_db = UNTAR_KRAKEN_FOCUSED.out.untar.map { it[1] } // Do not extract the meta map which is emitted by default
            ch_versions = ch_versions.mix(UNTAR_KRAKEN_FOCUSED.out.versions.first())
        } else {
            ch_kraken_focused_db = Channel.value(file(params.kraken_focused_db))
        }

    //
    // MODULE: Identify the instrument ID
    //
    INSTRUMENTID (
        reads
    )
    ch_versions = ch_versions.mix(INSTRUMENTID.out.versions.first())

    //
    // MODULE: Run Bowtie2_build to create a reference index
    //

    BOWTIE2_BUILD(
        [ [], file(params.references) ] // Add empty meta map before the reference file path
    )
    ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)

    //
    // MODULE: Create Blast database from reference sequences
    //
    BLAST_MAKEBLASTDB (
        [ [id:"blastdb/"], file(params.references) ] // Add meta map before the reference file path
    )
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC_RAW (
        reads
    )
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions)

    //
    // MODULE: Trim reads with either Cutadapt or Fastp
    //
    if (params.trimmer == "fastp") {
        FASTP(
            reads
                .map { meta, fastq ->
                    def n = fastq[0].countFastq() // Count fastq reads in the R1 fastq file
                    return [meta, fastq, n] // Add the count as the last element in the tuple
                }
                .filter { _meta, _fastq, n -> n > 0 } // Filter out empty fastq files because Fastp does not produce correct gzipped file if input fastq is empty
                .map { meta, fastq, _n -> [ meta, fastq, [] ] }, // Remove the count and add empty element to get the structure correct for Fastp
            false,
            false,
            false
        )
        ch_trimmed_reads = FASTP.out.reads
        ch_versions = ch_versions.mix(FASTP.out.versions)
    } else if (params.trimmer == "cutadapt") {
        CUTADAPT(
            reads
        )
        ch_trimmed_reads = CUTADAPT.out.reads
        ch_versions = ch_versions.mix(CUTADAPT.out.versions)
    }

    //
    // MODULE: Run FastQC on trimmed reads
    //
    FASTQC_TRIM (
        ch_trimmed_reads
    )
    ch_versions = ch_versions.mix(FASTQC_TRIM.out.versions.first())

    //
    // MODULE: Remove low complexity reads with Prinseq++ (optional)
    //

    // Initialise a default empty channel
    ch_prinseq_log = Channel.empty()
    if (params.filter_low_complexity) {
        // NOTE:
        // In some cases there are empty fastq files after trimming. Remove these before PRINSEQPLUSPLUS
        ch_prinseq = ch_trimmed_reads
            .map { meta, fastq ->
                def n = fastq[0].countFastq() // Count fastq reads in the R1 fastq file
                return [meta, fastq, n] // Add the count as the last element in the tuple
            }
            .filter { _meta, _fastq, n -> n > 0 } // Filter out empty fastq files
            .map { meta, fastq, _n -> [meta, fastq] } // Remove the count to get the channel structure correct for PRINSEQPLUSPLUS

        PRINSEQPLUSPLUS (
            ch_prinseq
        )
        ch_versions = ch_versions.mix(PRINSEQPLUSPLUS.out.versions.first())

        // Kraken2 input is from PRINSEQPLUSPLUS
        ch_kraken_input = PRINSEQPLUSPLUS.out.good_reads

        // For MultiQC later
        ch_prinseq_log = PRINSEQPLUSPLUS.out.log.collect { it[1] }.ifEmpty([])
    } else {
        // Skip PRINSEQPLUSPLUS and use trimmed reads as input to Kraken2
        ch_kraken_input = ch_trimmed_reads
    }

    //
    // MODULE: Run Kraken2 to classify reads
    //
    KRAKEN2_KRAKEN2 (
        ch_kraken_input,
        ch_kraken_all_db,
        false,
        false
    )
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first())

    //
    // MODULE: Run Kraken2 to identify target viral reads
    //
    KRAKEN2_FOCUSED (
        ch_kraken_input,
        ch_kraken_focused_db,
        params.save_output_fastqs,
        params.save_reads_assignment
    )

    //
    // MODULE: Run Spades to assemble classified reads
    //

    // Create input read channel for SPADES.
    // A tuple with meta, paired Illumina reads, and empty elements for pacbio and nanopore reads
    // Filter in case of no classified reads emitted by KRAKEN2_FOCUSED
    ch_reads = KRAKEN2_FOCUSED.out.classified_reads_fastq
        .map { meta, fastq ->
            def n = fastq[0].countFastq() // Count fastq reads in the R1 fastq file
            return [meta, fastq, n] // Add the count as the last element in the tuple
        }
        .filter { _meta, _fastq, n -> n > 0 } // Filter out empty fastq files
        .map { meta, fastq, _n -> [ meta, fastq, [], [] ] } // Recreate the channel structure correct for SPADES

    SPADES (
        ch_reads,
        [], // Empty input channel. Can be used to specify hmm profile
        []  // Empty input channel. Placeholder for separate specification of reads.
    )
    ch_versions = ch_versions.mix(SPADES.out.versions.first())

    //
    // MODULE: Blast assembled contigs against viral references.
    //
    // NOTE:
    // In some cases there is an empty contig file produced by Spades. Filter out these
    ch_blastn = SPADES.out.contigs
        .map { meta, contigs ->
        def n = contigs.countFasta() // Count fasta records
        return [meta, contigs, n] // Add the count as the last element in the tuple
    }
    .filter { _meta, _contigs, n -> n > 0 } // Filter out empty fasta files
    .map { meta, contigs, _n -> [meta, contigs] } // Return the count to get the channel structure correct for BLASTN_BLASTN
    BLAST_BLASTN (
        ch_blastn,
        BLAST_MAKEBLASTDB.out.db,
        [], // taxidlist - empty, no taxonomic filtering
        "", // taxids - empty string, no taxonomic filtering
        false // negative_tax - false, not using negative filtering
    )
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())

    //
    // MODULE: Parse blast output
    //
    ch_blastparse = BLAST_BLASTN.out.txt.join(SPADES.out.contigs) // Create input channel that holds val(meta), path(blast_out), path(contigs)
    BLASTPARSE (
        ch_blastparse,
        file(params.references),
        params.agens
    )
    ch_versions = ch_versions.mix(BLASTPARSE.out.versions.first())

    //
    // SUBWORKFLOW: Detect cross-sample contamination via all-vs-all BLAST
    //
    if (!params.skip_contamination_check) {
        CONTAMINATION_CHECK(
            SPADES.out.contigs,
            Channel.empty(),  // fastp JSONs — not wired in main pipeline
            Channel.empty()   // GLUE JSONs  — not wired in main pipeline
        )
        ch_versions = ch_versions.mix(CONTAMINATION_CHECK.out.versions)
    }

    //
    // MODULE: Map classified reads against all references
    //
    BOWTIE2_ALIGN (
        KRAKEN2_FOCUSED.out.classified_reads_fastq,
        BOWTIE2_BUILD.out.index,
        [ [], file(params.references) ], // Add empty meta map and reference fasta for CRAM support
        false, // Do not save unmapped reads
        true // Sort bam file
    )
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())
    ch_aligned = BOWTIE2_ALIGN.out.bam

    //
    // SUBWORKFLOW: Get mapping statistics with duplicates included
    //
    GET_MAPPING_STATS_WITHDUP (
        ch_aligned, // Channel: [ val(meta), [ bam ] ]
        Channel.value(file(params.references)).map { [ [:], it ] } // Add empty meta map before the reference file path
    )
    versions = GET_MAPPING_STATS_WITHDUP.out.versions

    // Remove duplicate reads
    SAMTOOLS_SORMADUP (
        ch_aligned,
        [ [], file(params.references) ]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORMADUP.out.versions.first())

    //
    // MODULE: Collect insert size metrics with duplicates removed
    //
    PICARD_COLLECTINSERTSIZEMETRICS (
        SAMTOOLS_SORMADUP.out.bam
    )
    ch_versions = ch_versions.mix(PICARD_COLLECTINSERTSIZEMETRICS.out.versions.first())

    //
    // SUBWORKFLOW: Get mapping statistics with duplicates removed
    //
    GET_MAPPING_STATS_MARKDUP (
        SAMTOOLS_SORMADUP.out.bam,
        Channel.value(file(params.references)).map { [ [:], it ] } // Add empty meta map before the reference file path
    )
    versions = GET_MAPPING_STATS_MARKDUP.out.versions

    //
    // MODULE: Identify the two references with most mapped reads, duplicates included
    //
    ch_parsefirstmapping = GET_MAPPING_STATS_WITHDUP.out.idxstats.join(GET_MAPPING_STATS_WITHDUP.out.tsv) // val(meta), path(idxstats), path(tsv)
        .filter { _meta, _idxstats, tsv ->
            tsv.size() > 0 // Filter out empty tsv files
        }

    PARSEFIRSTMAPPING (
        // Join idxstats and depth on the meta map
        ch_parsefirstmapping,
        file(params.references),
        file("${projectDir}/bin/genotype_utils.R", checkIfExists: true)
    )

    //
    // MODULE: De-novo subtype rescue (Phase 10, denovo-subtype-rescue, D-09..D-12)
    //
    // RESCUE_EVALUATION sits BETWEEN the de-novo BLAST evidence (BLASTPARSE) and the
    // candidate-mapping builder. It may REPLACE a mapped candidate's reference with the
    // de-novo-derived reference when the mapped subtype disagrees with the de-novo top hit
    // and the contig clears the four quality floors, and forces that candidate's
    // confirmation_status to 'pass' so the builder routes it to JOINT_MAPPING.
    //
    // De novo assembly always runs, so BLASTPARSE.out.support is always defined.
    ch_blastparse_support = BLASTPARSE.out.support

    // Build the RESCUE_EVALUATION input tuple (meta, candidates_csv, support_csv,
    // cand_fastas) by joining on meta.id. PARSEFIRSTMAPPING.out.candidate_fasta
    // is tuple(meta, parsefirstmapping_csv, cand_fastas) -- extract cand_fastas. The
    // support leg uses remainder:true (D-10) so a sample with no de-novo contig
    // (empty BLASTPARSE support) does not drop samples; the R script's typed-empty
    // guard handles the missing file.
    ch_rescue_input = PARSEFIRSTMAPPING.out.candidates
        .join(PARSEFIRSTMAPPING.out.candidate_fasta, remainder: true)       // meta, candidates_csv, parsefirstmapping_csv?, cand_fastas?
        .join(ch_blastparse_support, remainder: true)                       // ..., support_csv?
        .map { meta, candidates_csv, _parsefirstmapping_csv, cand_fastas, support_csv ->
            // remainder:true fills absent legs with null. The module's path() inputs accept []
            // for a missing optional file; normalize null -> [] so staging never NPEs.
            tuple(
                meta,
                candidates_csv,
                support_csv ?: [],
                cand_fastas ?: []
            )
        }

    RESCUE_EVALUATION (
        ch_rescue_input,
        file(params.references)
    )
    ch_versions = ch_versions.mix(RESCUE_EVALUATION.out.versions.first())

    //
    // SUBWORKFLOW: Competitive joint mapping (D-01/D-02/D-03 — replaces TARGETED_MAPPING)
    //
    // JOINT_MAPPING takes ONE element per SAMPLE — NOT a per-candidate flatMap. The
    // per-candidate fan-out now happens INSIDE JOINT_MAPPING, AFTER the combined dedup BAM is
    // split by reference (D-03). So here we build a simple per-sample tuple:
    //   tuple(meta, cand_fastas_list, classified_reads, candidates_csv)
    // and pass the full collected candidate-FASTA list directly (no splitCsv / no
    // confirmation_status filter at this level — that all moves into the subworkflow).
    //
    // Join the per-sample inputs by meta.id: the RESCUE_EVALUATION candidates CSV, the
    // collected candidate FASTAs (RESCUE_EVALUATION.out.candidate_fasta), and the classified
    // reads. candidate_fasta is `optional: true` (a no-candidate sample emits nothing, a
    // single-candidate sample emits only `_cand1.fa`), so its join uses `remainder: true` —
    // otherwise a missing optional emit would silently DROP the whole sample. remainder:true
    // pads the absent side with a single null (not a tuple), giving a VARIABLE-arity join
    // output; normalize to a FIXED 2-tuple (meta, fastas_or_empty) in a .map first so the
    // downstream shape is stable regardless of the optional emit. A single bare FASTA is
    // normalized to a list so the subworkflow's rank-indexed lookup is uniform (mirrors the
    // legacy line-472 normalization). candidate_rank stays a STRING throughout — never
    // .toInteger() in channel logic (NA would crash, Pitfall 4).
    ch_joint_mapping = RESCUE_EVALUATION.out.candidates
        .join(RESCUE_EVALUATION.out.candidate_fasta, remainder: true)      // meta, candidates_csv, fasta_list?
        .map { tup ->
            // tup = [meta, candidates_csv, fasta_list?]. remainder:true gives [meta, csv, null]
            // when no FASTA matched; a matched item gives [meta, csv, fasta_list].
            def meta           = tup[0]
            def candidates_csv = tup[1]
            def fasta_list     = (tup.size() > 2) ? tup[2] : null
            // Normalize: no-candidate sample -> [] ; single bare FASTA -> [fasta] ; list -> as-is.
            def fastas = (fasta_list == null) ? [] : (fasta_list instanceof List ? fasta_list : [fasta_list])
            tuple(meta, candidates_csv, fastas)
        }
        .join(KRAKEN2_FOCUSED.out.classified_reads_fastq)                  // meta, candidates_csv, fastas, classified_reads
        .map { meta, candidates_csv, fastas, classified_reads ->
            // JOINT_MAPPING take: tuple(meta, cand_fastas, reads, candidates_csv)
            tuple(meta, fastas, classified_reads, candidates_csv)
        }

    JOINT_MAPPING(
        ch_joint_mapping, // val(meta), path(cand_fastas), path(reads), path(candidates_csv)
    )
    ch_versions = ch_versions.mix(JOINT_MAPPING.out.versions)

    //
    // MODULE: Run GLUE genotyping and resistance annotation for HCV
    //
    if (!params.skip_hcvglue) {
        ch_glue_bams = JOINT_MAPPING.out.aligned
            .filter { meta, _bam -> (meta.candidate_nodup_reads ?: 0) >= params.glue_min_reads }
            .collect({ it[1] })
            .collect()

        HCVGLUE (
            ch_glue_bams,
            params.hcvglue_threshold
        )
        ch_versions = ch_versions.mix(HCVGLUE.out.versions)

        // Collect all glue reports and parse them
        HCV_GLUE_PARSER (
            HCVGLUE.out.GLUE_json.collect()
        )
        ch_versions = ch_versions.mix(HCV_GLUE_PARSER.out.versions)
    }

    //
    // MODULE: Summarize
    //
    // Create channel with this structure: path(stats), path(depth), path(blast), path(json)
    // Collect all the files in separate channels for clarity. Don't need the meta
    ch_sequence_id      = INSTRUMENTID.out.id.collect({it[1]})
    if (params.trimmer == "fastp") {
        ch_trimmed_reads         = FASTP.out.log.collect({it[1]})
    } else if (params.trimmer == "cutadapt") {
        ch_trimmed_reads         = CUTADAPT.out.log.collect({it[1]})
    }
    ch_classified_reads = KRAKEN2_FOCUSED.out.report.collect({it[1]})
    // Phase 7 (ASUP-02): stage the Phase-6 long-format *.candidates.csv alongside
    // the legacy *.parsefirstmapping.csv into parsefirst_mapping/ so summarize.R
    // can read it for the genotype-level assembly-support join.
    // D-08: stage the RESCUE_EVALUATION candidates CSV (carrying rescued_from /
    // rescue_trigger) instead of the raw PARSEFIRSTMAPPING candidates, so summarize.R
    // reads the rescue audit columns. The legacy *.parsefirstmapping.csv leg is unchanged.
    ch_summarize_first_mapping = PARSEFIRSTMAPPING.out.csv.collect({it[1]}).mix(RESCUE_EVALUATION.out.candidates.collect({it[1]})).collect()
    // Phase 11 (JMAP-03, D-07/D-09/D-15): read counts now come from SAMTOOLS_IDXSTATS on the
    // COMBINED BAM (pre- and post-dedup), ONE idxstats file per sample carrying every candidate
    // reference as a row. SAMTOOLS_STATS is removed. The variable names ch_stats_withdup /
    // ch_stats_markdup are KEPT so the SUMMARIZE call signature is unchanged; their content is
    // now idxstats TSV (.withdup.idxstats / .nodup.idxstats), staged into stats_withdup/ and
    // stats_markdup/. Plan 03 migrates summarize.R's two STATS loops to the idxstats format.
    ch_stats_withdup    = JOINT_MAPPING.out.idxstats_withdup.collect({it[1]})
    ch_stats_markdup    = JOINT_MAPPING.out.idxstats_nodup.collect({it[1]})
    ch_depth            = JOINT_MAPPING.out.depth.collect({it[1]})
    // De novo / BLAST evidence (PLUMB-01/PLUMB-02): collect the parsed BLASTPARSE
    // CSVs (*.blastparse.csv) and the per-contig table (*_blast_out.csv) into one
    // staged channel. De novo assembly always runs, so BLASTPARSE.out is always
    // defined; .ifEmpty([]) still handles a batch where every sample produced no
    // contig (empty denovo/ staging dir -> NA de novo columns + no dropped rows).
    // Phase 7 (ASUP-02): also stage the per-subtype *.assembly_support.csv into
    // denovo/ so summarize.R can join it to candidates at genotype level.
    ch_denovo = BLASTPARSE.out.csv.collect({it[1]}).mix(BLASTPARSE.out.blast_res.collect({it[1]})).mix(BLASTPARSE.out.support.collect({it[1]})).collect().ifEmpty([])
    if (params.agens == "HCV" && !params.skip_hcvglue) {
        ch_glue = HCV_GLUE_PARSER.out.GLUE_summary
    } else {
        ch_glue = []
    }
    ch_variation = JOINT_MAPPING.out.variation.collect()
    ch_consensus_distance = JOINT_MAPPING.out.consensus_distance.collect({it[1]})

    SUMMARIZE (
        workflow.manifest.version,
        workflow.manifest.name,
        file(params.input),
        ch_trimmed_reads.collect(),
        ch_classified_reads.collect(),
        ch_summarize_first_mapping,
        ch_stats_withdup.collect(),
        ch_stats_markdup.collect(),
        ch_depth.collect(),
        ch_denovo,
        ch_glue,
        ch_sequence_id.collect(),
        ch_variation.collect(),
        ch_consensus_distance.collect(),
        file("${projectDir}/bin/genotype_utils.R"),
        file("${projectDir}/bin/denovo_confirm.R"),
        file("${projectDir}/bin/denovo_layer.R"),
        file("${projectDir}/bin/assembly_support_join.R"),
        file("${projectDir}/bin/classify_roles.R"),
    )
    ch_versions = ch_versions.mix(SUMMARIZE.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: "software_mqc_versions.yml",
            sort: true,
            newLine: true
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    def summary_params = paramsSummaryMap(workflow)
    workflow_summary    = WorkflowHCVTyper.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowHCVTyper.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]))
    if (params.trimmer == "cutadapt") {
        ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT.out.log.collect{it[1]}.ifEmpty([]))
    } else if (params.trimmer == "fastp") {
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))
    }
    ch_multiqc_files = ch_multiqc_files.mix(ch_prinseq_log)
    ch_multiqc_files = ch_multiqc_files.mix(ch_trimmed_reads.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIM.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_KRAKEN2.out.report.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_FOCUSED.out.report.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(SUMMARIZE.out.mqc.collect())
    if (!params.skip_contamination_check) {
        ch_multiqc_files = ch_multiqc_files.mix(CONTAMINATION_CHECK.out.mqc)
    }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
