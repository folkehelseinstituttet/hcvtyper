/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowViralseq.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK                       } from '../subworkflows/local/input_check'
include { BAM_MARKDUPLICATES_SAMTOOLS       } from '../subworkflows/nf-core/bam_markduplicates_samtools/main'
include { TARGETED_MAPPING as MAJOR_MAPPING } from '../subworkflows/local/targeted_mapping'
include { TARGETED_MAPPING as MINOR_MAPPING } from '../subworkflows/local/targeted_mapping'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
//include { fromSamplesheet                    } from 'plugin/nf-validation'
include { BOWTIE2_BUILD                      } from '../modules/nf-core/bowtie2/build/main'
include { BLAST_MAKEBLASTDB                  } from '../modules/nf-core/blast/makeblastdb/main'
include { FASTQC as FASTQC_RAW               } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM              } from '../modules/nf-core/fastqc/main'
include { CUTADAPT                           } from '../modules/nf-core/cutadapt/main'
include { HCV_GLUE             } from '../modules/local/hcvglue'
include { MULTIQC                            } from '../modules/nf-core/multiqc/main'
include { KRAKEN2_KRAKEN2                    } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_FOCUSED } from '../modules/nf-core/kraken2/kraken2/main'
include { SPADES                             } from '../modules/nf-core/spades/main'
include { BLAST_BLASTN                       } from '../modules/nf-core/blast/blastn/main'
include { BOWTIE2_ALIGN                      } from '../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_WITHDUP } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MARKDUP } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_IDXSTATS as SAMTOOLS_IDXSTATS_WITHDUP } from '../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_IDXSTATS as SAMTOOLS_IDXSTATS_MARKDUP } from '../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_DEPTH as SAMTOOLS_DEPTH_WITHDUP } from '../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_DEPTH as SAMTOOLS_DEPTH_MARKDUP } from '../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_STATS                     } from '../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_SORMADUP                  } from '../modules/nf-core/samtools/sormadup/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS        } from '../modules/nf-core/custom/dumpsoftwareversions/main'

//
// Local modules
//
include { INSTRUMENT_ID                       } from '../modules/local/instrument_id'
include { BLASTPARSE                          } from '../modules/local/blastparse.nf'
include { TANOTI_ALIGN                        } from '../modules/local/tanoti.nf'
include { PARSEFIRSTMAPPING                   } from '../modules/local/parsefirstmapping.nf'
include { GLUEPARSE as HCV_GLUE_PARSER_MAJOR  } from '../modules/local/glueparse'
include { GLUEPARSE as HCV_GLUE_PARSER_MINOR  } from '../modules/local/glueparse'
include { SUMMARIZE_HCV as SUMMARIZE          } from '../modules/local/summarize_hcv'
include { SORT_IDXSTATS                       } from '../modules/local/idxstats_sort.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO: Make subworkflows for all refs, major and minor

// Info required for completion email and summary
def multiqc_report = []

workflow HCV_ILLUMINA {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // MODULE: Identify the instrument ID
    //
    INSTRUMENT_ID (
        INPUT_CHECK.out.reads
    )

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
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())

    //
    // MODULE: Trim reads with Cutadapt
    //
    CUTADAPT(
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(CUTADAPT.out.versions)


    //
    // MODULE: Run FastQC on trimmed reads
    //
    FASTQC_TRIM (
        CUTADAPT.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC_TRIM.out.versions.first())

    // NOTE:
    // In some cases there are empty fastq files after trimming. Remove these before Kraken2
    ch_kraken = CUTADAPT.out.reads
        .map { meta, fastq ->
            n = fastq[0].countFastq() // Count fastq reads in the R1 fastq file
            return [meta, fastq, n] // Add the count as the last element in the tuple
        }
        .filter { n > 1 } // Filter out empty fastq files
        .map { meta, fastq, n -> [meta, fastq] } // Return the count to get the channel structure correct for BBMAP_BBNORM

    //
    // MODULE: Run Kraken2 to classify reads
    //
    KRAKEN2_KRAKEN2 (
        ch_kraken,
        Channel.value(file(params.kraken_all_db)),
        false,
        false
    )
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first().ifEmpty(null))

    //
    // MODULE: Run Kraken2 to identify target viral reads
    //
    KRAKEN2_FOCUSED (
        ch_kraken,
        Channel.value(file(params.kraken_focused)),
        params.save_output_fastqs,
        params.save_reads_assignment
    )

    //
    // MODULE: Run Spades to assemble classified reads
    //

    // Create input read channel for SPADES.
    // A tuple with meta, paired Illumina reads, and empty elements for pacbio and nanopore reads
    ch_reads = KRAKEN2_FOCUSED.out.classified_reads_fastq.map { meta, fastq -> [ meta, fastq, [], [] ] }
    if (!params.skip_assembly) {
            SPADES (
                ch_reads,
                [], // Empty input channel. Can be used to specify hmm profile
                []  // Empty input channel. Placeholder for separate specification of reads.
            )
            ch_versions = ch_versions.mix(SPADES.out.versions.first())

            //
            // MODULE: Blast assembled scaffolds against viral references.
            //
            BLAST_BLASTN (
                SPADES.out.scaffolds,
                BLAST_MAKEBLASTDB.out.db
            )
            ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())

            //
            // MODULE: Parse blast output
            //
            // Create input channel that holds val(meta), path(blast_out), path(scaffolds)
            // TODO: Decide the major genotype present (I guess the best blast hit, maybe with a few criteria) and the potential minor.
            // Then create an output channel with this info that can be put into MAJOR_MAPPING and MINOR_MAPPING.
            // See Issue on this
            ch_blastparse = BLAST_BLASTN.out.txt.join(SPADES.out.scaffolds)
            BLASTPARSE (
                ch_blastparse,
                file(params.references),
                params.agens
            )
            ch_versions = ch_versions.mix(BLASTPARSE.out.versions.first())
    }

    //
    // MODULE: Map classified reads against all references
    //
    if (params.mapper == "bowtie2") {
        BOWTIE2_ALIGN (
            KRAKEN2_FOCUSED.out.classified_reads_fastq,
            BOWTIE2_BUILD.out.index,
            false, // Do not save unmapped reads
            true // Sort bam file
        )
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())
        ch_aligned = BOWTIE2_ALIGN.out.aligned
    }
    else if (params.mapper == "tanoti") {
        TANOTI_ALIGN (
            KRAKEN2_FOCUSED.out.classified_reads_fastq,
            [ [], file(params.references) ], // Add empty meta map before the reference file path
            true, // Sort bam file
            params.tanoti_stringency_1
        )
        ch_versions = ch_versions.mix(TANOTI_ALIGN.out.versions.first())
        ch_aligned = TANOTI_ALIGN.out.aligned
    }

    SAMTOOLS_INDEX_WITHDUP (
        ch_aligned
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_WITHDUP.out.versions.first())

    SAMTOOLS_IDXSTATS_WITHDUP (
        ch_aligned.join(SAMTOOLS_INDEX_WITHDUP.out.bai) // val(meta), path(bam), path(bai)
    )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS_WITHDUP.out.versions.first())

    SORT_IDXSTATS (
        SAMTOOLS_IDXSTATS_WITHDUP.out.idxstats
    )
    ch_versions = ch_versions.mix(SORT_IDXSTATS.out.versions)

    SAMTOOLS_DEPTH_WITHDUP (
        ch_aligned,
        [ [], []] // Passing empty channels instead of an interval file
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH_WITHDUP.out.versions.first())

    //
    // MODULE: Identify the two references with most mapped reads, duplicates included
    //
    ch_parsefirstmapping = SAMTOOLS_IDXSTATS_WITHDUP.out.idxstats.join(SAMTOOLS_DEPTH_WITHDUP.out.tsv) // val(meta), path(idxstats), path(tsv)
        .filter { meta, idxstats, tsv ->
            tsv.size() > 0 // Filter out empty tsv files
        }

    PARSEFIRSTMAPPING (
        // Join idxstats and depth on the meta map
        ch_parsefirstmapping,
        file(params.references)
    )

    // Remove duplicate reads
    SAMTOOLS_SORMADUP (
        ch_aligned,
        [ [], file(params.references) ]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORMADUP.out.versions.first())

    SAMTOOLS_INDEX_MARKDUP (
        SAMTOOLS_SORMADUP.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_MARKDUP.out.versions.first())

    SAMTOOLS_IDXSTATS_MARKDUP (
        SAMTOOLS_SORMADUP.out.bam.join(SAMTOOLS_INDEX_MARKDUP.out.bai) // val(meta), path(bam), path(bai)
    )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS_MARKDUP.out.versions.first())

    SAMTOOLS_DEPTH_MARKDUP (
        SAMTOOLS_SORMADUP.out.bam,
        [ [], []] // Passing empty channels instead of an interval file
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH_MARKDUP.out.versions.first())

    SAMTOOLS_STATS (
        SAMTOOLS_SORMADUP.out.bam.join(SAMTOOLS_INDEX_MARKDUP.out.bai), // val(meta), path(bam), path(bai)
        Channel.value(file(params.references)).map { [ [:], it ] } // Add empty meta map before the reference file path
    )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    //
    // SUBWORKFLOW: Map reads against the majority reference
    //
    if (params.strategy == "mapping") {
        ch_major_mapping = PARSEFIRSTMAPPING.out.major_fasta.join(KRAKEN2_FOCUSED.out.classified_reads_fastq)
    } else if (params.strategy == "denovo") {
        ch_major_mapping = BLASTPARSE.out.major_fasta.join(KRAKEN2_FOCUSED.out.classified_reads_fastq)
    }

    MAJOR_MAPPING(
        ch_major_mapping, // val(meta), path(fasta), path(reads)
    )
    ch_versions = ch_versions.mix(MAJOR_MAPPING.out.versions)

    //
    // SUBWORKFLOW: Map reads against a potential minority reference
    //

    // Create input channel for mapping against a subset of the references
    if (params.strategy == "mapping") {
        ch_join = PARSEFIRSTMAPPING.out.minor_fasta.join(KRAKEN2_FOCUSED.out.classified_reads_fastq) // meta, fasta, reads
        ch_join_2 = ch_join.join(PARSEFIRSTMAPPING.out.csv) // meta, fasta, reads, csv

        // Create a new channel with the structure tuple val(meta), path(fasta), path(reads)
        // The meta will contain all the elements from meta and the csv file. meta, reads
        ch_map_minor = ch_join_2
            .map { meta, fasta, reads, csv ->
            def elements = csv.splitCsv( header: true, sep:',')
            return [meta + elements[0], fasta, reads]
            }

        // Filter on read nr and coverage
        // This will result in a channel with values that meet the read nr and coverage criteria
        ch_map_minor_filtered = ch_map_minor
        .filter { entry ->
            def mappedReads = entry[0]['total_mapped_reads'].toInteger()
            def minorCov = entry[0]['minor_cov'].toInteger()
            mappedReads > params.minAgensRead && minorCov > params.minAgensCov
        }
    } else if (params.strategy == "denovo") {
        ch_join = BLASTPARSE.out.minor_fasta.join(KRAKEN2_FOCUSED.out.classified_reads_fastq) // meta, fasta, reads
        ch_join_2 = ch_join.join(BLASTPARSE.out.csv) // meta, fasta, reads, csv

        // Create a new channel with the structure tuple val(meta), path(reads)
        // The meta will contain all the elements from meta and the csv file. meta, reads
        ch_map_minor = ch_join_2
            .map { meta, fasta, reads, csv ->
            def elements = csv.splitCsv( header: true, sep:',')
            return [meta + elements[0], fasta, reads]
            }

        // Filter on read nr and coverage
        // This will result in a channel with values that meet the read nr and coverage criteria
        ch_map_minor_filtered = ch_map_minor
        .filter { entry ->
            def minorLength = entry[0]['minor_length'].toInteger()
            minorLength > params.minDenovoLength
        }
    }

    MINOR_MAPPING (
        ch_map_minor_filtered, // val(meta), path(fasta), path(reads)
    )

    //
    // MODULE: Run GLUE genotyping and resistance annotation for HCV
    //
    if (params.agens == "HCV" && !params.skip_hcvglue) {
        HCV_GLUE (
            MAJOR_MAPPING.out.aligned.collect({it[1]}).mix(MINOR_MAPPING.out.aligned.collect({it[1]}))
        )
        ch_versions = ch_versions.mix(HCV_GLUE.out.versions)

        //HCV_GLUE_MINOR (
        //    MINOR_MAPPING.out.aligned.collect{ it[1] }
        //)
        //ch_versions = ch_versions.mix(HCV_GLUE_MINOR.out.versions)

        // Collect all glue reports and parse them
        HCV_GLUE_PARSER (
            HCV_GLUE.out.GLUE_json.collect{ it[1] },
            "major"
        )
        ch_versions = ch_versions.mix(HCV_GLUE_PARSER_MAJOR.out.versions)

        //HCV_GLUE_PARSER_MINOR (
        //    HCV_GLUE_MINOR.out.GLUE_json.collect{ it[1] },
        //    "minor"
       // )
        //ch_versions = ch_versions.mix(HCV_GLUE_PARSER_MINOR.out.versions)
    }

    //
    // MODULE: Summarize
    //
    // Create channel with this structure: path(stats), path(depth), path(blast), path(json)
    // Collect all the files in separate channels for clarixty. Don't need the meta
    ch_sequence_id      = INSTRUMENT_ID.out.id.collect({it[1]})
    ch_cutadapt         = CUTADAPT.out.log.collect({it[1]})
    ch_classified_reads = KRAKEN2_FOCUSED.out.report.collect({it[1]})
    ch_stats_withdup    = MAJOR_MAPPING.out.stats_withdup.collect({it[1]}).mix(MINOR_MAPPING.out.stats_withdup.collect({it[1]}))
    ch_stats_markdup    = MAJOR_MAPPING.out.stats_markdup.collect({it[1]}).mix(MINOR_MAPPING.out.stats_markdup.collect({it[1]}))
    ch_depth            = MAJOR_MAPPING.out.depth.collect({it[1]}).mix(MINOR_MAPPING.out.depth.collect({it[1]}))
    if (!params.skip_assembly) {
        ch_blast = BLAST_BLASTN.out.txt.collect({it[1]})
    } else {
        ch_blast = file("dummy_file")
    }
    if (params.agens == "HCV" && !params.skip_hcvglue) {
        ch_glue = HCV_GLUE_PARSER.out.GLUE_summary
    } else {
        ch_glue = file("dummy_file")
    }

    SUMMARIZE (
        params.tanoti_stringency_1,
        params.tanoti_stringency_2,
        ch_cutadapt.collect(),
        ch_classified_reads.collect(),
        ch_stats_withdup.collect(),
        ch_stats_markdup.collect(),
        ch_depth.collect(),
        ch_blast,
        ch_glue,
        ch_sequence_id.collect()
    )

    //
    // MODULE: Dump software versions
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowViralseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowViralseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collectFile(name: 'software_mqc_versions.yml'))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT.out.log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIM.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_KRAKEN2.out.report.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(SUMMARIZE.out.mqc.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )
    multiqc_report = MULTIQC.out.report.toList()

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
