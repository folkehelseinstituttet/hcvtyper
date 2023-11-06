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
include { INPUT_CHECK                    } from '../subworkflows/local/input_check'
include { MAJOR_MAPPING                  } from '../subworkflows/local/major_mapping'
include { MAJOR_MAPPING as MINOR_MAPPING } from '../subworkflows/local/major_mapping'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { BOWTIE2_BUILD                      } from '../modules/nf-core/bowtie2/build/main'
include { BLAST_MAKEBLASTDB                  } from '../modules/nf-core/blast/makeblastdb/main'
include { FASTQC                             } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM              } from '../modules/nf-core/fastqc/main'
include { CUTADAPT                           } from '../modules/nf-core/cutadapt/main'
include { MULTIQC                            } from '../modules/nf-core/multiqc/main'
include { KRAKEN2_KRAKEN2                    } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_FOCUSED } from '../modules/nf-core/kraken2/kraken2/main'
include { SPADES                             } from '../modules/nf-core/spades/main'
include { BLAST_BLASTN                       } from '../modules/nf-core/blast/blastn/main'
include { BLASTPARSE                         } from '../modules/local/blastparse.nf'
include { BOWTIE2_ALIGN                      } from '../modules/nf-core/bowtie2/align/main' 
include { PARSEFIRSTMAPPING                  } from '../modules/local/parsefirstmapping.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS        } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO: Make subworkflows for all refs, major and minor

// Info required for completion email and summary
def multiqc_report = []

workflow VIRALSEQ {

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
    // MODULE: Run Bowtie2_build to create a reference index
    //
    BOWTIE2_BUILD(
        file(params.references)
    )
    ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)

    //
    // MODULE: Create Blast database from reference sequences
    //
    BLAST_MAKEBLASTDB (
        file(params.references)
    )
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

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

    //
    // MODULE: Run Kraken2 to classify reads
    //
    //Channel.value(file(params.kraken_all_db)).view()
    KRAKEN2_KRAKEN2 (
        CUTADAPT.out.reads,
        Channel.value(file(params.kraken_all_db)),
        false,
        false
    )
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first().ifEmpty(null))

    //
    // MODULE: Run Kraken2 to identify target viral reads
    //
    KRAKEN2_FOCUSED (
        CUTADAPT.out.reads,
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
    SPADES (
        ch_reads,
        [], // Empty input channel. Can be used to specify hmm profile
        []  // Empty input channel. Placeholder for separate speficication of reads. 
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
    ch_blastparse = BLAST_BLASTN.out.txt.join(SPADES.out.scaffolds)
    BLASTPARSE (
        ch_blastparse,
        file(params.references),
        params.agens
    )
    ch_versions = ch_versions.mix(BLASTPARSE.out.versions.first())

    //
    // MODULE: Map classified reads against all references
    //
    if (params.mapper == "bowtie2") {
        BOWTIE2_ALIGN (
            KRAKEN2_FOCUSED.out.classified_reads_fastq,
            BOWTIE2_BUILD.out.index,
            false, // Do not save unmapped reads
            true, // Sort bam file
            "first_mapping" // name for output files
        )
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())
    } 
    else if (params.mapper == "tanoti") {
        TANOTI (
            KRAKEN2_FOCUSED.out.classified_reads_fastq,
            [],
            false, // Do not save unmapped reads
            true, // Sort bam file
            "first_mapping" // name for output files
        )
        ch_versions = ch_versions.mix(TANOTI.out.versions.first())
    }

    //
    // MODULE: Identify the two references with most mapped reads
    //
    // Join idxstats and depth on the meta to ensure no sample mixup
    ch_parse = BOWTIE2_ALIGN.out.idxstats.join(BOWTIE2_ALIGN.out.depth)
    PARSEFIRSTMAPPING (
        ch_parse,
        file(params.references)
    )

    //
    // SUBWORKFLOW: Map reads against the majority reference
    //
    // Create input channel for BOWTIE2_BUILD in the major map channel
    //ch_build_major = PARSEFIRSTMAPPING.out.major_fasta.map { it[1] } // Extract the fasta file
    ch_major_mapping = PARSEFIRSTMAPPING.out.major_fasta.join(KRAKEN2_FOCUSED.out.classified_reads_fastq)
    //ch_major_mapping.view()
    MAJOR_MAPPING (
        ch_major_mapping,
        "majority"
    )
    ch_versions = ch_versions.mix(MAJOR_MAPPING.out.versions)
/*
    //
    // SUBWORKFLOW: Map reads against a potential minority reference
    //

    // Create input channel for mapping against a subset of the references    
    
    // The thing is that I need to ensure that the nr of reads mapped and the coverage 
    // comes from the same sample as the reads that will be mapped against the minority.
    // I.e., that the the meta of the values is the same as the meta of the reads.
    // Join the two channels, which will join on the meta. 
    // This will now have the structure: [ meta ], [ reads ], [ csv ]
    ch_join = CUTADAPT.out.reads.join(PARSEFIRSTMAPPING.out.csv)

    // Create a new channel with the structure tuple val(meta), path(reads)
    // The meta will contain all the elements from meta and the csv file
    ch_map_minor = ch_join
        .map { meta, reads, csv -> 
        def elements = csv.splitCsv( header: true, sep:',')
        return [meta + elements[0], reads] 
        }
    
    // This will extract the minor_reads and minor_cov keys
    ch_map_minor.map { meta, reads -> [ meta.subMap( ['minor_reads','minor_cov'] ) ] }.view()
    ch_map_minor.map { meta, reads -> [ meta.subMap( ['minor_reads'] ) ] }.view()
    println ch_map_minor.meta.minor_reads
    // if (ch_map_minor.map { meta, reads -> [ meta.subMap( ['minor_reads'] ) ] }.toInteger() > 50) {
    //     println "HEI"
    // }
    // if (ch_map_minor.collect { list -> list[2]} > params.minAgensRead & yy > cov) {
    //     MINOR_MAPPING (
    //     ch_build_major,
    //     ch_map_minor,
    //     "minority"
    // )
    // }


    

    

    // Next step: choose whether to map against the two references with most mapped reads,
    // or to map against the references with the "best" blast hits.
    // I guess I should make an input channel that can be used for the same bowtie process as above.
    // Either use Nextflow/groovy code here, or create a local process that spits out the 
    // correct channel and takes a parameter that specifies which strategy to use.
    // Should I identify minority also with de novo? The problem is that there are very often several
    // good blast hits. 
*/
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
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
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
