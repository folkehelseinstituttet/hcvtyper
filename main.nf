#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    niph/viralseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/niph/viralseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def primer_set         = ''
def primer_set_version = 0
if (params.platform == 'nanopore') {
    primer_set          = 'artic'
    primer_set_version  = params.primer_set_version
    params.artic_scheme = WorkflowMain.getGenomeAttribute(params, 'scheme', log, primer_set, primer_set_version)
}

params.fasta         = WorkflowMain.getGenomeAttribute(params, 'fasta'     , log, primer_set, primer_set_version)
params.gff           = WorkflowMain.getGenomeAttribute(params, 'gff'       , log, primer_set, primer_set_version)
params.bowtie2_index = WorkflowMain.getGenomeAttribute(params, 'bowtie2'   , log, primer_set, primer_set_version)
params.primer_bed    = WorkflowMain.getGenomeAttribute(params, 'primer_bed', log, primer_set, primer_set_version)

params.nextclade_dataset           = WorkflowMain.getGenomeAttribute(params, 'nextclade_dataset'          , log, primer_set, primer_set_version)
params.nextclade_dataset_name      = WorkflowMain.getGenomeAttribute(params, 'nextclade_dataset_name'     , log, primer_set, primer_set_version)
params.nextclade_dataset_reference = WorkflowMain.getGenomeAttribute(params, 'nextclade_dataset_reference', log, primer_set, primer_set_version)
params.nextclade_dataset_tag       = WorkflowMain.getGenomeAttribute(params, 'nextclade_dataset_tag'      , log, primer_set, primer_set_version)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp } from 'plugin/nf-validation'

// Print help message if needed
if (params.help) {
    def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
    def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
    def String command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv --genome GRCh37 -profile docker"
    log.info logo + paramsHelp(command) + citation + NfcoreTemplate.dashedLine(params.monochrome_logs)
    System.exit(0)
}

// Validate input parameters
if (params.validate_params) {
    validateParameters()
}

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.platform == 'illumina' & params.agens == 'HCV') {
    include { VIRALSEQ } from './workflows/viralseq'
} else if (params.platform == 'illumina' & params.agens == 'ROV') {
    include { ROV_ILLUMINA } from './workflows/rov_illumina'
} else if (params.platform == 'nanopore' & params.agens == 'HBV') {
    include { HBV_NANOPORE } from './workflows/hbv_nanopore'
}



//
// WORKFLOW: Run main niph/viralseq analysis pipeline
//
workflow NIPH_VIRALSEQ {

    //
    // WORKFLOW: HCV genome assembly and analysis from Illumina capture data
    //
    if (params.platform == 'illumina' & params.agens == 'HCV') {
        VIRALSEQ ()

    //
    // WORKFLOW: ROV genome assembly and analysis from Illumina data
    //
    } else if (params.platform == 'illumina' & params.agens == 'ROV') {
        ROV_ILLUMINA ()

    //
    // WORKFLOW: HBV genome assembly and analysis from Gunther PCR and Nanopore data
    //
    } else if (params.platform == 'nanopore' & params.agens == 'HBV') {
        HBV_NANOPORE ()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NIPH_VIRALSEQ ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
