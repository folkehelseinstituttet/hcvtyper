#!/usr/bin/env nextflow

/*
 * Standalone contamination-check pipeline.
 * Runs CONTAMINATION_CHECK subworkflow on a directory of *.contigs.fa.gz files.
 *
 * Usage:
 *   nextflow run contamination_check_standalone.nf \
 *       -c contamination_runs/nextflow.config \
 *       --contigs_dir /path/to/spades/ \
 *       [--fastp_dir  /path/to/fastp/] \
 *       [--glue_dir   /path/to/hcvglue/] \
 *       --outdir      contamination_runs/<run_id>/results \
 *       -work-dir     contamination_runs/<run_id>/work
 */

nextflow.enable.dsl = 2

include { CONTAMINATION_CHECK } from './subworkflows/local/contamination_check/main'

params.contigs_dir = null
params.fastp_dir   = null  // optional: directory of *.fastp.json files
params.glue_dir    = null  // optional: directory of *.major.major.*.nodup.json files
params.outdir      = 'results'

workflow {
    if (!params.contigs_dir) {
        error "Please supply --contigs_dir pointing to a folder of *.contigs.fa.gz files"
    }

    ch_contigs = Channel
        .fromPath("${params.contigs_dir}/*.contigs.fa.gz")
        .map { fa ->
            def sample_id = fa.name.replaceAll(/\.contigs\.fa\.gz$/, '')
            [ [id: sample_id], fa ]
        }

    // Collect fastp JSON files if a directory was provided
    ch_fastp = params.fastp_dir
        ? Channel.fromPath("${params.fastp_dir}/*.fastp.json")
        : Channel.empty()

    // Collect GLUE major JSON files if a directory was provided
    ch_glue = params.glue_dir
        ? Channel.fromPath("${params.glue_dir}/*_major.major.nodup.json")
        : Channel.empty()

    CONTAMINATION_CHECK(ch_contigs, ch_fastp, ch_glue)
}
