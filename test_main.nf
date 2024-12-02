#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQSCREEN_FASTQSCREEN } from './modules/nf-core/fastqscreen/fastqscreen/main'

params.reads = "$projectDir/data/*_{1,2}.fastq.gz"
params.outdir = 'results'
params.igenomes_base = 's3://ngi-igenomes/igenomes/'
params.genome = 'GRCh37'
params.bowtie2_index

// Create input channel
Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .map { tuple([ id: it[0] ], it[1]) }
    .set { input_reads }

// Download and prepare FastQ Screen database
process prepare_fastq_screen_db {
    conda "wget"

    input:
    path(fasta)

    output:
    path "fastq_screen_db"

    script:
    """
    mkdir -p fastq_screen_db
    cat <<EOF > fastq_screen_db/fastq_screen.conf
    # FastQ Screen configuration file
    BOWTIE2 bowtie2
    DATABASE     Human   /home/user/data/genomics/homo_sapiens/genome/chr21/sequence/genome
    EOF
    """
}

workflow {
    prepare_fastq_screen_db(params.genome)
    FASTQSCREEN_FASTQSCREEN ( input_reads, prepare_fastq_screen_db.out )
}
