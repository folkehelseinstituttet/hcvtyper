process INSTRUMENT_ID {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
       'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
       'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*sequencerID.tsv"), emit: id

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Get the first line of the fastq file
    sqid=\$(gzip -cd ${reads[0]} | awk 'FNR <= 1')
    echo \${sqid} > ${prefix}_sequencerID.tsv
    """
}
