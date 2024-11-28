process VIGOR {

    label 'process_medium'
    errorStrategy 'terminate'
    containerOptions '-u $(id -u):$(id -g)'

    conda "YOUR-TOOL-HERE"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'docker.io/jonbra/vigor4:1.0' }"

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("*vigor4_out.gff3"), emit: gff3
    tuple val(meta), path("*contigs.fasta")  , emit: contigs
    path  "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    gzip -cd $contigs > ${prefix}_contigs.fasta
    /home/vigor4/vigor4/bin/vigor4 \\
        -i ${prefix}_contigs.fasta \\
        -o ${meta.id}_vigor4_out \\
        -d rtva

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vigor4: \$( /home/vigor4/vigor4/bin/vigor4 --version )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
