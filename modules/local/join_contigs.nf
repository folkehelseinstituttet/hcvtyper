process JOIN_CONTIGS {
    tag "$meta.id"
    label 'process_single'

    conda "YOUR-TOOL-HERE"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'docker.io/pegi3s/biopython:1.78' }"

    input:
    tuple val(id), val(meta),  path(fasta)

    output:
    tuple val(id), path("*.fasta"), emit: contig, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat *.fasta > ${id.sample}_contigs.fasta
    """

}
