process COMBINE_GFF_FASTA {
    tag "$meta.id"
    label 'process_single'

    conda "YOUR-TOOL-HERE"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'docker.io/pegi3s/biopython:1.78' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*collected_gffs.fasta"), emit: collected_gffs, optional: true
    path  "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat $fasta > ${meta.id}_collected_gffs.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        BioPython: \$(python -c "import Bio; print(Bio.__version__)")
        Python: \$(python --version 2>&1 | cut -d' ' -f2)
    END_VERSIONS
    """

}
