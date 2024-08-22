process PREPARE_MAFFT {

    tag "$meta.id"
    label 'process_single'

    conda "YOUR-TOOL-HERE"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'docker.io/pegi3s/biopython:1.78' }"

    input:
    tuple val(meta),  path(fasta)
    tuple val(meta2), path(references)

    output:
    tuple val(meta), path("*_merged.fasta"), emit: fasta, optional: true
    tuple val(meta), env(gene)             , emit: gene , optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Untar the different reference fiiles
    tar -xzf $references

    # Loop through gene fastas, extract gene name from fasta file,
    # and merge with correspooding reference file
    for fasta_file in $prefix*.fasta
    do
        # OLD: gene=\$(echo \$fasta_file | cut -d'.' -f2 | cut -d'_' -f1)
        gene=\$(echo \$fasta_file | cut -d'.' -f2)

        # Concatenate the fasta file with the corresponding gene references
        cat \$fasta_file \$(ls References_\$gene.fasta) > ${prefix}.\${gene}_merged.fasta
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        BioPython: \$(python -c "import Bio; print(Bio.__version__)")
        Python: \$(python --version 2>&1 | cut -d' ' -f2)
    END_VERSIONS
    """
}
