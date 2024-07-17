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

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Untar the different reference fiiles
    tar -xzf $references

    # Loop through gene fastas, extract gene name from fasta file,
    # and merge with correspooding reference file
    for fasta_file in $prefix*highest_cov.fasta
    do
        gene=\$(echo \$fasta_file | cut -d'.' -f2 | cut -d'_' -f1)

        # Concatenate the fasta file with the corresponding gene references
        cat \$fasta_file \$(ls References_\$gene.fasta) > ${prefix}.\${gene}_merged.fasta
    done
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
