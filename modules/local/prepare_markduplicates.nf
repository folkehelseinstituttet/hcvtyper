process PREPARE_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.18--h50ea8bc_1' :
        'community.wave.seqera.io/library/samtools_biopython:292d1d97907ae775' }"

    input:
    tuple val(meta), path(bam), path(contigs)

    output:
    tuple val(meta), path("*NODE*.bam"), path("*.fasta"), emit: bam_contig
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Extract the contig name from the header
    contig=\$(samtools head ${bam} | awk 'NR==2 {print \$2}' | sed 's/^SN://')
    echo \$contig > ${prefix}.contig

    # Then extract the corresponding contig from the contigs file
    prepare_markduplicates.py $contigs \$contig

    mv ${prefix}.bam ${prefix}.\$contig.bam

    # Rename the bam file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        BioPython: \$(python -c "import Bio; print(Bio.__version__)")
        Python: \$(python --version 2>&1 | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}


