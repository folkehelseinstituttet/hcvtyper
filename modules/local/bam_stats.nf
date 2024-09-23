process BAM_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.18--h50ea8bc_1' :
        'biocontainers/samtools:1.18--h50ea8bc_1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.stats"), emit: stats
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools head ${bam} > ${prefix}.head

    # Extract the contig name from the header
    contig1=\$(awk 'NR==2 {print \$2}' ${prefix}.head | sed 's/^SN://')
    tmp=\$(cat ${prefix}.head | head -3 | tail -1)
    # Cut on tab delimiter and extract the 9th field, then remove the first 15 characters
    tmp2=\$(echo \$tmp | cut -f9 -d' ' | cut -c 11-)
    contig2=\$(echo \$tmp2 | cut -c 6-)
    #contig2=\$(awk 'NR==3 {match(\$0, /.*bowtie2\\/([^ ]+)/, arr); print arr[1]}' ${prefix}.head)

    echo \$contig1 >> ${prefix}.tmp
    echo \$contig2 >> ${prefix}.tmp

    # Check if contig1 and contig2 are identical
    # Not using at the moment
    if [ "\$contig1" = "\$contig1" ]; then
        samtools \\
            stats \\
            --threads ${task.cpus} \\
            ${bam} \\
            > ${prefix}.\$contig1.stats
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
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


