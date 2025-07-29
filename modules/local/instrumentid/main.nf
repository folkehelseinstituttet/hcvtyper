process INSTRUMENTID {
    tag "$meta.id"
    label 'process_single'

    // Environment with the seqkit software. Created using seqera containers.
    // URLs:
    // Docker image: https://wave.seqera.io/view/builds/bd-03b4774218b4b7ef_1?_gl=1*1mobsod*_gcl_au*MTkxMjgxNTMwMi4xNzUzNzczOTQz
    // Singularity image: https://wave.seqera.io/view/builds/bd-9a5d37887d7c4e09_1?_gl=1*1mobsod*_gcl_au*MTkxMjgxNTMwMi4xNzUzNzczOTQz
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/90/902e5f138e9145c41c3e7848a9d8b7b7f8a21335ac76ac811b96cabcb3d277ad/data':
        'community.wave.seqera.io/library/seqkit:2.10.0--03b4774218b4b7ef' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*sequencerID.tsv"), emit: id
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Get the first line of the fastq file
    sqid=\$(seqkit head -n 1 ${reads[0]} | awk 'FNR <= 1')
    echo \${sqid} > ${prefix}.sequencerID.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    // TODO nf-core: If the module doesn't use arguments ($args), you SHOULD remove:
    //               - The definition of args `def args = task.ext.args ?: ''` above.
    //               - The use of the variable in the script `echo $args ` below.
    """
    echo $args

    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        instrumentid: \$(instrumentid --version)
    END_VERSIONS
    """
}
