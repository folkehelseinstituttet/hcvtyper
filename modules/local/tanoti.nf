process TANOTI_ALIGN {
    tag "$meta.id"
    label "process_low"
    label "process_long"

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'docker.io/jonbra/viral_haplo:1.3' }"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(references)
    val   sort_bam
    val   stringency

    output:
    tuple val(meta), path("*.{bam}"), emit: aligned
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def args2 = task.ext.args2 ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    def samtools_command = sort_bam ? 'sort' : 'view'
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension = (args2 ==~ extension_pattern) ? (args2 =~ extension_pattern)[0][2].toLowerCase() : "bam"

    """
    gzip -cd ${reads[0]} > ${prefix}_R1.fastq
    gzip -cd ${reads[1]} > ${prefix}_R2.fastq

    tanoti \\
        -r $references \\
        -i ${prefix}_R1.fastq ${prefix}_R2.fastq \\
        -o ${prefix}.sam \\
        -p 1 -u 0 -m ${stringency}

    samtools $samtools_command $args2 -@ $task.cpus -o ${prefix}.bam ${prefix}.sam


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tanoti: Latest commit on master branch https://github.com/vbsreenu/Tanoti/tree/master
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args2 = task.ext.args2 ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension = (args2 ==~ extension_pattern) ? (args2 =~ extension_pattern)[0][2].toLowerCase() : "bam"
    def create_unmapped = ""

    """
    touch ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tanoti: Latest commit on master branch https://github.com/vbsreenu/Tanoti/tree/master
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}
