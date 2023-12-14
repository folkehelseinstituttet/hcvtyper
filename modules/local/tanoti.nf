process TANOTI_ALIGN {
    tag "$meta.id"
    label "process_low" // process_high

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'docker.io/jonbra/viral_haplo:1.3' }"

    input:
    tuple val(meta), path(reads)
    path (references)
    val   save_unaligned
    val   sort_bam
    val reference
    val prefix2
    val stringency

    output:
    tuple val(meta), path("*markdup.bam")      , emit: aligned
    tuple val(meta), path("*withdup.bam")      , emit: aligned_withdup
    //tuple val(meta), path("*.log")            , emit: log
    tuple val(meta), path("*fastq.gz")        , emit: fastq, optional:true
    tuple val(meta), path("*markdup.idxstats")       , emit: idxstats
    tuple val(meta), path("*withdup.stats")          , emit: stats_withdup
    tuple val(meta), path("*markdup.stats")          , emit: stats_markdup
    tuple val(meta), path("*markdup.coverage.txt.gz"), emit: depth
    path  "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def args2 = task.ext.args2 ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    def unaligned = ""
    def reads_args = ""
    if (meta.single_end) {
        unaligned = save_unaligned ? "--un-gz ${prefix}.unmapped.fastq.gz" : ""
        reads_args = "-U ${reads}"
    } else {
        unaligned = save_unaligned ? "--un-conc-gz ${prefix}.unmapped.fastq.gz" : ""
        reads_args = "${reads[0]} ${reads[1]}"
    }

    def samtools_command = sort_bam ? 'sort' : 'view'
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension = (args2 ==~ extension_pattern) ? (args2 =~ extension_pattern)[0][2].toLowerCase() : "bam"

    """
    gzip -cd ${reads[0]} > ${prefix}_R1.fastq
    gzip -cd ${reads[1]} > ${prefix}_R2.fastq

    tanoti \\
        -r $references \\
        -i tmp_R1.fastq tmp_R2.fastq \\
        -o ${prefix}.${reference}.${prefix2}.tmp.sam \\
        -p 1 -u 0 -m ${stringency}

    samtools $samtools_command $args2 -@ $task.cpus -o ${prefix}.${reference}.${prefix2}.withdup.bam ${prefix}.${reference}.${prefix2}.tmp.sam

    # Create stats file for summary later with duplicates included
    samtools stats ${prefix}.${reference}.${prefix2}.withdup.bam > ${prefix}.${reference}.${prefix2}.withdup.stats

    # Remove duplicate reads
    samtools sort -n ${prefix}.${reference}.${prefix2}.withdup.bam \\
      | samtools fixmate -m - - \
      | samtools sort -O BAM \
      | samtools markdup --no-PG -r - ${prefix}.${reference}.${prefix2}.markdup.bam

    # Creating file with coverage per site
    samtools depth -aa -d 1000000 ${prefix}.${reference}.${prefix2}.markdup.bam | gzip > ${prefix}.${reference}.${prefix2}.markdup.coverage.txt.gz

    # Summarize reads mapped per reference
    samtools idxstats ${prefix}.${reference}.${prefix2}.markdup.bam > ${prefix}.${reference}.${prefix2}.markdup.idxstats

    # Create stats file for summary later
    samtools stats ${prefix}.${reference}.${prefix2}.markdup.bam > ${prefix}.${reference}.${prefix2}.markdup.stats

    if [ -f ${prefix}.unmapped.fastq.1.gz ]; then
        mv ${prefix}.unmapped.fastq.1.gz ${prefix}.unmapped_1.fastq.gz
    fi

    if [ -f ${prefix}.unmapped.fastq.2.gz ]; then
        mv ${prefix}.unmapped.fastq.2.gz ${prefix}.unmapped_2.fastq.gz
    fi

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

    """
    touch ${prefix}.${extension}
    touch ${prefix}.bowtie2.log
    touch ${prefix}.unmapped_1.fastq.gz
    touch ${prefix}.unmapped_2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

}
