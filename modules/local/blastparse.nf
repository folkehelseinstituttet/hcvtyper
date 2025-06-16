process BLASTPARSE {
    tag "$meta.id"
    label 'process_single'

    conda "r=2.4.4 r-tidyverse=1.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '':
        'docker.io/jonbra/tidyverse_seqinr:2.0' }"

    input:
    tuple val(meta), path(blast_out), path(scaffolds)
    path(references)
    val(agens)

    output:
    tuple val(meta), path("*scaffolds.fa")  , emit: scaffolds
    tuple val(meta), path('*blast_out.csv') , emit: blast_res
    tuple val(meta), path("*major.fa")      , emit: major_fasta
    tuple val(meta), path("*minor.fa")      , emit: minor_fasta, optional: true
    tuple val(meta), path("*blastparse.csv"), emit: csv
    tuple val(meta), path("*.png")          , emit: png
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def agens  = task.ext.agens ?: "${agens}"

    """
    blast_parse.R $prefix $blast_out $scaffolds $references $agens

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
      tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
      seqinr: \$(Rscript -e "library(seqinr); cat(as.character(packageVersion('seqinr')))")
    END_VERSIONS
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
