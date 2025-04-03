process SUMMARIZE_HCV {

    label 'process_single'
    errorStrategy 'terminate'

    conda "YOUR-TOOL-HERE"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'community.wave.seqera.io/library/r-gridextra_r-png_r-seqinr_r-tidyverse:3536dd50a17de0ab' }"

    input:
    val stringency_1
    val stringency_2
    path 'cutadapt/'
    path 'kraken_classified/'
    path 'stats_withdup/'
    path 'stats_markdup/'
    path 'depth/'
    path 'blast/'
    path 'glue/'
    path 'id/'
    path 'variation/'

    output:
    path '*long.csv'        , emit: summary
    path '*mqc.csv'         , emit: mqc
    path '*LW_import.tsv'   , emit: lw
    path '*png'             , emit: png
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    summarize.R $stringency_1 $stringency_2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
      tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

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
