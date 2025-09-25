process SUMMARIZE {

    label 'process_medium'
    errorStrategy 'terminate'

    // Environment with R tidyverse and seqinr packages from the conda-forge channel. Created using seqera containers.
    // URLs:
    // Docker image: https://wave.seqera.io/view/builds/bd-3536dd50a17de0ab_1?_gl=1*16bm7ov*_gcl_au*MTkxMjgxNTMwMi4xNzUzNzczOTQz
    // Singularity image: https://wave.seqera.io/view/builds/bd-88101835c4571845_1?_gl=1*5trzpp*_gcl_au*MTkxMjgxNTMwMi4xNzUzNzczOTQz
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2a/2a1764abd77b9638883a202b96952a48f46cb0ee6c4f65874b836b9455a674d1/data':
        'community.wave.seqera.io/library/r-gridextra_r-png_r-seqinr_r-tidyverse:3536dd50a17de0ab' }"

    input:
    path samplesheet
    val stringency_1
    val stringency_2
    path 'trimmed/'
    path 'kraken_classified/'
    path 'parsefirst_mapping/'
    path 'stats_withdup/'
    path 'stats_markdup/'
    path 'depth/'
    path 'blast/'
    path 'glue/'
    path 'id/'
    path 'variation/'

    output:
    path 'Summary.csv'      , emit: summary
    path '*mqc.csv'         , emit: mqc
    path '*LW_import.tsv'   , emit: lw
    path '*png'             , emit: png
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    summarize.R \\
        $samplesheet \\
        $stringency_1 \\
        $stringency_2 \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
        seqinr: \$(Rscript -e "library(seqinr); cat(as.character(packageVersion('seqinr')))")
        gridExtra: \$(Rscript -e "library(gridExtra); cat(as.character(packageVersion('gridExtra')))")
        png: \$(Rscript -e "library(png); cat(as.character(packageVersion('png')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

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
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
        seqinr: \$(Rscript -e "library(seqinr); cat(as.character(packageVersion('seqinr')))")
        gridextra: \$(Rscript -e "library(gridextra); cat(as.character(packageVersion('gridextra')))")
        png: \$(Rscript -e "library(png); cat(as.character(packageVersion('png')))")
    END_VERSIONS
    """
}
