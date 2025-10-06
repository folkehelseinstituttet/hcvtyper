process PLOT_BAMVARIATION {
    tag "$meta.id"
    label 'process_single'

    // Environment with R tidyverse and samtools packages from the conda-forge and bioconda channels. Created using seqera containers.
    // URLs:
    // Docker image: https://wave.seqera.io/view/builds/bd-7bb8dfe78e0951cb_1?_gl=1*1vz3710*_gcl_au*MTM5MTA4NDk2NS4xNzUzNjg2MzUxLjY0MTQxNDc2Ni4xNzU2MzA2NTExLjE3NTYzMDY1MjU.
    // Singularity image: https://wave.seqera.io/view/builds/bd-b5576c3fd661445b_1?_gl=1*1s6azpo*_gcl_au*MTM5MTA4NDk2NS4xNzUzNjg2MzUxLjY0MTQxNDc2Ni4xNzU2MzA2NTExLjE3NTYzMDY1MjU.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-rsamtools_r-tidyverse:b5576c3fd661445b':
        'community.wave.seqera.io/library/bioconductor-rsamtools_r-tidyverse:7bb8dfe78e0951cb' }"

    input:
    tuple val(meta), path(bam)

    output:
    path("*.png")      , emit: png, optional: true
    path("*.tsv")      , emit: tsv, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    plot_bam_variation.R "${bam}" "${prefix}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
        Rsamtools: \$(Rscript -e "library(Rsamtools); cat(as.character(packageVersion('Rsamtools')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Print args safely (Groovy interpolates before bash -u)
    echo "${args}"

    # Minimal, deterministic outputs matching emit patterns
    printf "pos\tref\talt\tdepth\talt_count\tfreq\n" > ${prefix}.variation.tsv
    : > ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
        Rsamtools: \$(Rscript -e "library(Rsamtools); cat(as.character(packageVersion('Rsamtools')))")
    END_VERSIONS
    """
}
