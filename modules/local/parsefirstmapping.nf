process PARSEFIRSTMAPPING {
    tag "$meta.id"
    label 'process_single'

    conda "YOUR-TOOL-HERE"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'docker.io/jonbra/tidyverse_seqinr:2.0' }"

    input:
    tuple val(meta), path(idxstats), path(depth)
    path(references)

    output:
    tuple val(meta), path("*.csv")    , emit: csv, optional: true
    tuple val(meta), path("*major.fa"), emit: major_fasta, optional: true
    tuple val(meta), path("*minor.fa"), emit: minor_fasta, optional: true
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    summarize_mapping_to_all_references.R ${idxstats} ${depth} ${prefix} ${references}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
      tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
      seqinr: \$(Rscript -e "library(seqinr); cat(as.character(packageVersion('seqinr')))")
    END_VERSIONS
    """
}
