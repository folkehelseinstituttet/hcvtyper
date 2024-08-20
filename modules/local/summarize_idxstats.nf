process SUMMARIZE_IDXSTATS {
    label 'process_single'

    conda "YOUR-TOOL-HERE"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'docker.io/jonbra/tidyverse_seqinr:2.0' }"

    input:
    tuple val(meta), path(idxstats)

    output:
    tuple val(meta), path("*.csv"), emit: csv

    script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    """



}
