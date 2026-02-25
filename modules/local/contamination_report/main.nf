process CONTAMINATION_REPORT {
    label 'process_single'

    // Reuse the BLASTPARSE container — R + tidyverse + jsonlite
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b5/b5d7b977f4b94903794f65bbd852248eef87ce608c5eef34d605afabf514f397/data':
        'community.wave.seqera.io/library/r-seqinr_r-tidyverse:5358395134867368' }"

    input:
    path blast_txt // tabular BLAST output (outfmt 6 with custom columns)

    output:
    path "*.contamination_pairs.tsv" , emit: tsv
    path "*.contamination_heatmap.png", emit: heatmap
    path "*.contamination_mqc.json"  , emit: mqc
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "contamination"
    """
    contamination_report.R \\
        ${blast_txt} \\
        ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
        jsonlite: \$(Rscript -e "library(jsonlite); cat(as.character(packageVersion('jsonlite')))")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "contamination"
    """
    touch ${prefix}.contamination_pairs.tsv
    touch ${prefix}.contamination_heatmap.png
    touch ${prefix}.contamination_mqc.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
        jsonlite: \$(Rscript -e "library(jsonlite); cat(as.character(packageVersion('jsonlite')))")
    END_VERSIONS
    """
}
