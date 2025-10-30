process PARSEFIRSTMAPPING {
    tag "$meta.id"
    label 'process_single'

    // Environment with R tidyverse and seqinr packages from the conda-forge channel. Created using seqera containers.
    // URLs:
    // Docker image: https://wave.seqera.io/view/builds/bd-5358395134867368_1?_gl=1*1vaclhd*_gcl_au*MTkxMjgxNTMwMi4xNzUzNzczOTQz
    // Singularity image: https://wave.seqera.io/view/builds/bd-0225dab2b8112adf_1?_gl=1*111m8r6*_gcl_au*MTkxMjgxNTMwMi4xNzUzNzczOTQz
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b5/b5d7b977f4b94903794f65bbd852248eef87ce608c5eef34d605afabf514f397/data':
        'community.wave.seqera.io/library/r-seqinr_r-tidyverse:5358395134867368' }"

    input:
    tuple val(meta), path(idxstats), path(depth)
    path(references)

    output:
    tuple val(meta), path("*.csv"), path("*major.fa"), emit: major_mapping, optional: true
    tuple val(meta), path("*.csv"), path("*minor.fa"), emit: minor_mapping, optional: true
    tuple val(meta), path("*.csv"),                    emit: csv,           optional: true
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    summarize_mapping_to_all_references.R \\
        ${idxstats} \\
        ${depth} \\
        ${prefix} \\
        ${references} \\
        $args

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
    """
    # Safe echo (interpolated by Groovy)
    echo "${args}"

    # Create outputs matching the declared emit patterns
    printf "reference,mapped_reads,unmapped_reads,depth_mean\n" > ${prefix}.mapping_summary.csv
    printf "3a_D17763,8079,1,94.3\n" >> ${prefix}.mapping_summary.csv
    printf "4k_EU392173,40,0,5.3\n" >> ${prefix}.mapping_summary.csv

    # Optional FASTA outputs (touch to create empty files)
    : > ${prefix}.major.fa
    : > ${prefix}.minor.fa

    # Stable versions file
    cat > versions.yml <<'YAML'
    "${task.process}":
      r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
      tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
      seqinr: \$(Rscript -e "library(seqinr); cat(as.character(packageVersion('seqinr')))")
    END_VERSIONS
    """
}
