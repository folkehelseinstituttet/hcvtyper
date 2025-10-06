process BLASTPARSE {
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
    tuple val(meta), path(blast_out), path(contigs)
    path(references)
    val(agens)

    output:
    tuple val(meta), path("*contigs.fa")    , emit: contigs    , optional: true
    tuple val(meta), path('*blast_out.csv') , emit: blast_res
    tuple val(meta), path("*major.fa")      , emit: major_fasta, optional: true
    tuple val(meta), path("*minor.fa")      , emit: minor_fasta, optional: true
    tuple val(meta), path("*blastparse.csv"), emit: csv
    tuple val(meta), path("*.png")          , emit: png
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def agens  = task.ext.agens ?: "${agens}"

    """
    blast_parse.R \\
        $prefix \\
        $blast_out \\
        $contigs \\
        $references \\
        $agens \\
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
    echo "${args}"

    # Deterministic stub outputs matching declared outputs
    printf "qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore\n" > ${prefix}.blast_out.csv
    printf "id,header\n" > ${prefix}.blastparse.csv
    : > ${prefix}.png

    # Optional outputs (safe to keep empty)
    : > ${prefix}.contigs.fa
    : > ${prefix}.major.fa
    : > ${prefix}.minor.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
      tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
      seqinr: \$(Rscript -e "library(seqinr); cat(as.character(packageVersion('seqinr')))")
    END_VERSIONS
    """
}
