process CONSENSUS_DISTANCE {
    tag "$meta.id"
    label 'process_low'

    // Uses the same container as SUMMARIZE (has r-seqinr)
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2a/2a1764abd77b9638883a202b96952a48f46cb0ee6c4f65874b836b9455a674d1/data':
        'community.wave.seqera.io/library/r-gridextra_r-png_r-seqinr_r-tidyverse:3536dd50a17de0ab' }"

    input:
    tuple val(meta), path(consensus)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*.consensus_distance.tsv"), emit: tsv
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    consensus_distance.R \\
        $consensus \\
        $reference \\
        ${prefix}.consensus_distance.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        seqinr: \$(Rscript -e "library(seqinr); cat(as.character(packageVersion('seqinr')))")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat > ${prefix}.consensus_distance.tsv << 'EOF'
sample\treference\tsimilarity_pct\tn_differences\talignment_length\tconsensus_length
${prefix}\tref\t99.5\t10\t2000\t1800
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        seqinr: \$(Rscript -e "library(seqinr); cat(as.character(packageVersion('seqinr')))")
    END_VERSIONS
    """
}
