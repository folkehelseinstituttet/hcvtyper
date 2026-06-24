process CONSENSUS_DISTANCE {
    tag "$meta.id"
    label 'process_low'

    // Environment with Bioconductor Biostrings and pwalign packages. Created using seqera containers.
    // TODO: submit environment.yml to Wave (https://wave.seqera.io) to obtain the updated container URLs
    //       after adding bioconductor-biostrings and bioconductor-pwalign.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-pwalign:TO_BE_UPDATED':
        'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-pwalign:TO_BE_UPDATED' }"

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
        bioconductor-biostrings: \$(Rscript -e "library(Biostrings); cat(as.character(packageVersion('Biostrings')))")
        bioconductor-pwalign: \$(Rscript -e "library(pwalign); cat(as.character(packageVersion('pwalign')))")
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
        bioconductor-biostrings: \$(Rscript -e "library(Biostrings); cat(as.character(packageVersion('Biostrings')))")
        bioconductor-pwalign: \$(Rscript -e "library(pwalign); cat(as.character(packageVersion('pwalign')))")
    END_VERSIONS
    """
}
