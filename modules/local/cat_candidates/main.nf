process CAT_CANDIDATES {
    tag "$meta.id"
    label 'process_single'

    // Uses only cat — any standard unix container works; reuse BLAST container
    // (same convention as CAT_FILTERED_CONTIGS) so no new image is introduced.
    conda "bioconda::blast=2.17.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0c/0c86cbb145786bf5c24ea7fb13448da5f7d5cd124fd4403c1da5bc8fc60c2588/data':
        'community.wave.seqera.io/library/blast:2.17.0--d4fb881691596759' }"

    input:
    tuple val(meta), path(fastas) // PER-SAMPLE list of passing candidate FASTAs (meta-keyed)

    output:
    tuple val(meta), path("*.combined.fa"), emit: fasta
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Per-sample concatenation: cat of one FASTA is identity (single-candidate degenerate
    // case, D-14). Distinct from CAT_FILTERED_CONTIGS which collects ACROSS samples.
    """
    cat ${fastas} > ${prefix}.combined.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version 2>&1 | head -1 | sed 's/^.*(GNU coreutils) //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.combined.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version 2>&1 | head -1 | sed 's/^.*(GNU coreutils) //')
    END_VERSIONS
    """
}
