process CONTIG_FILTER_RENAME {
    tag "$meta.id"
    label 'process_single'

    // Reuse the BLAST container — has Python 3 and standard unix tools
    conda "bioconda::blast=2.17.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0c/0c86cbb145786bf5c24ea7fb13448da5f7d5cd124fd4403c1da5bc8fc60c2588/data':
        'community.wave.seqera.io/library/blast:2.17.0--d4fb881691596759' }"

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("*.filtered.fa"), emit: filtered
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def min_length = task.ext.args   ?: "1000"
    def sample_id  = meta.id
    def decompress = contigs.getExtension() == "gz" ? "gunzip -c ${contigs}" : "cat ${contigs}"
    """
    ${decompress} | awk -v sample="${sample_id}" -v min_len="${min_length}" '
        /^>/ {
            if (header != "" && length(seq) >= min_len) {
                print ">" sample "|" header
                print seq
            }
            header = substr(\$0, 2)
            seq    = ""
            next
        }
        { seq = seq \$0 }
        END {
            if (header != "" && length(seq) >= min_len) {
                print ">" sample "|" header
                print seq
            }
        }
    ' > ${prefix}.filtered.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filtered.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1)
    END_VERSIONS
    """
}
