process CONTIG_DUPLICATES {
    tag "$meta.id"
    label 'process_medium'

    // Same container as the nf-core bowtie2 modules (bowtie2 + samtools + pigz).
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b41b403e81883126c3227fc45840015538e8e2212f13abc9ae84e4b98891d51c/data' :
        'community.wave.seqera.io/library/bowtie2_htslib_samtools_pigz:edeb13799090a2a6' }"

    // Per-contig PCR-duplicate rate of the sample's own reads.
    //
    // The classified reads are mapped back to the sample's own SPAdes contigs and
    // duplicates are marked by fragment coordinates (samtools fixmate -m + markdup).
    // Within one library every molecule is amplified and sequenced alike, so a genuine
    // minor strain shares the library's duplicate rate; reads that arrived from another
    // library after its PCR (index hopping, demultiplexing bleed, pooled carry-over) are
    // a sparse draw from a large pool and are almost all unique. summarize.R compares a
    // minor's duplicate rate with the dominant genotype's.
    //
    // Output: <prefix>.contig_dups.tsv — contig, fragments (read-1 primary alignments),
    // duplicates. Only aligned contigs appear.

    input:
    tuple val(meta), path(contigs), path(reads)

    output:
    tuple val(meta), path("*.contig_dups.tsv"), emit: tsv
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def decompress = contigs.getExtension() == "gz" ? "gunzip -c ${contigs}" : "cat ${contigs}"
    """
    ${decompress} > contigs.fa
    bowtie2-build --threads $task.cpus -q contigs.fa contigs

    bowtie2 \\
        -p $task.cpus \\
        --no-unal \\
        -x contigs \\
        -1 ${reads[0]} -2 ${reads[1]} \\
        $args \\
        2> ${prefix}.contig_dups.bowtie2.log \\
        | samtools fixmate -m - - \\
        | samtools sort -@ $task.cpus -o fixmate.bam -

    samtools markdup -@ $task.cpus fixmate.bam markdup.bam

    # One line per contig: read-1 primary alignments (= fragments) and how many are
    # marked duplicate (flag 0x400). 0x904 drops unmapped, secondary and supplementary.
    printf "contig\\tfragments\\tduplicates\\n" > ${prefix}.contig_dups.tsv
    samtools view -F 0x904 -f 0x40 markdup.bam \\
        | awk -F '\\t' -v OFS='\\t' '{ n[\$3]++; if (int(\$2 / 1024) % 2 == 1) d[\$3]++ }
                                  END { for (c in n) print c, n[c], d[c] + 0 }' \\
        | sort -k2,2nr >> ${prefix}.contig_dups.tsv

    rm -f fixmate.bam markdup.bam contigs.fa contigs.*.bt2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "contig\\tfragments\\tduplicates\\n" > ${prefix}.contig_dups.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: 2.5.4
        samtools: 1.21
    END_VERSIONS
    """
}
