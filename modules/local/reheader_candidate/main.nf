process REHEADER_CANDIDATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    // Per-candidate region-split BAM (reads already restricted to candidate_ref) + that
    // candidate's reference FASTA.
    //
    // The `samtools view <region>` split keeps EVERY @SQ line in the BAM header, so a
    // per-candidate BAM still advertises all candidate references. Downstream per-candidate
    // tools then misbehave: `samtools depth -aa` emits all-zero rows for the other references
    // and summarize.R's `ref_length <- nrow(cov)` / cv_evenness span every reference's
    // positions (Pitfall 2), and bam_coverage.R sees >1 distinct reference and aborts.
    //
    // Rebuild a SINGLE-reference BAM: re-encode the read records (whose RNAME in SAM text is
    // the reference NAME, not the header index) against a header carrying only the candidate
    // reference's @SQ. The SAM round-trip renumbers RNAME indices correctly — a plain
    // `samtools reheader` would corrupt them whenever the kept reference is not first in the
    // original @SQ order (e.g. the 2nd candidate). Result: a true single-reference BAM
    // identical in spirit to the legacy per-candidate independent mapping (JMAP-02).
    tuple val(meta), path(bam), path(fasta)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.reheadered"
    if ("${bam}" == "${prefix}.bam") error "Input and output names are the same, use task.ext.prefix to disambiguate!"
    """
    # Single @SQ line (name + length) from the candidate reference FASTA index.
    samtools faidx ${fasta}
    {
        # Preserve the original @HD (sort order), then a single @SQ, then non-@SQ/@HD header
        # lines (e.g. @PG, @RG), then the read records re-encoded against this header.
        # grep -E instead of grep -P: the patterns use no Perl-specific features and -E is
        # POSIX standard, supported by BusyBox grep (which rejects -P with a fatal error).
        samtools view -H ${bam} | grep -E '^@HD' || true
        awk 'BEGIN{FS="\\t"; OFS="\\t"}{print "@SQ", "SN:"\$1, "LN:"\$2}' ${fasta}.fai
        samtools view -H ${bam} | grep -vE '^@HD|^@SQ' || true
        samtools view ${bam}
    } | samtools view -b -o ${prefix}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.reheadered"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
