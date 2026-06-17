process RESCUE_EVALUATION {
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
    tuple val(meta), path(candidates_csv), path(blastparse_csv), path(support_csv), path(cand_fastas)
    path(references)

    output:
    tuple val(meta), path("*.candidates.csv"), emit: candidates
    tuple val(meta), path("*_cand*.fa"),       emit: candidate_fasta, optional: true
    path "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Stage pass-through per-rank FASTAs into the work dir first (D-07): the
    # Rscript writes ONLY the replaced rescue FASTA, so unchanged ranks must be
    # copied here so the *_cand*.fa emit collects BOTH. `cp -n` never clobbers,
    # and the rescue FASTA is written AFTER (below) so it replaces the same-rank
    # pass-through by matching the _cand{rank}. basename.
    for f in ${cand_fastas}; do
        [ -e "\$f" ] && cp -n "\$f" "./\$(basename \$f)" 2>/dev/null || true
    done

    rescue_evaluation.R \\
        $prefix \\
        $candidates_csv \\
        $blastparse_csv \\
        $support_csv \\
        $references \\
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

    # Stub corrected candidates CSV: the 8-col PARSEFIRSTMAPPING header PLUS the
    # two rescue audit columns (rescued_from, rescue_trigger) + 2 rows so a
    # -stub-run fans out into cand_1/cand_2. Must byte-match the *.candidates.csv
    # emit glob.
    printf "sample,candidate_rank,candidate_ref,candidate_subtype,candidate_genotype,candidate_reads,candidate_cov,confirmation_status,rescued_from,rescue_trigger\n" > ${prefix}.candidates.csv
    printf "${prefix},1,3a_D17763,3a,3,8079,94,pass,NA,NA\n" >> ${prefix}.candidates.csv
    printf "${prefix},2,4k_EU392173,4k,4,40,5,below_threshold,NA,NA\n" >> ${prefix}.candidates.csv

    # Per-rank cand FASTA outputs (one per stub candidates.csv row). Filenames
    # match the declared emit glob "*_cand*.fa" (underscore before cand).
    : > ${prefix}.stubref_cand1.fa
    : > ${prefix}.stubref_cand2.fa

    # Stable versions file. Plain echo lines (no heredoc) so the output is
    # immune to Groovy script-indent vs bash tab-strip mismatch. Deterministic.
    echo '"${task.process}":' > versions.yml
    echo '  r-base: stub' >> versions.yml
    echo '  tidyverse: stub' >> versions.yml
    echo '  seqinr: stub' >> versions.yml
    """
}
