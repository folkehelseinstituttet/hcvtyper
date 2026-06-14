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
    path(genotype_utils)

    output:
    // Phase 9 (COMPAT-02 / D-01 / D-04): N-FASTA candidate emit — one entry per
    // ranked candidate (`*_cand*.fa` written by summarize_mapping_to_all_references.R),
    // replacing the fixed two-slot major_mapping/minor_mapping pair. The Plan-02
    // workflow fan-out joins this by full meta and iterates per rank.
    // The `*.parsefirstmapping.csv` glob is pinned to the specific filename so it
    // never also captures the long-format `*.candidates.csv` written alongside it
    // (T-06-05 glob collision).
    tuple val(meta), path("*.parsefirstmapping.csv"), path("*_cand*.fa"), emit: candidate_fasta, optional: true
    tuple val(meta), path("*.parsefirstmapping.csv"),                    emit: csv,             optional: true
    // Phase 6 (REFSEL-01): long-format candidate table — one row per neutrally-
    // ranked candidate. The routable channel the Plan-03 workflow fan-out consumes.
    tuple val(meta), path("*.candidates.csv"),                           emit: candidates,    optional: true
    path "versions.yml",                                                 emit: versions

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
        ${params.minRead} \\
        ${params.minCov} \\
        ${params.n_candidates} \\
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

    # Stub CSV must match the real script's output contract: filename
    # (${prefix}.parsefirstmapping.csv) and the full 10-column header so a
    # -stub-run of the workflow lifts `sample`/`minor_call` into the meta map
    # (the downstream `id == sample` assert and `minor_call == 'yes'` filter).
    printf "sample,total_mapped_reads,major_ref,major_reads,major_cov,minor_ref,minor_reads,minor_cov,minor_call,gate_flag\n" > ${prefix}.parsefirstmapping.csv
    printf "${prefix},8119,3a_D17763,8079,94,4k_EU392173,40,5,no,ok\n" >> ${prefix}.parsefirstmapping.csv

    # Stub long-format candidates CSV (REFSEL-01 contract). Full header + 2 rows so
    # a -stub-run of the workflow fans out into 2 candidates (cand_1/cand_2), matching
    # the default n_candidates=2 two-slot topology consumed by the Plan-03 fan-out.
    printf "sample,candidate_rank,candidate_ref,candidate_subtype,candidate_genotype,candidate_reads,candidate_cov,confirmation_status\n" > ${prefix}.candidates.csv
    printf "${prefix},1,3a_D17763,3a,3,8079,94,pass\n" >> ${prefix}.candidates.csv
    printf "${prefix},2,4k_EU392173,4k,4,40,5,below_threshold\n" >> ${prefix}.candidates.csv

    # Optional per-rank cand FASTA outputs (one per stub candidates.csv row).
    : > ${prefix}.cand1.fa
    : > ${prefix}.cand2.fa

    # Stable versions file. Plain echo lines (no heredoc) so the output is
    # immune to Groovy script-indent stripping vs bash <<- tab-stripping — the
    # mismatch that left a literal END_VERSIONS and stray indentation in the
    # stub versions.yml. Deterministic, byte-stable across runs.
    echo '"${task.process}":' > versions.yml
    echo '  r-base: stub' >> versions.yml
    echo '  tidyverse: stub' >> versions.yml
    echo '  seqinr: stub' >> versions.yml
    """
}
