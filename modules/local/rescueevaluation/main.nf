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
    // cand_fastas are staged into a SUBDIR (input_fastas/) via stageAs so their
    // basenames do NOT collide with the *_cand*.fa OUTPUT names. Nextflow excludes
    // input files from output matching BY NAME, so a pass-through copied back to the
    // top level under its own name would otherwise be shadowed by the input and the
    // candidate_fasta emit would collect NOTHING (#10-03 integration bug).
    tuple val(meta), path(candidates_csv), path(blastparse_csv), path(support_csv), path(cand_fastas, stageAs: 'input_fastas/*')
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
    # Stage pass-through per-rank FASTAs as REAL top-level files (D-07): the Rscript
    # writes ONLY the replaced rescue FASTA, so unchanged ranks must be materialised
    # here so the *_cand*.fa emit collects BOTH. The inputs live in input_fastas/
    # (stageAs above), so copying each up to the top level under its own basename
    # produces a GENUINE task output whose name is not shadowed by an input. When a
    # rank IS rescued, the rescue FASTA embeds the NEW ref name in its basename
    # ({prefix}.{rescue_ref}_cand{rank}.fa) and therefore does NOT overwrite the
    # pass-through ({prefix}.{orig_ref}_cand{rank}.fa); rescue_evaluation.R deletes
    # that stale pass-through itself, so only one FASTA per rank survives. (A stale
    # pass-through left behind would duplicate an @SQ line in the combined per-sample
    # reference and crash BOWTIE2_BUILD / samtools sort.)
    for f in ${cand_fastas}; do
        if [ -e "\$f" ]; then
            cp -L "\$f" "./\$(basename \$f)"
        fi
    done

    # Skip-assembly path (D-10): blastparse_csv / support_csv arrive as EMPTY path
    # inputs ([]), so the staged variable expands to an empty string and the
    # positional thresholds would shift into the missing slots ("Usage:" error).
    # Materialise a header-only placeholder for any empty leg so the 10 positional
    # args stay aligned; rescue_evaluation.R's read_csv_guarded() treats a zero-row
    # file as a typed-empty frame -> candidates pass through, rescue columns NA.
    bp='${blastparse_csv}'
    sup='${support_csv}'
    if [ -z "\$bp" ]; then
        bp="${prefix}.EMPTY.blastparse.csv"
        : > "\$bp"
    fi
    if [ -z "\$sup" ]; then
        sup="${prefix}.EMPTY.assembly_support.csv"
        : > "\$sup"
    fi

    rescue_evaluation.R \\
        $prefix \\
        $candidates_csv \\
        "\$bp" \\
        "\$sup" \\
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
    printf "sample,candidate_rank,candidate_ref,candidate_subtype,candidate_genotype,candidate_reads,candidate_cov,confirmation_status,rescued_from,rescue_trigger\n" > ${prefix}.rescued.candidates.csv
    printf "${prefix},1,3a_D17763,3a,3,8079,94,pass,NA,NA\n" >> ${prefix}.rescued.candidates.csv
    printf "${prefix},2,4k_EU392173,4k,4,40,5,below_threshold,NA,NA\n" >> ${prefix}.rescued.candidates.csv

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
