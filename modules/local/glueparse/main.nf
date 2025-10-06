process GLUEPARSE {

    label 'process_single'

    // Environment with R tidyverse and jsonlite packages from the conda-forge channel. Created using seqera containers.
    // URLs:
    // Docker image: https://wave.seqera.io/view/builds/bd-fcb57a72bf127ac8_1?_gl=1*ipusgc*_gcl_au*MTkxMjgxNTMwMi4xNzUzNzczOTQz
    // Singularity image: https://wave.seqera.io/view/builds/bd-739b63e9fb54f431_1?_gl=1*420hyj*_gcl_au*MTkxMjgxNTMwMi4xNzUzNzczOTQz
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a2/a29ad4a97f1b927d4c34550beb44729a7ce4c95be62ede9a96153939c0878b79/data':
        'community.wave.seqera.io/library/r-jsonlite_r-tidyverse:2.0.0--fcb57a72bf127ac8' }"

    input:
    path '*'

    output:
    path("*.tsv")      , emit: GLUE_summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    GLUE_json_parser.R \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
        jsonlite: \$(Rscript -e "library(jsonlite); cat(as.character(packageVersion('jsonlite')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    # Safe echo under -u (interpolated by Groovy)
    echo "${args}"

    # Deterministic minimal output matching emit pattern
    # Write the expected TSV header
    printf "Sample\tReference\tMajor_minor\tGLUE_genotype\tGLUE_subtype\tglecaprevir\tglecaprevir_mut\tglecaprevir_mut_short\tgrazoprevir\tgrazoprevir_mut\tgrazoprevir_mut_short\tparitaprevir\tparitaprevir_mut\tparitaprevir_mut_short\tvoxilaprevir\tvoxilaprevir_mut\tvoxilaprevir_mut_short\tNS34A\tNS34A_short\tdaclatasvir\tdaclatasvir_mut\tdaclatasvir_mut_short\telbasvir\telbasvir_mut\telbasvir_mut_short\tledipasvir\tledipasvir_mut\tledipasvir_mut_short\tombitasvir\tombitasvir_mut\tombitasvir_mut_short\tpibrentasvir\tpibrentasvir_mut\tpibrentasvir_mut_short\tvelpatasvir\tvelpatasvir_mut\tvelpatasvir_mut_short\tNS5A\tNS5A_short\tdasabuvir\tdasabuvir_mut\tdasabuvir_mut_short\tsofosbuvir\tsofosbuvir_mut\tsofosbuvir_mut_short\tNS5B\tNS5B_short\tHCV project version\tGLUE engine version\tPHE drug resistance extension version\n" > GLUE_summary.tsv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
        jsonlite: \$(Rscript -e "library(jsonlite); cat(as.character(packageVersion('jsonlite')))")
    END_VERSIONS
    """
}
