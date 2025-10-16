process SUMMARIZE {

    label 'process_medium'
    errorStrategy 'terminate'

    // Environment with R tidyverse and seqinr packages from the conda-forge channel. Created using seqera containers.
    // URLs:
    // Docker image: https://wave.seqera.io/view/builds/bd-3536dd50a17de0ab_1?_gl=1*16bm7ov*_gcl_au*MTkxMjgxNTMwMi4xNzUzNzczOTQz
    // Singularity image: https://wave.seqera.io/view/builds/bd-88101835c4571845_1?_gl=1*5trzpp*_gcl_au*MTkxMjgxNTMwMi4xNzUzNzczOTQz
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2a/2a1764abd77b9638883a202b96952a48f46cb0ee6c4f65874b836b9455a674d1/data':
        'community.wave.seqera.io/library/r-gridextra_r-png_r-seqinr_r-tidyverse:3536dd50a17de0ab' }"

    // Set environment variables for summarize.R
    env "PIPELINE_NAME",      "${workflow.manifest.name}"
    env "PIPELINE_VERSION",   "${workflow.manifest.version}"                 // Primary version source
    env "PIPELINE_REVISION",  "${workflow.revision ?: ''}"                  // Secondary (branch/tag info)
    env "PIPELINE_COMMIT",    "${workflow.commitId ?: ''}"                  // Tertiary (commit info)

    input:
    val version
    val name
    path samplesheet
    val stringency_1
    val stringency_2
    path 'trimmed/'
    path 'kraken_classified/'
    path 'parsefirst_mapping/'
    path 'stats_withdup/'
    path 'stats_markdup/'
    path 'depth/'
    path 'blast/'
    path 'glue/'
    path 'id/'
    path 'variation/'

    output:
    path 'Summary.csv'      , emit: summary
    path '*mqc.csv'         , emit: mqc
    path '*png'             , emit: png
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    summarize.R \\
        $samplesheet \\
        $stringency_1 \\
        $stringency_2 \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
        seqinr: \$(Rscript -e "library(seqinr); cat(as.character(packageVersion('seqinr')))")
        gridExtra: \$(Rscript -e "library(gridExtra); cat(as.character(packageVersion('gridExtra')))")
        png: \$(Rscript -e "library(png); cat(as.character(packageVersion('png')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    # Safe echo (interpolated by Groovy)
    echo "${args}"

    # Create realistic Summary.csv with proper header and sample data
    cat > Summary.csv << 'EOF'
sampleName,total_raw_reads,total_trimmed_reads,total_classified_reads,total_mapped_reads,fraction_mapped_reads_vs_median,Major_genotype_mapping,Major_reference,Minor_genotype_mapping,Minor_reference,major_typable,minor_typable,Reads_withdup_mapped_major,Reads_nodup_mapped_major,Percent_reads_mapped_of_trimmed_with_dups_major,Major_cov_breadth_min_5,Major_cov_breadth_min_10,percent_mapped_reads_major_firstmapping,Reads_withdup_mapped_minor,Reads_nodup_mapped_minor,Percent_reads_mapped_of_trimmed_with_dups_minor,Minor_cov_breadth_min_5,Minor_cov_breadth_min_10,percent_mapped_reads_minor_firstmapping,sequencer_id,Reads_nodup_mapped_first_mapping,Major_cov_breadth_min_1,Minor_cov_breadth_min_1,Major_avg_depth,Minor_avg_depth,Reference,GLUE_genotype,GLUE_subtype,glecaprevir,glecaprevir_mut,glecaprevir_mut_short,grazoprevir,grazoprevir_mut,grazoprevir_mut_short,paritaprevir,paritaprevir_mut,paritaprevir_mut_short,voxilaprevir,voxilaprevir_mut,voxilaprevir_mut_short,NS34A,NS34A_short,daclatasvir,daclatasvir_mut,daclatasvir_mut_short,elbasvir,elbasvir_mut,elbasvir_mut_short,ledipasvir,ledipasvir_mut,ledipasvir_mut_short,ombitasvir,ombitasvir_mut,ombitasvir_mut_short,pibrentasvir,pibrentasvir_mut,pibrentasvir_mut_short,velpatasvir,velpatasvir_mut,velpatasvir_mut_short,NS5A,NS5A_short,dasabuvir,dasabuvir_mut,dasabuvir_mut_short,sofosbuvir,sofosbuvir_mut,sofosbuvir_mut_short,NS5B,NS5B_short,HCV project version,GLUE engine version,PHE drug resistance extension version,script_name_stringency
Test_1,40000,34446,31482,31478,1.2804002521914215,NA,NA,NA,NA,NO,NO,NA,NA,NA,NA,NA,95.59,NA,NA,NA,NA,NA,3.55,@SRR24174266.1 1/1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,HCVTyper (version unknown)
EOF

    # Create identical summary_mqc.csv (same data, MultiQC format)
    cat > summary_mqc.csv << 'EOF'
sampleName,total_raw_reads,total_trimmed_reads,total_classified_reads,total_mapped_reads,fraction_mapped_reads_vs_median,Major_genotype_mapping,Major_reference,Minor_genotype_mapping,Minor_reference,major_typable,minor_typable,Reads_withdup_mapped_major,Reads_nodup_mapped_major,Percent_reads_mapped_of_trimmed_with_dups_major,Major_cov_breadth_min_5,Major_cov_breadth_min_10,percent_mapped_reads_major_firstmapping,Reads_withdup_mapped_minor,Reads_nodup_mapped_minor,Percent_reads_mapped_of_trimmed_with_dups_minor,Minor_cov_breadth_min_5,Minor_cov_breadth_min_10,percent_mapped_reads_minor_firstmapping,sequencer_id,Reads_nodup_mapped_first_mapping,Major_cov_breadth_min_1,Minor_cov_breadth_min_1,Major_avg_depth,Minor_avg_depth,Reference,GLUE_genotype,GLUE_subtype,glecaprevir,glecaprevir_mut,glecaprevir_mut_short,grazoprevir,grazoprevir_mut,grazoprevir_mut_short,paritaprevir,paritaprevir_mut,paritaprevir_mut_short,voxilaprevir,voxilaprevir_mut,voxilaprevir_mut_short,NS34A,NS34A_short,daclatasvir,daclatasvir_mut,daclatasvir_mut_short,elbasvir,elbasvir_mut,elbasvir_mut_short,ledipasvir,ledipasvir_mut,ledipasvir_mut_short,ombitasvir,ombitasvir_mut,ombitasvir_mut_short,pibrentasvir,pibrentasvir_mut,pibrentasvir_mut_short,velpatasvir,velpatasvir_mut,velpatasvir_mut_short,NS5A,NS5A_short,dasabuvir,dasabuvir_mut,dasabuvir_mut_short,sofosbuvir,sofosbuvir_mut,sofosbuvir_mut_short,NS5B,NS5B_short,HCV project version,GLUE engine version,PHE drug resistance extension version,script_name_stringency
Test_1,40000,34446,31482,31478,1.2804002521914215,NA,NA,NA,NA,NO,NO,NA,NA,NA,NA,NA,95.59,NA,NA,NA,NA,NA,3.55,@SRR24174266.1 1/1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,HCVTyper (version unknown)
EOF

    # Create a minimal PNG plot placeholder
    : > summary_plot.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
        seqinr: \$(Rscript -e "library(seqinr); cat(as.character(packageVersion('seqinr')))")
        gridExtra: \$(Rscript -e "library(gridExtra); cat(as.character(packageVersion('gridExtra')))")
        png: \$(Rscript -e "library(png); cat(as.character(packageVersion('png')))")
    END_VERSIONS
    """
}
