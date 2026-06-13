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

    input:
    val version
    val name
    path samplesheet
    path 'trimmed/'
    path 'kraken_classified/'
    path 'parsefirst_mapping/'
    path 'stats_withdup/'
    path 'stats_markdup/'
    path 'depth/'
    path 'denovo/'
    path 'glue/'
    path 'id/'
    path 'variation/'
    path 'consensus_distance/'
    path(genotype_utils)
    path(denovo_confirm)
    path(denovo_layer)
    path(assembly_support_join)

    output:
    path 'Summary.csv'      , emit: summary
    path '*mqc.tsv'         , emit: mqc
    path '*png'             , emit: png
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    summarize.R \\
        $samplesheet \\
        $version \\
        $name \\
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
sampleName,total_raw_reads,total_trimmed_reads,total_classified_reads,total_mapped_reads,fraction_mapped_reads_vs_median,Major_genotype_mapping,Major_reference,Minor_genotype_mapping,Minor_reference,major_typable,minor_typable,Reads_withdup_mapped_major,Reads_nodup_mapped_major,Percent_reads_mapped_of_trimmed_with_dups_major,Major_cov_breadth_min_5,Major_cov_breadth_min_10,percent_mapped_reads_major_firstmapping,Reads_withdup_mapped_minor,Reads_nodup_mapped_minor,Percent_reads_mapped_of_trimmed_with_dups_minor,Minor_cov_breadth_min_5,Minor_cov_breadth_min_10,percent_mapped_reads_minor_firstmapping,sequencer_id,Reads_nodup_mapped_first_mapping,Major_cov_breadth_min_1,Minor_cov_breadth_min_1,Major_avg_depth,Minor_avg_depth,Reference,GLUE_genotype,GLUE_subtype,glecaprevir,glecaprevir_mut,glecaprevir_mut_short,grazoprevir,grazoprevir_mut,grazoprevir_mut_short,paritaprevir,paritaprevir_mut,paritaprevir_mut_short,voxilaprevir,voxilaprevir_mut,voxilaprevir_mut_short,NS34A,NS34A_short,daclatasvir,daclatasvir_mut,daclatasvir_mut_short,elbasvir,elbasvir_mut,elbasvir_mut_short,ledipasvir,ledipasvir_mut,ledipasvir_mut_short,ombitasvir,ombitasvir_mut,ombitasvir_mut_short,pibrentasvir,pibrentasvir_mut,pibrentasvir_mut_short,velpatasvir,velpatasvir_mut,velpatasvir_mut_short,NS5A,NS5A_short,dasabuvir,dasabuvir_mut,dasabuvir_mut_short,sofosbuvir,sofosbuvir_mut,sofosbuvir_mut_short,NS5B,NS5B_short,HCV project version,GLUE engine version,PHE drug resistance extension version,script_name_stringency,denovo_major_ref,denovo_major_contig_length,denovo_minor_ref,denovo_minor_contig_length,minor_denovo_status
Test_1,40000,34446,31482,31478,1.2804002521914215,NA,NA,NA,NA,NO,NO,NA,NA,NA,NA,NA,95.59,NA,NA,NA,NA,NA,3.55,@SRR24174266.1 1/1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,HCVTyper (version unknown),NA,NA,NA,NA,NA
EOF

    # Create identical summary_mqc.tsv (same data, MultiQC format — TSV avoids quoting issues)
    printf 'sampleName\ttotal_raw_reads\ttotal_trimmed_reads\ttotal_classified_reads\ttotal_mapped_reads\tfraction_mapped_reads_vs_median\tMajor_genotype_mapping\tMajor_reference\tMinor_genotype_mapping\tMinor_reference\tmajor_typable\tminor_typable\tReads_withdup_mapped_major\tReads_nodup_mapped_major\tPercent_reads_mapped_of_trimmed_with_dups_major\tMajor_cov_breadth_min_5\tMajor_cov_breadth_min_10\tpercent_mapped_reads_major_firstmapping\tReads_withdup_mapped_minor\tReads_nodup_mapped_minor\tPercent_reads_mapped_of_trimmed_with_dups_minor\tMinor_cov_breadth_min_5\tMinor_cov_breadth_min_10\tpercent_mapped_reads_minor_firstmapping\tsequencer_id\tReads_nodup_mapped_first_mapping\tMajor_cov_breadth_min_1\tMinor_cov_breadth_min_1\tMajor_avg_depth\tMinor_avg_depth\tReference\tGLUE_genotype\tGLUE_subtype\tglecaprevir\tglecaprevir_mut\tglecaprevir_mut_short\tgrazoprevir\tgrazoprevir_mut\tgrazoprevir_mut_short\tparitaprevir\tparitaprevir_mut\tparitaprevir_mut_short\tvoxilaprevir\tvoxilaprevir_mut\tvoxilaprevir_mut_short\tNS34A\tNS34A_short\tdaclatasvir\tdaclatasvir_mut\tdaclatasvir_mut_short\telbasvir\telbasvir_mut\telbasvir_mut_short\tledipasvir\tledipasvir_mut\tledipasvir_mut_short\tombitasvir\tombitasvir_mut\tombitasvir_mut_short\tpibrentasvir\tpibrentasvir_mut\tpibrentasvir_mut_short\tvelpatasvir\tvelpatasvir_mut\tvelpatasvir_mut_short\tNS5A\tNS5A_short\tdasabuvir\tdasabuvir_mut\tdasabuvir_mut_short\tsofosbuvir\tsofosbuvir_mut\tsofosbuvir_mut_short\tNS5B\tNS5B_short\tHCV project version\tGLUE engine version\tPHE drug resistance extension version\tscript_name_stringency\tdenovo_major_ref\tdenovo_major_contig_length\tdenovo_minor_ref\tdenovo_minor_contig_length\tminor_denovo_status\n' > summary_mqc.tsv
    printf 'Test_1\t40000\t34446\t31482\t31478\t1.2804002521914215\tNA\tNA\tNA\tNA\tNO\tNO\tNA\tNA\tNA\tNA\tNA\t95.59\tNA\tNA\tNA\tNA\tNA\t3.55\t@SRR24174266.1 1/1\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tHCVTyper (version unknown)\tNA\tNA\tNA\tNA\tNA\n' >> summary_mqc.tsv

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
