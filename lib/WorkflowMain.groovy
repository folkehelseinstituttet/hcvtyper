//
// This file holds several functions specific to the main.nf workflow in the niph/viralseq pipeline
//

import nextflow.Nextflow

class WorkflowMain {

    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            // TODO nf-core: Add Zenodo DOI for pipeline after first release
            //"* The pipeline\n" +
            //"  https://doi.org/10.5281/zenodo.XXXXXXX\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }


    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {

        // Print workflow version and exit on --version
        if (params.version) {
            String workflow_version = NfcoreTemplate.version(workflow)
            log.info "${workflow.manifest.name} ${workflow_version}"
            System.exit(0)
        }

        // Check that a -profile or Nextflow config has been provided to run the pipeline
        NfcoreTemplate.checkConfigProvided(workflow, log)

        // Check that conda channels are set-up correctly
        if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
            Utils.checkCondaChannels(log)
        }

        // Check AWS batch settings
        NfcoreTemplate.awsBatch(workflow, params)

        // Check input has been provided
        if (!params.input) {
            Nextflow.error("Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'")
        }
    }

        //
    // Get attribute from genome config file e.g. fasta
    //
    public static String getGenomeAttribute(params, attribute, log, primer_set='', primer_set_version=0) {
        def val = ''
        def support_link =  " The default genome config used by the pipeline can be found here:\n" +
                            "   - https://github.com/nf-core/configs/blob/master/conf/pipeline/viralrecon/genomes.config\n\n" +
                            " If you would still like to blame us please come and find us on nf-core Slack:\n" +
                            "   - https://nf-co.re/viralrecon#contributions-and-support\n" +
                            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
            def genome_map = params.genomes[ params.genome ]
            if (primer_set) {
                if (genome_map.containsKey('primer_sets')) {
                    genome_map = genome_map[ 'primer_sets' ]
                    if (genome_map.containsKey(primer_set)) {
                        genome_map = genome_map[ primer_set ]
                        primer_set_version = primer_set_version.toString()
                        if (genome_map.containsKey(primer_set_version)) {
                            genome_map = genome_map[ primer_set_version ]
                        } else {
                            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                                " --primer_set_version '${primer_set_version}' not found!\n\n" +
                                " Currently, the available primer set version keys are: ${genome_map.keySet().join(", ")}\n\n" +
                                " Please check:\n" +
                                "   - The value provided to --primer_set_version (currently '${primer_set_version}')\n" +
                                "   - The value provided to --primer_set (currently '${primer_set}')\n" +
                                "   - The value provided to --genome (currently '${params.genome}')\n" +
                                "   - Any custom config files provided to the pipeline.\n\n" + support_link
                            System.exit(1)
                        }
                    } else {
                        log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                            " --primer_set '${primer_set}' not found!\n\n" +
                            " Currently, the available primer set keys are: ${genome_map.keySet().join(", ")}\n\n" +
                            " Please check:\n" +
                            "   - The value provided to --primer_set (currently '${primer_set}')\n" +
                            "   - The value provided to --genome (currently '${params.genome}')\n" +
                            "   - Any custom config files provided to the pipeline.\n\n" + support_link
                        System.exit(1)
                    }
                } else {
                    log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                        " Genome '${params.genome}' does not contain any primer sets!\n\n" +
                        " Please check:\n" +
                        "   - The value provided to --genome (currently '${params.genome}')\n" +
                        "   - Any custom config files provided to the pipeline.\n\n" + support_link
                    System.exit(1)
                }
            }
            if (genome_map.containsKey(attribute)) {
                val = genome_map[ attribute ]
            } else if (params.genomes[ params.genome ].containsKey(attribute)) {
                val = params.genomes[ params.genome ][ attribute ]
            }
        }
        return val
    }
}
