include { CREATE_JPG              } from '../../modules/local/create_jpg'
include { IQTREE                  } from '../../modules/nf-core/iqtree/main'
include { MAFFT                   } from '../../modules/nf-core/mafft/main'
include { MAFFT as MAFFT_PAIRWISE } from '../../modules/nf-core/mafft/main'
include { PARSE_PHYLOGENY         } from '../../modules/local/parse_phylogeny'
include { PARSE_PHYLOGENY_2       } from '../../modules/local/parse_phylogeny_2'
include { PREPARE_MAFFT           } from '../../modules/local/prepare_mafft'

workflow MAFFT_IQTREE_GENOTYPE {

    take:
    ch_vigorparse   // channel: [ val(meta), path(gene_fasta) ]
    ch_references // channel: [ val(meta), path(gene_references) ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Combine the highes covered gene contig with the corresponding reference dataset for MAFFT input
    //
    PREPARE_MAFFT(
        ch_vigorparse,
        ch_references
    )

    // Transpose the PREPARE_MAFFT output to split the concatenated individual gene fasta files into different channels.
    // This is to be able to parallelize the MAFFT and IQTREE steps for each gene.
    ch_temp  = PREPARE_MAFFT.out.fasta.transpose()

    // Prepare a mafft input channel that has the gene name in the meta map. This is for renaming the output files.
    // Add the gene name from the fasta file names to the meta map like this:
    // [[id:, single_end:, gene:], file path]
    ch_mafft = ch_temp.map { item ->
        def meta = item[0] // Original meta map
        def filePath = item[1].toString() // File path as string

        // Use getName to get the filename from the file path
        def fileName = new File(filePath).getName()

        // Extract the gene name from the file path
        def geneName = fileName.split("_")[0].split("\\.")[1]

        // Add the gene name to the meta map
        def updatedMeta = meta.clone() // Clone the original meta map to avoid modifying the original
        updatedMeta.gene = geneName // Add the gene name

        // Emit the updated item
        return [updatedMeta, filePath]
    }

    //
    // MODULE: Align gene sequences with MAFFT
    //
    MAFFT(
        ch_mafft,
        [ [:], [] ],
        [ [:], [] ],
        [ [:], [] ],
        [ [:], [] ],
        [ [:], [] ],
        false
    )
    ch_versions = ch_versions.mix(MAFFT.out.versions.first())

    // Add empty element to MAFFT.out.fas to comply with IQTREE input
    ch_iqtree = MAFFT.out.fas.map {
        meta ->
        return meta + [[]]
    }

    //
    // MODULE: IQTREE
    //
    IQTREE(
        ch_iqtree,
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(IQTREE.out.versions.first())

    // Here we need to join all the output from a single sample before further parsing

    //
    // MODULE: Parse phylogeny to genotype the sample sequence
    //

    // Add the gene name to the meta of ch_vigorparse ([id:, single_end:, gene:])
    // First transpose the channel to separate the different "highest_cov" fasta files
    ch_temp_2 = ch_vigorparse.transpose()
    ch_parse_phylogeny_temp = ch_temp_2.map { item ->
        def meta = item[0] // Original meta map
        def filePath = item[1].toString() // File path as string

        // Use getName to get the filename from the file path
        def fileName = new File(filePath).getName()

        // Extract the gene name from the file path
        def geneName = fileName.split("_")[0].split("\\.")[1]

        // Add the gene name to the meta map
        def updatedMeta = meta.clone() // Clone the original meta map to avoid modifying the original
        updatedMeta.gene = geneName // Add the gene name

        // Emit the updated item
        return [updatedMeta, filePath]
    }

    // Then join the ch_parse_phylogeny with the IQTREE output to have the phylogeny and the highest cov fasta for the same gene
    ch_parse_phylogeny = ch_parse_phylogeny_temp.join(IQTREE.out.phylogeny)
    PARSE_PHYLOGENY(
        ch_parse_phylogeny,
        ch_references
    )
    ch_versions = ch_versions.mix(PARSE_PHYLOGENY.out.versions.first())

    //
    // MODULE: PURPOSE?
    //
    MAFFT_PAIRWISE(
        PARSE_PHYLOGENY.out.percentcalc_fasta,
        [ [:], [] ],
        [ [:], [] ],
        [ [:], [] ],
        [ [:], [] ],
        [ [:], [] ],
        false
    )
    ch_versions = ch_versions.mix(MAFFT_PAIRWISE.out.versions.first())

    //
    // MODULE: PURPOSE?
    //
    //NB! RENAME PROCESS
    // Need to join input channels to keep the same sample and genes together. Join on the meta map.
    // ch_parse_phylogeny_2 = ch_parse_phylogeny_temp.join(MAFFT_PAIRWISE.out.fas)
    ch_parse_phylogeny_2 = MAFFT_PAIRWISE.out.fas.join(PARSE_PHYLOGENY.out.ratio)
    PARSE_PHYLOGENY_2(
        ch_parse_phylogeny_2
    )
    ch_versions = ch_versions.mix(PARSE_PHYLOGENY_2.out.versions.first())

    CREATE_JPG(
        IQTREE.out.phylogeny
    )
    ch_versions = ch_versions.mix(CREATE_JPG.out.versions.first())

    emit:
    genotype = PARSE_PHYLOGENY_2.out.genotyping
    ratio    = PARSE_PHYLOGENY.out.ratio
    header   = PARSE_PHYLOGENY.out.header
    fasta    = PARSE_PHYLOGENY.out.percentcalc_fasta
    aligned  = MAFFT.out.fas
    versions = ch_versions

}
