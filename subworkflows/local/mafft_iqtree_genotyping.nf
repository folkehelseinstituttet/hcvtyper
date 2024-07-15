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

    // Prepare a mafft input channel that has the gene name in the meta map. This is for renaming the output files.
    // Add the gene name from the fasta file names to the meta map like this:
    // [[meta.id, meta.single_end, meta.gene], file path]
    ch_mafft = PREPARE_MAFFT.out.fasta.map { item ->
        def meta = item[0] // Original meta map
        def filePath = item[1].toString() // File path as string

        // Use getName to get the filename from the file path
        def fileName = new File(filePath).getName()

        // Extract the gene name from the file path
        def geneName = fileName.split("_")[0]

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

    //
    // MODULE: Parse phylogeny to genotype the sample sequence
    //
    PARSE_PHYLOGENY(
        IQTREE.out.phylogeny,
        ch_references,
        ch_vigorparse
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
    PARSE_PHYLOGENY_2(
        MAFFT_PAIRWISE.out.fas,
        PARSE_PHYLOGENY.out.ratio
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
