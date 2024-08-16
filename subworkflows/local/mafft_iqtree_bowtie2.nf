include { PREPARE_BOWTIE2_BUILD                } from '../../modules/local/prepare_bowtie2_build'
include { BOWTIE2_BUILD                        } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN                        } from '../../modules/nf-core/bowtie2/align/main'
include { CREATE_JPG                           } from '../../modules/local/create_jpg'
include { CALCULATE_PAIRWISE_ALIGNMENT_METRICS } from '../../modules/local/calculate_pairwise_alignment_metrics'
include { IQTREE                               } from '../../modules/nf-core/iqtree/main'
include { MAFFT                                } from '../../modules/nf-core/mafft/main'
include { MAFFT as MAFFT_PAIRWISE              } from '../../modules/nf-core/mafft/main'
include { PARSE_PHYLOGENY                      } from '../../modules/local/parse_phylogeny'
include { EXTRACT_COMBINE_SEQS                 } from '../../modules/local/extract_combine_seqs'
include { COLLECT_GENOTYPE_INFO                } from '../../modules/local/collect_genotype_info'
include { PREPARE_MAFFT                        } from '../../modules/local/prepare_mafft'

workflow MAFFT_IQTREE_BOWTIE2 {

    take:
    ch_vigorparse   // channel: [ val(meta), path(gene_fasta) ]
    ch_references // channel: [ val(meta), path(gene_references) ]
    ch_reads

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Combine all contigs from a given segment with the corresponding reference dataset for MAFFT input
    //
    //
    // MODULE: Combine the highest covered gene contig with the corresponding reference dataset for MAFFT input
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

    //
    // MODULE: Create a jpg image of the phylogeny
    //
    CREATE_JPG(
        IQTREE.out.phylogeny
    )
    ch_versions = ch_versions.mix(CREATE_JPG.out.versions.first())

    // Here we need to join all the output from a single sample before further parsing

    //
    // MODULE: Parse the phylogeny to identify the contigs, their nearest referencees and clade identity.
    // The process can handle multiple contigs in the tree
    //

    // Channel only needs meta and treefile. The meta has the gene name added at this point.
    PARSE_PHYLOGENY(
        IQTREE.out.phylogeny
    )
    ch_versions = ch_versions.mix(PARSE_PHYLOGENY.out.versions.first())

    //
    // MODULE: Extract and combine sequences
    //

    // Add the gene name to the meta of ch_vigorparse ([id:, single_end:, gene:])
    // First transpose the channel to separate the different fasta files for each segment from extract from gff.
    ch_temp_2 = ch_vigorparse.transpose()
    ch_extract_combine_seqs_temp = ch_temp_2.map { item ->
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

    // Then join the ch_parse_phylogeny with the PARSE_PHYLOGENY output to have the phylogeny and the segment fasta for the same gene
    // PARSE_PHYLOGENY.out.parse_phylo has this structure: val(id:, single_end:, gene:), path(gff_extract_fasta for single gene)
    ch_extract_combine_seqs = ch_extract_combine_seqs_temp.join(PARSE_PHYLOGENY.out.parse_phylo) // val(meta), path(gff_extract_fasta), path(parse_phylo)

    EXTRACT_COMBINE_SEQS(
        ch_extract_combine_seqs,
        ch_references
    )

    //
    // MODULE: Align each contig to its nearest reference sequence in the tree
    //

    // EXTRACT_COMBINE_SEQS.out.combined_fasta can contain multiple combined fasta files if there are several contigs matching the same gene.
    // THE MAFFT module needs only one fasta file at a time, so we need to transpose the channel to separate the different fasta files.
    ch_mafft_pairwise = EXTRACT_COMBINE_SEQS.out.combined_fasta.transpose()
    MAFFT_PAIRWISE(
        ch_mafft_pairwise,
        [ [:], [] ],
        [ [:], [] ],
        [ [:], [] ],
        [ [:], [] ],
        [ [:], [] ],
        false
    )
    ch_versions = ch_versions.mix(MAFFT_PAIRWISE.out.versions.first())

    //
    // MODULE: CALCULATE PAIRWISE ALIGNMENT METRICS
    //
    CALCULATE_PAIRWISE_ALIGNMENT_METRICS(
        MAFFT_PAIRWISE.out.fas
    )

    //
    // MODULE: Map reads to the different contigs.
    //

    PREPARE_BOWTIE2_BUILD(
        ch_mafft_pairwise
    )
    // ch_build = ch_mafft_pairwise
    //     .splitFasta(record: [header: true, seqString: true]) // Split fasta into records
    //     .filter { meta, record -> record.header =~ /^NODE.*/ }
    // ch_build.view()
    BOWTIE2_BUILD(
        PREPARE_BOWTIE2_BUILD.out.contig
    )

    // BOWTIE2_ALIGN(
    //
    // )


    emit:
    alignment_metrics = CALCULATE_PAIRWISE_ALIGNMENT_METRICS.out.metrics
    parse_phylo       = PARSE_PHYLOGENY.out.parse_phylo
    fasta             = EXTRACT_COMBINE_SEQS.out.combined_fasta
    aligned           = MAFFT.out.fas
    versions          = ch_versions

}
