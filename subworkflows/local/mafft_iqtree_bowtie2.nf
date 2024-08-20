include { PREPARE_BOWTIE2_BUILD                } from '../../modules/local/prepare_bowtie2_build'
include { JOIN_CONTIGS                         } from '../../modules/local/join_contigs'
include { BOWTIE2_BUILD                        } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN                        } from '../../modules/nf-core/bowtie2/align/main'
include { CREATE_JPG                           } from '../../modules/local/create_jpg'
include { CALCULATE_PAIRWISE_ALIGNMENT_METRICS } from '../../modules/local/calculate_pairwise_alignment_metrics'
include { SAMTOOLS_INDEX as INDEX_WITHDUP      } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as INDEX_MARKDUP      } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS as STATS_WITHDUP      } from '../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_STATS as STATS_MARKDUP      } from '../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_IDXSTATS as SAMTOOLS_IDXSTATS_WITHDUP                    } from '../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_IDXSTATS as SAMTOOLS_IDXSTATS_MARKDUP                    } from '../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_DEPTH                       } from '../../modules/nf-core/samtools/depth/main'
include { IQTREE                               } from '../../modules/nf-core/iqtree/main'
include { MAFFT                                } from '../../modules/nf-core/mafft/main'
include { MAFFT as MAFFT_PAIRWISE              } from '../../modules/nf-core/mafft/main'
include { PARSE_PHYLOGENY                      } from '../../modules/local/parse_phylogeny'
include { EXTRACT_COMBINE_SEQS                 } from '../../modules/local/extract_combine_seqs'
include { COLLECT_GENOTYPE_INFO                } from '../../modules/local/collect_genotype_info'
include { PREPARE_MAFFT                        } from '../../modules/local/prepare_mafft'
include { SUMMARIZE_IDXSTATS as SUMMARIZE_IDXSTATS_WITHDUP } from '../../modules/local/summarize_idxstats'
include { SUMMARIZE_IDXSTATS as SUMMARIZE_IDXSTATS_MARKDUP } from '../../modules/local/summarize_idxstats'

include { BAM_MARKDUPLICATES_SAMTOOLS          } from '../../subworkflows/nf-core/bam_markduplicates_samtools/main'

workflow MAFFT_IQTREE_BOWTIE2 {

    take:
    ch_vigorparse   // channel: [ val(meta), path(gene_fasta) ]
    ch_references // channel: [ val(meta), path(gene_references) ]
    ch_classified_reads

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
    // This process simply extracts the contig fasta file from the fasta file that comes out of the pairwise mafft which also includes the reference sequence.
    // We could use groovy code, but this will enter the fasta file as a string and not a file.
    // It's also possible to save as a file, but then it's hard to control the file names as I like to use the file names to keep track of samples and genes.  And also to collect files later.
    //    ch_mafft_pairwise
    //        .splitFasta(record: [header: true, seqString: true]) // Split fasta into records
    //        .filter { meta, record -> record.header =~ /^NODE.*/ }
    //        .collectFile(name: 'contig.fasta', newLine: true) { ">${it.get(1).header}\n${it.get(1).seqString}"}
    PREPARE_BOWTIE2_BUILD(
        ch_mafft_pairwise
    )

    // Create a process that will take all the fasta files per sample id. Then simply concatenate them and output together with the sample id.
//    JOIN_CONTIGS(
//        // Collect all the contig fastas per sample
//        PREPARE_BOWTIE2_BUILD.out.contig
//            .map {
//                meta, fasta -> [
//                    meta.subMap( ['id'] ), fasta // Keep only id. I.e. remove the "single_end" and "gene" keys from the meta map
//                    ]
//           }
//            .groupTuple(by: 0) // Group by the meta map which only holds the sample id
//    )

    BOWTIE2_BUILD(
        //JOIN_CONTIGS.out.contig
        PREPARE_BOWTIE2_BUILD.out.contig
    )

    // TODO: Ensure that classified reads are from the same sample as Bowtie2 build output
    // TODO: It's probably better to map to each contig separately. This is because one read may map to several contigs.
    // TODO: Create two input channels by splitting after joining classified reads and bowtie2 index output

    // Can  I join on the sample id only
    // First remove the keys "single_end" and "gene" from the meta map of BOWTIE2_BUILD.out.index, then join on the sample id
//    BOWTIE2_BUILD.out.index.map {
//        meta, path -> [
//            meta.subMap( ['id','single_end'] ), path // Keep only "id" and "single_end" keys for joining
//        ]
//    }.join(ch_classified_reads) // Join with the classified reads channel on the meta map
//    .multiMap { meta, index, reads ->
//        index: [ meta, index ]
//        reads: [ meta, reads ]
//    }.set { ch_bowtie2_align }

    //ch_bowtie2_align.index.view()
    //ch_bowtie2_align.reads.view()

    // Combine the classified reads with the bowtie2 index output.
    // This is a few to many combination, as the classified reads channel (one per sample) needs to be combined with many contig channels
    BOWTIE2_BUILD.out.index.map {
        meta, path -> [
            meta.subMap( ['id','single_end'] ), path // Keep only "id" and "single_end" keys for joining
        ]
        }.combine(ch_classified_reads, by: 0)
        .multiMap { meta, index, reads ->
            index: [ meta, index ]
            reads: [ meta, reads ]
        }.set { ch_bowtie2_align }

    // Now define two new channels: one for the classified reads and one for the contigs. They should both contain the meta map


    BOWTIE2_ALIGN (
        ch_bowtie2_align.reads,
        ch_bowtie2_align.index,
        false, // Do not save unmapped reads
        true // Sort bam file
    )
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

    // Generate stats file with duplicates included
    INDEX_WITHDUP (
        BOWTIE2_ALIGN.out.aligned
    )

   // STATS_WITHDUP (
    //    BOWTIE2_ALIGN.out.aligned.join(INDEX_WITHDUP.out.bai), // val(meta), path(bam), path(bai)
    //    PREPARE_BOWTIE2_BUILD.out.contig // val(meta), path(fasta)
   // )
    SAMTOOLS_IDXSTATS_WITHDUP (
        BOWTIE2_ALIGN.out.aligned.join(INDEX_WITHDUP.out.bai) // val(meta), path(bam), path(bai)
    )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS_WITHDUP.out.versions.first())

    // Summarize IDXSTATS OUTPUT HERE:
    // SHOULD BE STATS FILE MAYBE? SAMTOOLS STATS?
    //SUMMARIZE_IDXSTATS_WITHDUP (
    //    SAMTOOLS_IDXSTATS_WITHDUP.out.idxstats
    //)

    // Remove duplicate reads
    BAM_MARKDUPLICATES_SAMTOOLS(
        BOWTIE2_ALIGN.out.aligned,
        PREPARE_BOWTIE2_BUILD.out.contig
    )
    ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_SAMTOOLS.out.versions.first())

    INDEX_MARKDUP (
        BAM_MARKDUPLICATES_SAMTOOLS.out.bam
    )
    ch_versions = ch_versions.mix(INDEX_MARKDUP.out.versions.first())

    SAMTOOLS_IDXSTATS_MARKDUP (
        BAM_MARKDUPLICATES_SAMTOOLS.out.bam.join(INDEX_MARKDUP.out.bai) // val(meta), path(bam), path(bai)
    )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS_MARKDUP.out.versions.first())

    // Summarize IDXSTATS OUTPUT HERE. CREATE OUTPUT FILES THAT HAS THE SAMPLE NAME, GENE NAME, CONTIG NAME AND REFERENCE NAME
    //SUMMARIZE_IDXSTATS_MARKDUP (
    //    SUMMARIZE_IDXSTATS_MARKDUP.out.idxstats
    //)

    SAMTOOLS_DEPTH (
        BAM_MARKDUPLICATES_SAMTOOLS.out.bam,
        [ [], []] // Passing empty channels instead of an interval file
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    STATS_MARKDUP (
        BAM_MARKDUPLICATES_SAMTOOLS.out.bam.join(INDEX_MARKDUP.out.bai), // val(meta), path(bam), path(bai)
        PREPARE_BOWTIE2_BUILD.out.contig // val(meta), path(fasta)
    )
    ch_versions = ch_versions.mix(STATS_MARKDUP.out.versions.first())



    emit:
    alignment_metrics = CALCULATE_PAIRWISE_ALIGNMENT_METRICS.out.metrics
    parse_phylo       = PARSE_PHYLOGENY.out.parse_phylo
    fasta             = EXTRACT_COMBINE_SEQS.out.combined_fasta
    aligned           = MAFFT.out.fas
    bam_nodups        = BAM_MARKDUPLICATES_SAMTOOLS.out.bam
    stats_withdup     = SAMTOOLS_IDXSTATS_WITHDUP.out.idxstats
    stats_markdup     = SAMTOOLS_IDXSTATS_MARKDUP.out.idxstats
    depth             = SAMTOOLS_DEPTH.out.tsv
    versions          = ch_versions
}
