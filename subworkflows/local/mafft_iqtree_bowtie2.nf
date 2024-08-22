//
// Align with mafft, create phylogeny with iqtree, map with Bowtie2, sort, stats and markdup using Samtools
//

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
include { SAMTOOLS_DEPTH                       } from '../../modules/nf-core/samtools/depth/main'
include { IQTREE                               } from '../../modules/nf-core/iqtree/main'
include { MAFFT                                } from '../../modules/nf-core/mafft/main'
include { MAFFT as MAFFT_PAIRWISE              } from '../../modules/nf-core/mafft/main'
include { PARSE_PHYLOGENY                      } from '../../modules/local/parse_phylogeny'
include { EXTRACT_COMBINE_SEQS                 } from '../../modules/local/extract_combine_seqs'
include { COLLECT_GENOTYPE_INFO                } from '../../modules/local/collect_genotype_info'
include { PREPARE_MAFFT                        } from '../../modules/local/prepare_mafft'
include { SUMMARIZE_STATS as SUMMARIZE_STATS_WITHDUP } from '../../modules/local/summarize_stats'
include { SUMMARIZE_STATS as SUMMARIZE_STATS_MARKDUP } from '../../modules/local/summarize_stats'
include { SUMMARIZE_DEPTH } from '../../modules/local/summarize_depth'

include { BAM_MARKDUPLICATES_SAMTOOLS          } from '../../subworkflows/nf-core/bam_markduplicates_samtools/main'

workflow MAFFT_IQTREE_BOWTIE2 {

    take:
    ch_vigorparse   // channel: [ val(meta), path(gene_fasta) ]
    ch_classified_reads
    ch_references // channel: [ val(meta), path(gene_references) ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Combine all contigs from a given segment with the corresponding reference dataset for MAFFT input
    //
    PREPARE_MAFFT(
        ch_vigorparse,
        ch_references
    )
    ch_versions = ch_versions.mix(PREPARE_MAFFT.out.versions.first())

    //
    // MODULE: Align gene sequences with MAFFT
    //

    // NOTE:
    // First transpose the PREPARE_MAFFT output to split the concatenated individual gene fasta files into different channels.
    // This is to be able to parallelize the MAFFT and IQTREE steps for each gene.
    // Then prepare a mafft input channel that has the gene name in the meta map. This is for renaming the output files.
    // Add the gene name from the fasta file names to the meta map like this: [[id:, single_end:, gene:], file path]
    ch_mafft = PREPARE_MAFFT.out.fasta.transpose()
        .map { meta, fasta ->

        def filePath = fasta.toString() // File path as string

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

    //
    // MODULE: IQTREE
    //
    IQTREE(
        MAFFT.out.fas.map { meta -> return meta + [[]] }, // Add empty element to MAFFT.out.fas to comply with IQTREE input
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

    //
    // MODULE: Parse the phylogeny to identify the contigs, their nearest referencees and clade identity. The process can handle multiple contigs in the tree
    //
    PARSE_PHYLOGENY(
        IQTREE.out.phylogeny
    )
    ch_versions = ch_versions.mix(PARSE_PHYLOGENY.out.versions.first())

    //
    // MODULE: Extract and combine sequences
    //

    // NOTE:
    // Add the gene name to the meta of ch_vigorparse ([id:, single_end:, gene:])
    // First transpose the channel to separate the different fasta files for each segment from extract from gff.
    //ch_temp_2 = ch_vigorparse.transpose()
    //ch_extract_combine_seqs_temp = ch_temp_2.map { item ->
    ch_extract_combine_seqs_temp = ch_vigorparse.transpose()
        .map { meta, contigs ->
        //def meta = item[0] // Original meta map
        //def filePath = item[1].toString() // File path as string
        def filePath = contigs.toString() // File path as string

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

    // NOTE:
    // Then join the ch_parse_phylogeny with the PARSE_PHYLOGENY output to have the phylogeny and the segment fasta for the same gene
    // PARSE_PHYLOGENY.out.parse_phylo has this structure: val(id:, single_end:, gene:), path(gff_extract_fasta for single gene)
    EXTRACT_COMBINE_SEQS(
        ch_extract_combine_seqs_temp.join(PARSE_PHYLOGENY.out.parse_phylo), // val(meta), path(gff_extract_fasta), path(parse_phylo),
        ch_references
    )
    ch_versions = ch_versions.mix(EXTRACT_COMBINE_SEQS.out.versions.first())

    //
    // MODULE: Align each contig to its nearest reference sequence in the tree using mafft
    //

    // NOTE:
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
    ch_versions = ch_versions.mix(CALCULATE_PAIRWISE_ALIGNMENT_METRICS.out.versions.first())

    //
    // MODULE: Map reads to the different contigs.
    //

    // NOTE:
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
    ch_versions = ch_versions.mix(PREPARE_BOWTIE2_BUILD.out.versions.first())

    //
    // MODULE: Create a bowtie2 index
    //
    BOWTIE2_BUILD(
        PREPARE_BOWTIE2_BUILD.out.contig
    )
    ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions.first())

    //
    // MODULE: Map reads to the contigs with BOWTIE2
    //

    // NOTE:
    // Ensure that classified reads are from the same sample as Bowtie2 build output
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
        }.set { ch_bowtie2_align } // Output into channel ch_bowtie2_align

    BOWTIE2_ALIGN (
        ch_bowtie2_align.reads,
        ch_bowtie2_align.index,
        false, // Do not save unmapped reads
        true // Sort bam file
    )
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

    //
    // MODULE: Index the bam file with samtools index
    //
    INDEX_WITHDUP (
        BOWTIE2_ALIGN.out.aligned
    )
    ch_versions = ch_versions.mix(INDEX_WITHDUP.out.versions.first())

    //
    // MODULE: Run samtools stats on the bam files with duplicates included
    //

    // NOTE:
    // Ensure no sample mixup for the STATS module. Merging output from PREPARE_BOWTIE2_BUILD, BOWTIE2_ALIGN and INDEX_WITHDUP
    PREPARE_BOWTIE2_BUILD.out.contig.map {
        meta, contig -> [
            meta.subMap( ['id','single_end'] ), contig // Keep only "id" and "single_end" keys for joining
        ]
        }
        .join(BOWTIE2_ALIGN.out.aligned, by: 0) // Join on meta with the bam file from BOWTIE2_ALIGN
        .join(INDEX_WITHDUP.out.bai, by: 0) // Join on meta with the bai file from INDEX_WITHDUP
        // Split into two channels for input into the STATS module
        .multiMap { meta, contig, bam, bai ->
            bam_bai: [ meta, bam, bai ]
            contig: [ meta, contig ]
        }.set { ch_stats }

    STATS_WITHDUP (
        ch_stats.bam_bai, // val(meta), path(bam), path(bai)
        ch_stats.contig // val(meta), path(fasta)
    )
    ch_versions = ch_versions.mix(STATS_WITHDUP.out.versions.first())

    //
    // MODULE: Summarize the samtools stats output
    //
    SUMMARIZE_STATS_WITHDUP (
        STATS_WITHDUP.out.stats
    )
    ch_versions = ch_versions.mix(SUMMARIZE_STATS_WITHDUP.out.versions.first())

    //
    // MODULE: Remove duplicate reads
    //

    //
    // SUBWORKFLOW: Collate, fixmate, sort and markdup using Samtools
    //

    // NOTE:
    // Ensure that the two input channels are from the same sample
    PREPARE_BOWTIE2_BUILD.out.contig.map {
        meta, contig -> [
            meta.subMap( ['id','single_end'] ), contig // Keep only "id" and "single_end" keys for joining
        ]
        }
        .join(BOWTIE2_ALIGN.out.aligned, by: 0) // Join on meta with the bam file from BOWTIE2_ALIGN
        // Split into two channels for input into the BAM_MARKDUPLICATES_SAMTOOLS subworkflow
        .multiMap { meta, contig, bam ->
            bam: [ meta, bam ]
            contig: [ meta, contig ]
        }.set { ch_bam_markdup_samtools }

    BAM_MARKDUPLICATES_SAMTOOLS(
        ch_bam_markdup_samtools.bam,
        ch_bam_markdup_samtools.contig
    )
    ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_SAMTOOLS.out.versions.first())

    //
    // MODULE: Index the bam file with samtools index
    //
    INDEX_MARKDUP (
        BAM_MARKDUPLICATES_SAMTOOLS.out.bam
    )
    ch_versions = ch_versions.mix(INDEX_MARKDUP.out.versions.first())

    //
    // MODULE: Get reads coverage per position with samtools depth
    //
    SAMTOOLS_DEPTH (
        BAM_MARKDUPLICATES_SAMTOOLS.out.bam,
        [ [], []] // Passing empty channels instead of an interval file
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    //
    // MODULE: Summarize samtools depth without duplicates
    //
    SUMMARIZE_DEPTH(
        SAMTOOLS_DEPTH.out.tsv
    )
    ch_versions = ch_versions.mix(SUMMARIZE_DEPTH.out.versions.first())

    // NOTE:
    // Ensure no sample mixup for the STATS module. Merging output from PREPARE_BOWTIE2_BUILD, BAM_MARKDUPLICATES_SAMTOOLS and INDEX_MARKDUP
    PREPARE_BOWTIE2_BUILD.out.contig.map {
        meta, contig -> [
            meta.subMap( ['id','single_end'] ), contig // Keep only "id" and "single_end" keys for joining
        ]
        }
        .join(BAM_MARKDUPLICATES_SAMTOOLS.out.bam, by: 0) // Join on meta with the bam file from BOWTIE2_ALIGN
        .join(INDEX_MARKDUP.out.bai, by: 0) // Join on meta with the bai file from INDEX_WITHDUP
        // Split into two channels for input into the STATS module
        .multiMap { meta, contig, bam, bai ->
            bam_bai: [ meta, bam, bai ]
            contig: [ meta, contig ]
        }.set { ch_stats_markdup }

    STATS_MARKDUP (
        ch_stats_markdup.bam_bai, // val(meta), path(bam), path(bai)
        ch_stats_markdup.contig // val(meta), path(fasta)
    )
    ch_versions = ch_versions.mix(STATS_MARKDUP.out.versions.first())

    //
    // MODULE: Summarize the samtools stats output
    //
    SUMMARIZE_STATS_MARKDUP (
        STATS_MARKDUP.out.stats
    )
    ch_versions = ch_versions.mix(SUMMARIZE_STATS_MARKDUP.out.versions.first())


    emit:
    alignment_metrics = CALCULATE_PAIRWISE_ALIGNMENT_METRICS.out.metrics
    parse_phylo       = PARSE_PHYLOGENY.out.parse_phylo
    fasta             = EXTRACT_COMBINE_SEQS.out.combined_fasta
    aligned           = MAFFT.out.fas
    stats_withdup     = SUMMARIZE_STATS_WITHDUP.out.csv
    stats_markdup     = SUMMARIZE_STATS_MARKDUP.out.csv
    bam_nodups        = BAM_MARKDUPLICATES_SAMTOOLS.out.bam
    depth             = SAMTOOLS_DEPTH.out.tsv
    versions          = ch_versions
}
