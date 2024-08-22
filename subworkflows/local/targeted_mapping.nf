include { BOWTIE2_BUILD     } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN     } from '../../modules/nf-core/bowtie2/align/main'
include { TANOTI_ALIGN      } from '../../modules/local/tanoti.nf'
include { SAMTOOLS_INDEX as INDEX_WITHDUP } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as INDEX_MARKDUP } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_IDXSTATS } from '../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_DEPTH    } from '../../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_STATS as STATS_WITHDUP } from '../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_STATS as STATS_MARKDUP } from '../../modules/nf-core/samtools/stats/main'
include { IVAR_CONSENSUS } from '../../modules/nf-core/ivar/consensus/main'

include { BAM_MARKDUPLICATES_SAMTOOLS } from '../nf-core/bam_markduplicates_samtools/main'

workflow TARGETED_MAPPING {
    take:
    ch_major_mapping    //tuple val(meta), path(fasta), path(reads)

    main:
    // Extract the val(meta), path(fasta)
    ch_build = ch_major_mapping
        .map { meta, fasta, reads ->
        return [meta, fasta]
    }

    // Extract only the fasta path
    ch_fasta = ch_major_mapping
        .map { meta, fasta, reads ->
        return [fasta]
    }

    // Extract val(meta), path(reads). Add the reference name to the meta map
    ch_align = ch_major_mapping
        .map { meta, fasta, reads ->
            new_meta = meta + [ reference: fasta.getBaseName().toString().split('\\.').last() ]
        return [new_meta, reads]
    }

    if (params.mapper == "bowtie2") {
        BOWTIE2_BUILD (
            ch_build // val(meta), path(fasta)
        )
        BOWTIE2_ALIGN (
            ch_align,
            BOWTIE2_BUILD.out.index,
            false, // Do not save unmapped reads
            true // Sort bam file
        )
        ch_aligned = BOWTIE2_ALIGN.out.aligned
        ch_versions = BOWTIE2_ALIGN.out.versions // channel: [ versions.yml ]
    }
    else if (params.mapper == "tanoti") {
        TANOTI_ALIGN (
            ch_align,
            ch_build,
            true, // Sort bam file
            params.tanoti_stringency_2
        )
        ch_aligned = TANOTI_ALIGN.out.aligned
        ch_versions = TANOTI_ALIGN.out.versions // channel: [ versions.yml ]
    }

    // Generate stats file with duplicates included
    INDEX_WITHDUP (
        ch_aligned
    )

    STATS_WITHDUP (
        ch_aligned.join(INDEX_WITHDUP.out.bai), // val(meta), path(bam), path(bai)
        ch_build // val(meta), path(fasta)

    )
    // Remove duplicate reads
    BAM_MARKDUPLICATES_SAMTOOLS (
        ch_aligned,
        ch_major_mapping.map { meta, fasta, _ -> [meta, fasta] } // Extract only the meta and fasta elements from the channel
    )
    ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_SAMTOOLS.out.versions.first())

    //
    // MODULE: Identify the two references with most mapped reads
    //
    INDEX_MARKDUP (
        BAM_MARKDUPLICATES_SAMTOOLS.out.bam
    )
    ch_versions = ch_versions.mix(INDEX_MARKDUP.out.versions.first())

    SAMTOOLS_IDXSTATS (
        BAM_MARKDUPLICATES_SAMTOOLS.out.bam.join(INDEX_MARKDUP.out.bai) // val(meta), path(bam), path(bai)
    )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions.first())

    SAMTOOLS_DEPTH (
        BAM_MARKDUPLICATES_SAMTOOLS.out.bam,
        [ [], []] // Passing empty channels instead of an interval file
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    STATS_MARKDUP (
        BAM_MARKDUPLICATES_SAMTOOLS.out.bam.join(INDEX_MARKDUP.out.bai), // val(meta), path(bam), path(bai)
        ch_build // val(meta), path(fasta)
    )
    ch_versions = ch_versions.mix(STATS_MARKDUP.out.versions.first())

    //
    // MODULE: Create consensus sequence
    //
    IVAR_CONSENSUS(
        BAM_MARKDUPLICATES_SAMTOOLS.out.bam,
        ch_fasta,
        false // Don't need the mpileup file
    )
    ch_versions = ch_versions.mix(IVAR_CONSENSUS.out.versions.first())

    emit:
    aligned = BAM_MARKDUPLICATES_SAMTOOLS.out.bam
    versions = ch_versions
    depth = SAMTOOLS_DEPTH.out.tsv
    stats_withdup = STATS_WITHDUP.out.stats
    stats_markdup = STATS_MARKDUP.out.stats
    consensus = IVAR_CONSENSUS.out.fasta
}
