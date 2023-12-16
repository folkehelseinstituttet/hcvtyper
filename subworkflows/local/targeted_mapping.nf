include { BOWTIE2_BUILD     } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN     } from '../../modules/nf-core/bowtie2/align/main'
include { TANOTI_ALIGN      } from '../../modules/local/tanoti.nf'
include { SAMTOOLS_INDEX    } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_IDXSTATS } from '../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_DEPTH    } from '../../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_STATS    } from '../../modules/nf-core/samtools/stats/main'

include { BAM_MARKDUPLICATES_SAMTOOLS } from '../nf-core/bam_markduplicates_samtools/main'

workflow TARGETED_MAPPING {
    take:
    ch_major_mapping    //tuple val(meta), path(fasta), path(reads)
    prefix2

    main:
    // Extract the val(meta), path(fasta)
    ch_build = ch_major_mapping
        .map { meta, fasta, reads ->
        return [meta, fasta]
        }

    // Extract val(meta), path(reads). Add the reference name to the meta map
    ch_align = ch_major_mapping
        .map {
        meta, fasta, reads ->
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

    // Remove duplicate reads
    BAM_MARKDUPLICATES_SAMTOOLS (
        ch_aligned,
        ch_major_mapping.map { it[1] }
    )
    ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_SAMTOOLS.out.versions.first())

    //
    // MODULE: Identify the two references with most mapped reads
    //
    SAMTOOLS_INDEX (
        BAM_MARKDUPLICATES_SAMTOOLS.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
    SAMTOOLS_IDXSTATS (
        BAM_MARKDUPLICATES_SAMTOOLS.out.bam.join(SAMTOOLS_INDEX.out.bai) // val(meta), path(bam), path(bai)
    )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions.first())

    SAMTOOLS_DEPTH (
        BAM_MARKDUPLICATES_SAMTOOLS.out.bam,
        [ [], []] // Passing empty channels instead of an interval file
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())
    SAMTOOLS_STATS (
        BAM_MARKDUPLICATES_SAMTOOLS.out.bam.join(SAMTOOLS_INDEX.out.bai), // val(meta), path(bam), path(bai)
        ch_build // val(meta), path(fasta)
    )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())
        // stats_markdup = BOWTIE2_ALIGN.out.stats_markdup
        // stats_withdup = BOWTIE2_ALIGN.out.stats_withdup
    emit:
    aligned = ch_aligned
    versions = ch_versions
    depth = SAMTOOLS_DEPTH.out.tsv
    stats_markdup = SAMTOOLS_STATS.out.stats
    // stats_withdup

}
