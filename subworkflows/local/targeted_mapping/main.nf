//
// Map reads against major or minor reference, get mapping stats, create consensus and create QC plots.
//

include { BOWTIE2_BUILD                   } from '../../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN                   } from '../../../modules/nf-core/bowtie2/align/main'
include { PLOTCOVERAGE                    } from '../../../modules/local/plotcoverage/main'
include { TANOTI_ALIGN                    } from '../../../modules/local/tanoti.nf'
include { SAMTOOLS_INDEX as INDEX_WITHDUP } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as INDEX_MARKDUP } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_IDXSTATS               } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_DEPTH                  } from '../../../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_STATS as STATS_WITHDUP } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_STATS as STATS_MARKDUP } from '../../../modules/nf-core/samtools/stats/main'
include { IVAR_CONSENSUS                  } from '../../../modules/nf-core/ivar/consensus/main'
include { PLOT_BAMVARIATION               } from '../../../modules/local/bamvariation'
include { SAMTOOLS_SORMADUP               } from '../../../modules/nf-core/samtools/sormadup/main'

workflow TARGETED_MAPPING {

    take:
    ch_major_mapping    //tuple val(meta), path(fasta), path(reads)

    main:

        // Enrich the meta map with the reference name before splitting the channel.
        // This ensures all branches (build, fasta, reads) share the same meta key so
        // that downstream joins and positional pairings are always in sync.
        ch_input = ch_major_mapping
        .map { meta, fasta, reads ->
            def new_meta = meta + [ reference: fasta.getBaseName().toString().split('\\.').last() ]
            tuple(new_meta, fasta, reads)
        }
        .multiMap { meta, fasta, reads ->
            build: [ meta, fasta ]
            fasta: [ fasta ]       // consumed by IVAR_CONSENSUS (positional, no meta needed)
            reads: [ meta, reads ]
        }

    if (params.mapper == "bowtie2") {
        BOWTIE2_BUILD (
            ch_input.build // val(meta), path(fasta)
        )

        // Join reads, index, and fasta by meta key before calling BOWTIE2_ALIGN.
        // BOWTIE2_BUILD emits index items in completion order (not submission order),
        // so a positional join would pair the wrong index with the wrong sample when
        // multiple samples are processed in parallel. Joining by meta key guarantees
        // that each BOWTIE2_ALIGN task always receives the correct matched triple.
        ch_aligned_input = ch_input.reads       // meta, reads
            .join( BOWTIE2_BUILD.out.index )   // meta, reads, index
            .join( ch_input.build )            // meta, reads, index, fasta

        BOWTIE2_ALIGN (
            ch_aligned_input.map { meta, reads, _index, _fasta -> [ meta, reads ] },
            ch_aligned_input.map { meta, _reads, index, _fasta -> [ meta, index ] },
            ch_aligned_input.map { meta, _reads, _index, fasta -> [ meta, fasta ] },
            false, // Do not save unmapped reads
            true   // Sort bam file
        )
        ch_aligned = BOWTIE2_ALIGN.out.bam
        ch_versions = BOWTIE2_ALIGN.out.versions // channel: [ versions.yml ]
    }
    else if (params.mapper == "tanoti") {
        TANOTI_ALIGN (
            ch_input.reads, // tuple val(meta), path(reads) — fasta is no longer bundled in here
            ch_input.build,
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
        ch_input.build // val(meta), path(fasta)

    )
    // Remove duplicate reads
    SAMTOOLS_SORMADUP (
        ch_aligned,
        ch_input.build // val(meta), path(fasta)
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORMADUP.out.versions.first())

    //
    // MODULE: Identify the two references with most mapped reads
    //
    INDEX_MARKDUP (
        SAMTOOLS_SORMADUP.out.bam
    )
    ch_versions = ch_versions.mix(INDEX_MARKDUP.out.versions.first())

    SAMTOOLS_IDXSTATS (
        SAMTOOLS_SORMADUP.out.bam.join(INDEX_MARKDUP.out.bai) // val(meta), path(bam), path(bai)
    )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions.first())

    SAMTOOLS_DEPTH (
        SAMTOOLS_SORMADUP.out.bam,
        [ [], []] // Passing empty channels instead of an interval file
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    STATS_MARKDUP (
        SAMTOOLS_SORMADUP.out.bam.join(INDEX_MARKDUP.out.bai), // val(meta), path(bam), path(bai)
        ch_input.build // val(meta), path(fasta)
    )
    ch_versions = ch_versions.mix(STATS_MARKDUP.out.versions.first())

    //
    // MODULE: Plot coverage from mapping
    //
    PLOTCOVERAGE (
        SAMTOOLS_DEPTH.out.tsv
    )
    ch_versions = ch_versions.mix(PLOTCOVERAGE.out.versions)

    //
    // MODULE: Plot variation in mapping file
    //
    PLOT_BAMVARIATION (
        SAMTOOLS_SORMADUP.out.bam
    )
    ch_versions = ch_versions.mix(PLOT_BAMVARIATION.out.versions)

    //
    // MODULE: Create consensus sequence
    //
    IVAR_CONSENSUS(
        SAMTOOLS_SORMADUP.out.bam,
        ch_input.fasta,
        false // Don't need the mpileup file
    )
    ch_versions = ch_versions.mix(IVAR_CONSENSUS.out.versions.first())

    emit:
    aligned = SAMTOOLS_SORMADUP.out.bam
    depth = SAMTOOLS_DEPTH.out.tsv
    stats_withdup = STATS_WITHDUP.out.stats
    stats_markdup = STATS_MARKDUP.out.stats
    consensus = IVAR_CONSENSUS.out.fasta
    variation = PLOT_BAMVARIATION.out.png

    versions = ch_versions                     // channel: [ versions.yml ]
}
