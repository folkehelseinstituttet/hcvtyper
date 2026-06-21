//
// Competitive joint mapping (replaces TARGETED_MAPPING).
//
// Per SAMPLE: concatenate all passing candidate FASTAs into one combined reference,
// build ONE Bowtie2 index, align all classified reads ONCE (Bowtie2 assigns each read to
// its single best location across the combined reference -> competitive partitioning, no
// conserved-region double-counting), then split the combined dedup BAM back into
// per-candidate BAMs that feed the existing per-candidate downstream steps unchanged.
//
// Read-count metrics now come from SAMTOOLS_IDXSTATS on the combined BAM, pre- and
// post-dedup (D-07/D-09). SAMTOOLS_STATS is removed (D-15) — summarize.R consumed only
// `reads mapped:`, which maps to idxstats column 3.
//

include { CAT_CANDIDATES                         } from '../../../modules/local/cat_candidates/main'
include { BOWTIE2_BUILD                          } from '../../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN                          } from '../../../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_INDEX as INDEX_COMBINED       } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as INDEX_DEDUP          } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_IDXSTATS as IDXSTATS_WITHDUP  } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_IDXSTATS as IDXSTATS_NODUP    } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_SORMADUP                      } from '../../../modules/nf-core/samtools/sormadup/main'
include { SAMTOOLS_VIEW                          } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_DEPTH                         } from '../../../modules/nf-core/samtools/depth/main'
include { IVAR_CONSENSUS                         } from '../../../modules/nf-core/ivar/consensus/main'
include { PLOTCOVERAGE                           } from '../../../modules/local/plotcoverage/main'
include { PLOT_BAMVARIATION                      } from '../../../modules/local/bamvariation'
include { CONSENSUS_DISTANCE                     } from '../../../modules/local/consensus_distance/main'

workflow JOINT_MAPPING {

    take:
    // ONE element PER SAMPLE (D-01): the collected passing candidate FASTAs, the classified
    // reads, and the candidates CSV used for the post-split per-candidate fan-out.
    ch_joint_mapping // tuple val(meta), path(cand_fastas), path(reads), path(candidates_csv)

    main:

    ch_versions = Channel.empty()

    // Split the per-sample take into the per-step shapes. cand_fastas + reads are needed for
    // the combined-mapping leg; candidates_csv (and the cand_fastas list) are carried to the
    // post-dedup fan-out so each candidate can pick its own reference FASTA by rank slot.
    ch_sample = ch_joint_mapping
        .multiMap { meta, cand_fastas, reads, candidates_csv ->
            cat:   [ meta, cand_fastas ]
            reads: [ meta, reads ]
            // Normalize cand_fastas to a list so the rank-indexed lookup in the fan-out is
            // uniform whether a sample has one FASTA (bare) or many (list).
            fanout: [ meta, candidates_csv, (cand_fastas instanceof List ? cand_fastas : [cand_fastas]) ]
        }

    //
    // MODULE: Concatenate all candidate FASTAs into one combined per-sample reference (D-04).
    //         cat of a single FASTA is identity (single-candidate degenerate case, D-14).
    //
    CAT_CANDIDATES (
        ch_sample.cat
    )
    ch_versions = ch_versions.mix(CAT_CANDIDATES.out.versions.first())

    //
    // MODULE: Build ONE combined Bowtie2 index per sample (D-05).
    //
    BOWTIE2_BUILD (
        CAT_CANDIDATES.out.fasta
    )
    ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions.first())

    //
    // MODULE: Competitive alignment — all reads mapped ONCE against the combined index (D-06).
    //         Bowtie2 default reporting (no -k/-a) gives one best alignment per read, so a read
    //         in a conserved region shared by two candidates lands at exactly one of them.
    //         save_unaligned=false + sort_bam=true exactly as TARGETED_MAPPING. The region split
    //         downstream naturally drops the unmapped `*` bin, so --no-unal is not required
    //         (Pitfall 5). ext.args (--very-sensitive-local) is supplied via config in Plan 02.
    //
    // Join reads + index + combined fasta by meta key (BOWTIE2_BUILD emits in completion order,
    // so a positional pairing would cross-pair indices across parallel samples).
    ch_aligned_input = ch_sample.reads                 // meta, reads
        .join( BOWTIE2_BUILD.out.index )               // meta, reads, index
        .join( CAT_CANDIDATES.out.fasta )              // meta, reads, index, combined.fa

    BOWTIE2_ALIGN (
        ch_aligned_input.map { meta, reads, _index, _fasta -> [ meta, reads ] },
        ch_aligned_input.map { meta, _reads, index, _fasta -> [ meta, index ] },
        ch_aligned_input.map { meta, _reads, _index, fasta -> [ meta, fasta ] },
        false, // Do not save unmapped reads
        true   // Sort bam file (required for index + region split)
    )
    ch_withdup = BOWTIE2_ALIGN.out.bam                 // combined sorted BAM (withdup)
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

    //
    // MODULE: Index the combined withdup BAM, then idxstats -> per-reference withdup counts (D-07).
    //
    INDEX_COMBINED (
        ch_withdup
    )
    ch_versions = ch_versions.mix(INDEX_COMBINED.out.versions.first())

    IDXSTATS_WITHDUP (
        ch_withdup.join(INDEX_COMBINED.out.bai)        // meta, bam, bai
    )
    ch_versions = ch_versions.mix(IDXSTATS_WITHDUP.out.versions.first())

    //
    // MODULE: Mark/remove duplicates on the combined BAM (D-08). Competitive mapping has already
    //         assigned each read to one candidate, so dedup is correctly scoped per reference
    //         sequence within the combined BAM.
    //
    SAMTOOLS_SORMADUP (
        ch_withdup,
        CAT_CANDIDATES.out.fasta                       // meta, combined.fa
    )
    ch_dedup = SAMTOOLS_SORMADUP.out.bam               // combined dedup BAM
    ch_versions = ch_versions.mix(SAMTOOLS_SORMADUP.out.versions.first())

    //
    // MODULE: Index the combined dedup BAM ONCE, then idxstats -> per-reference nodup counts (D-09).
    //         This same .bai is reused for the region split below (Open Q3 — two index passes total).
    //
    INDEX_DEDUP (
        ch_dedup
    )
    ch_versions = ch_versions.mix(INDEX_DEDUP.out.versions.first())

    IDXSTATS_NODUP (
        ch_dedup.join(INDEX_DEDUP.out.bai)             // meta, bam, bai
    )
    ch_versions = ch_versions.mix(IDXSTATS_NODUP.out.versions.first())

    //
    // Fan-out per candidate (D-11, Pattern 3): join the combined dedup BAM + its .bai with the
    // per-sample candidates CSV + the cand_fastas list, then flatMap one element per candidate
    // row. candidate_rank stays a STRING (never .toInteger() here — NA would crash, Pitfall 4),
    // id==sample is asserted, and the per-candidate reference FASTA is picked from cand_fastas by
    // the `_cand${rank}.` basename slot so IVAR_CONSENSUS / CONSENSUS_DISTANCE each get their own
    // reference. Only confirmation_status == 'pass' candidates are mapped (same gate as today).
    //
    ch_split = ch_dedup
        .join( INDEX_DEDUP.out.bai )                   // meta, bam, bai
        .join( ch_sample.fanout )                      // meta, bam, bai, candidates_csv, fastas
        .flatMap { meta, bam, bai, candidates_csv, fastas ->
            def rows = candidates_csv.splitCsv( header: true, sep: ',' )
            rows.collect { row ->
                def new_meta = meta + row              // carries candidate_rank (String), candidate_ref, confirmation_status
                assert new_meta.id == new_meta.sample : "Metadata mismatch: id=${new_meta.id}, sample=${new_meta.sample}"
                def rank  = new_meta.candidate_rank.toString()
                def fasta = fastas?.find { it.toString().contains("_cand${rank}.") }
                tuple(new_meta, bam, bai, fasta)
            }
        }
        // confirmation_status == 'pass' generalizes the legacy gate; the null-FASTA guard drops
        // below-threshold/absent candidates defensively (a passing candidate always has its FASTA).
        .filter { entry -> entry[0]['confirmation_status'] == 'pass' && entry[3] != null }

    // BAM-in / BAM-out region split. The positional region (candidate_ref) is supplied via
    // ext.args2 in Plan 02 config (lands after ${input} in the module command — NOT the -r
    // read-group flag, Pitfall 1). fasta input is [[],[]], qname is [], index_format is [].
    SAMTOOLS_VIEW (
        ch_split.map { meta, bam, bai, _fasta -> tuple(meta, bam, bai) }, // meta, bam, bai(index)
        [ [], [] ], // fasta (meta2, fasta) — not needed for BAM-in/BAM-out
        [],         // qname / readnames — not used
        []          // index_format — no on-the-fly index
    )
    ch_percand = SAMTOOLS_VIEW.out.bam                 // per-candidate dedup BAM
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())

    // Per-candidate reference FASTA channel, keyed by full meta (carries candidate_rank/ref),
    // so it joins back to the split BAM without cross-pairing same-sample candidates.
    ch_percand_fasta = ch_split.map { meta, _bam, _bai, fasta -> tuple(meta, fasta) }

    //
    // MODULE: Per-candidate depth — MUST run on the split (per-candidate) BAM, never the
    //         combined BAM (Pitfall 2: a combined depth file spans all references and corrupts
    //         breadth / cv_evenness).
    //
    SAMTOOLS_DEPTH (
        ch_percand,
        [ [], [] ] // Passing empty channels instead of an interval file
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    //
    // MODULE: Plot coverage from the per-candidate depth.
    //
    PLOTCOVERAGE (
        SAMTOOLS_DEPTH.out.tsv
    )
    ch_versions = ch_versions.mix(PLOTCOVERAGE.out.versions)

    //
    // MODULE: Plot variation in the per-candidate mapping file.
    //
    PLOT_BAMVARIATION (
        ch_percand
    )
    ch_versions = ch_versions.mix(PLOT_BAMVARIATION.out.versions)

    //
    // MODULE: Per-candidate consensus sequence. IVAR_CONSENSUS takes the reference FASTA as a
    //         bare positional path; supply this candidate's own reference FASTA.
    //
    IVAR_CONSENSUS (
        ch_percand,
        ch_percand_fasta.map { _meta, fasta -> fasta }, // bare path, positional
        false // Don't need the mpileup file
    )
    ch_versions = ch_versions.mix(IVAR_CONSENSUS.out.versions.first())

    //
    // MODULE: Compare per-candidate consensus to its own reference and compute distance.
    //
    CONSENSUS_DISTANCE (
        IVAR_CONSENSUS.out.fasta, // tuple val(meta), path(consensus.fa)
        ch_percand_fasta          // tuple val(meta), path(reference.fa)
    )
    ch_versions = ch_versions.mix(CONSENSUS_DISTANCE.out.versions.first())

    emit:
    aligned            = ch_percand                   // per-candidate dedup BAM
    idxstats_withdup   = IDXSTATS_WITHDUP.out.idxstats // combined BAM pre-dedup, per-reference counts
    idxstats_nodup     = IDXSTATS_NODUP.out.idxstats   // combined BAM post-dedup, per-reference counts
    depth              = SAMTOOLS_DEPTH.out.tsv
    consensus          = IVAR_CONSENSUS.out.fasta
    consensus_distance = CONSENSUS_DISTANCE.out.tsv
    variation          = PLOT_BAMVARIATION.out.png

    versions = ch_versions                            // channel: [ versions.yml ]
}
