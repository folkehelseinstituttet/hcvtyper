//
// Detect potential cross-sample contig contamination from de novo assembly.
//
// For each sample: filter contigs >= params.contamination_min_length bp and
// embed the sample ID in FASTA headers. All filtered contigs are then combined
// and run through a self-vs-self BLAST to detect near-identical sequences
// shared between samples.
//

include { CONTIG_FILTER_RENAME                          } from '../../../modules/local/contig_filter_rename/main'
include { CAT_FILTERED_CONTIGS                          } from '../../../modules/local/cat_filtered_contigs/main'
include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB_CONTIGS } from '../../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN      as BLAST_BLASTN_CONTIGS      } from '../../../modules/nf-core/blast/blastn/main'
include { CONTAMINATION_REPORT                          } from '../../../modules/local/contamination_report/main'

workflow CONTAMINATION_CHECK {

    take:
    ch_contigs // channel: [ val(meta), path(contigs) ]  — gzipped FASTA from SPADES

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Filter contigs >= min_length bp and prefix headers with sample ID
    //
    CONTIG_FILTER_RENAME(ch_contigs)
    ch_versions = ch_versions.mix(CONTIG_FILTER_RENAME.out.versions.first())

    // Discard samples that produced no long contigs (empty output files)
    // then collect all filtered FASTAs into a single list for concatenation.
    // This step collects all samples so needs to wait until SPADES is completed for all samples.
    ch_all_fas = CONTIG_FILTER_RENAME.out.filtered
        .filter { _meta, fa -> fa.size() > 0 }
        .map    { _meta, fa -> fa }
        .collect()

    //
    // MODULE: Concatenate all per-sample filtered FASTAs into one combined file
    //
    CAT_FILTERED_CONTIGS(ch_all_fas)
    ch_versions = ch_versions.mix(CAT_FILTERED_CONTIGS.out.versions)

    // Add a meta map so the nf-core BLAST modules accept the combined FASTA
    ch_combined = CAT_FILTERED_CONTIGS.out.fasta
        .map { fasta -> [ [ id: "all_samples_contigs" ], fasta ] }

    //
    // MODULE: Build a nucleotide BLAST database from the combined FASTA
    //
    BLAST_MAKEBLASTDB_CONTIGS(ch_combined)
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB_CONTIGS.out.versions)

    //
    // MODULE: All-vs-all BLAST — query combined FASTA against itself
    //
    BLAST_BLASTN_CONTIGS(
        ch_combined,
        BLAST_MAKEBLASTDB_CONTIGS.out.db,
        [],     // no taxidlist
        "",     // no inline taxids
        false   // no negative filtering
    )
    ch_versions = ch_versions.mix(BLAST_BLASTN_CONTIGS.out.versions)

    //
    // MODULE: Parse BLAST output — produce TSV, heatmap, and MultiQC JSON
    //
    CONTAMINATION_REPORT(
        BLAST_BLASTN_CONTIGS.out.txt.map { _meta, txt -> txt }
    )
    ch_versions = ch_versions.mix(CONTAMINATION_REPORT.out.versions)

    emit:
    tsv      = CONTAMINATION_REPORT.out.tsv      // path: contamination_pairs.tsv
    heatmap  = CONTAMINATION_REPORT.out.heatmap  // path: contamination_heatmap.png
    mqc      = CONTAMINATION_REPORT.out.mqc      // path: contamination_mqc.json
    versions = ch_versions                       // channel: versions.yml
}
