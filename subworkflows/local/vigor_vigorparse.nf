include { VIGOR               } from '../../modules/local/vigor4'
include { VIGOR_GFF_EXTRACT   } from '../../modules/local/vigorparse/ExtractFromGff'
include { VIGOR_HIGH_COVERAGE } from '../../modules/local/vigorparse/FindHighCoverage'

workflow VIGOR_VIGORPARSE {
    take:
    ch_vigor    //tuple val(meta), path(contigs)


    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Run Vigor4 on Spades contigs
    //
    VIGOR(
        ch_vigor
    )
    ch_versions = ch_versions.mix(VIGOR.out.versions.first())

    //
    // MODULE:Extract Rotavirus segments from Vigor output
    //
    VIGOR_GFF_EXTRACT(
        VIGOR.out.gff3.join(VIGOR.out.contigs) // Create a tuple channel with meta, gff3 and contigs
    )

    //
    // MODULE: Find contig with highest coverage per segment
    //
    VIGOR_HIGH_COVERAGE(
        // Create tuple channels with with meta and fasta for each fasta file from gff extract
        VIGOR_GFF_EXTRACT.out.gene_fasta.transpose()
    )
    ch_versions = ch_versions.mix(VIGOR.out.versions.first())

    emit:
    gene_fasta = VIGOR_HIGH_COVERAGE.out.gene_fasta
    versions = ch_versions

}
