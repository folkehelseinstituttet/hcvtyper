include { BOWTIE2_BUILD } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN } from '../../modules/nf-core/bowtie2/align/main' 
include { TANOTI_ALIGN  } from '../../modules/local/tanoti.nf'


workflow MAJOR_MAPPING {
    take:
    ch_major_mapping    //tuple val(meta), path(fasta), path(reads)
    prefix2

    main:
    // Extract the path(fasta)
    ch_build = ch_major_mapping.map { it[1] }
 
    // Extract val(meta), path(reads)
    ch_align = ch_major_mapping
        .map { meta, fasta, reads ->
        return [meta, reads]
        }
    
    // Extract only the name of the reference fasta mapped to. 
    // Splitting on '.' to remove the meta:id which is before the name.
    ref_name = ch_major_mapping.map { it[1].getBaseName().toString().split('\\.').last() }

    if (params.mapper == "bowtie2") {
            BOWTIE2_BUILD (
            ch_build
        )

        BOWTIE2_ALIGN (
            ch_align,
            BOWTIE2_BUILD.out.index,
            false,
            true,
            ref_name,
            prefix2
        )
        aligned = BOWTIE2_ALIGN.out.aligned
        depth = BOWTIE2_ALIGN.out.depth
        stats_markdup = BOWTIE2_ALIGN.out.stats_markdup
        stats_withdup = BOWTIE2_ALIGN.out.stats_withdup
        versions = BOWTIE2_ALIGN.out.versions // channel: [ versions.yml ]
    } 
    else if (params.mapper == "tanoti") {
        TANOTI_ALIGN (
            ch_align,
            ch_build,
            false,
            true,
            ref_name,
            prefix2,
            params.tanoti_stringency_2
        )
        aligned = TANOTI_ALIGN.out.aligned
        depth = TANOTI_ALIGN.out.depth
        stats_markdup = TANOTI_ALIGN.out.stats_markdup
        stats_withdup = TANOTI_ALIGN.out.stats_withdup
        versions = TANOTI_ALIGN.out.versions // channel: [ versions.yml ]
    }

    emit:
    aligned
    versions
    depth
    stats_markdup
    stats_withdup
    
}
