include { BOWTIE2_BUILD } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN } from '../../modules/nf-core/bowtie2/align/main' 

workflow MAJOR_MAPPING {
    take:
    ch_build
    ch_align
    prefix2

    main:
    if (params.mapper == "bowtie2") {
            BOWTIE2_BUILD (
            ch_build
        )

        BOWTIE2_ALIGN (
            ch_align,
            BOWTIE2_BUILD.out.index,
            false,
            true,
            prefix2
        )
        aligned = BOWTIE2_ALIGN.out.aligned
        depth = BOWTIE2_ALIGN.out.depth
        versions = BOWTIE2_ALIGN.out.versions // channel: [ versions.yml ]
    } 
    else if (params.mapper == "tanoti") {
        TANOTI (
            ch_align,
            false,
            true,
            prefix2
        )
        aligned = TANOTI.out.aligned
        depth = TANOTI.out.depth
        versions = TANOTI.out.versions // channel: [ versions.yml ]
    }

    emit:
    aligned
    versions
    
}