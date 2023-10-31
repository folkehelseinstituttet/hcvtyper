include { BOWTIE2_BUILD as BOWTIE2_BUILD_MAJOR } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN as BOWTIE2_MAJOR       } from '../../modules/nf-core/bowtie2/align/main' 

workflow MAJOR_MAPPING {
    take:
    ch_build_major
    ch_map_major

    main:
    if (params.mapper == "bowtie2") {
            BOWTIE2_BUILD_MAJOR (
            ch_build_major
        )

        BOWTIE2_MAJOR (
            ch_map_major,
            BOWTIE2_BUILD_MAJOR.out.index,
            false,
            true
        )
    } 
    else if (params.mapper == "tanoti") {
     TANOTI_MAJOR (
       ch_map_major
     )
    }

    
}