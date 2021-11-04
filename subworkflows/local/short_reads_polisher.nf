#!/usr/bin/env nextflow

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include {MINIMAP2_ALIGN_SHORT as MAPPING} from '../../modules/local/minimap_align_short'    addParams( options: modules['minimap_align_short']) 
include {SAMTOOLS_VIEW as VIEW} from '../../modules/nf-core/modules/samtools/view/main'     addParams( options: modules['samtools_view']) 
include {SAMTOOLS_SORT as SORT} from '../../modules/nf-core/modules/samtools/sort/main'     addParams( options: modules['samtools_sort']) 
include {SAMTOOLS_INDEX as INDEX} from '../../modules/nf-core/modules/samtools/index/main'  addParams( options: modules['samtools_index']) 
include {PILON } from '../../modules/local/pilon'                                           addParams(pilon_round: 1,  options: modules['pilon']) 
include {SEQSTATS as SEQSTATS_PILON} from '../../modules/local/seqstats' addParams( tool: '_pilon', options: modules['seqstats_assembly'])   

workflow RUN_PILON_POLISH {   

    take:
        short_reads
        contigs
    main:
        ch_versions = Channel.empty()

        MAPPING(short_reads, contigs)
        ch_versions = ch_versions.mix(MAPPING.out.versions.first())
        VIEW(MAPPING.out.sam)
        ch_versions = ch_versions.mix(VIEW.out.versions.first())

        SORT(VIEW.out.bam)
        INDEX(SORT.out.bam)
        PILON(SORT.out.bam, INDEX.out.bai, contigs)
        ch_versions = ch_versions.mix(PILON.out.versions.first())

        assembly = PILON.out.assembly 
        SEQSTATS_PILON(assembly)

    emit:
        assembly
        versions = ch_versions
        stats = SEQSTATS_PILON.out.stats

}

include {MINIMAP2_ALIGN_SHORT as MAPPING_STEP1} from '../../modules/local/minimap_align_short'      addParams( options: modules['minimap_align_short']) 

include {SAMTOOLS_VIEW as VIEW_STEP1} from '../../modules/nf-core/modules/samtools/view/main'       addParams( options: modules['samtools_view']) 
include {SAMTOOLS_SORT as SORT_STEP1} from '../../modules/nf-core/modules/samtools/sort/main'       addParams( options: modules['samtools_sort']) 
include {SAMTOOLS_INDEX as INDEX_STEP1} from '../../modules/nf-core/modules/samtools/index/main'    addParams( options: modules['samtools_index']) 
include {NEXTPOLISH as NEXTPOLISH1} from '../../modules/local/nextpolish'                           addParams(algorithm: 1,  options: modules['nextpolish']) 

include {MINIMAP2_ALIGN_SHORT as MAPPING_STEP2} from '../../modules/local/minimap_align_short'      addParams( options: modules['minimap_align_short']) 
include {SAMTOOLS_VIEW as VIEW_STEP2} from '../../modules/nf-core/modules/samtools/view/main'       addParams( options: modules['samtools_view']) 
include {SAMTOOLS_SORT as SORT_STEP2} from '../../modules/nf-core/modules/samtools/sort/main'       addParams( options: modules['samtools_sort']) 
include {SAMTOOLS_INDEX as INDEX_STEP2} from '../../modules/nf-core/modules/samtools/index/main'    addParams( options: modules['samtools_index']) 
include {NEXTPOLISH as NEXTPOLISH2} from '../../modules/local/nextpolish'                           addParams(algorithm: 2,  options: modules['nextpolish']) 
include {SEQSTATS as SEQSTATS_NEXTPOLISH} from '../../modules/local/seqstats' addParams( tool: '_nextpolish', options: modules['seqstats_assembly'])   

workflow RUN_NEXTPOLISH_POLISH {
    take:
        short_reads
        contigs
    main:
        ch_versions = Channel.empty()
        input = contigs
        // two iteration were recommended
        //for (int i = 0; i <2; i++) {
           // println("Hello World $i")

            MAPPING_STEP1(short_reads, input)
            ch_versions = ch_versions.mix(MAPPING_STEP1.out.versions.first())

            VIEW_STEP1(MAPPING_STEP1.out.sam)
            ch_versions = ch_versions.mix(VIEW_STEP1.out.versions.first())

            SORT_STEP1(VIEW_STEP1.out.bam)
            INDEX_STEP1(SORT_STEP1.out.bam)
            NEXTPOLISH1(SORT_STEP1.out.bam, INDEX_STEP1.out.bai, input)
            ch_versions = ch_versions.mix(NEXTPOLISH1.out.versions.first())

            input = NEXTPOLISH1.out.assembly 

            MAPPING_STEP2(short_reads, input)
            VIEW_STEP2(MAPPING_STEP2.out.sam)
            SORT_STEP2(VIEW_STEP2.out.bam)
            INDEX_STEP2(SORT_STEP2.out.bam)
            NEXTPOLISH2(SORT_STEP2.out.bam, INDEX_STEP2.out.bai, input)
            input = NEXTPOLISH2.out.assembly 
       //}
        assembly = input
        SEQSTATS_NEXTPOLISH(assembly)
    
    emit:
        assembly
        versions = ch_versions
        stats = SEQSTATS_NEXTPOLISH.out.stats
}

