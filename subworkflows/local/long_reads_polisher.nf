#!/usr/bin/env nextflow

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include {MINIMAP2_ALIGN_LONG as MAP1} from '../../modules/local/minimap_align_long'       addParams( options: modules['minimap_align_long']) 
include {MINIMAP2_ALIGN_LONG as MAP2} from '../../modules/local/minimap_align_long'       addParams( options: modules['minimap_align_long']) 
include {MINIMAP2_ALIGN_LONG as MAP3} from '../../modules/local/minimap_align_long'       addParams( options: modules['minimap_align_long']) 
include {MINIMAP2_ALIGN_LONG as MAP4} from '../../modules/local/minimap_align_long'       addParams( options: modules['minimap_align_long']) 

include {RACON as RACON1} from '../../modules/local/racon'                              addParams(racon_round: 1,  options: modules['racon']) 
include {RACON as RACON2} from '../../modules/local/racon'                              addParams(racon_round: 2,  options: modules['racon']) 
include {RACON as RACON3} from '../../modules/local/racon'                              addParams(racon_round: 3,  options: modules['racon']) 
include {RACON as RACON4} from '../../modules/local/racon'                              addParams(racon_round: 4,  options: modules['racon']) 

include {ASSEMBLY_STATS} from '../../modules/local/assembly_stats' addParams(options: modules['racon_assembly_stats'])     
workflow RUN_RACON_POLISH {   

    take:
        long_reads
        contigs
    main:
        ch_versions = Channel.empty()
        
        //1x
        MAP1(long_reads, contigs)
        ch_versions = ch_versions.mix(MAP1.out.versions.first())

        round1_paf = MAP1.out.paf
        RACON1(long_reads, round1_paf, contigs)
        ch_versions = ch_versions.mix(RACON1.out.versions.first())
        round1_asm = RACON1.out.assembly
        

        MAP2(long_reads, round1_asm)
        round2_paf = MAP2.out.paf
        RACON2(long_reads, round2_paf, round1_asm)
        round2_asm = RACON2.out.assembly
        

        MAP3(long_reads, round2_asm)
        round3_paf = MAP3.out.paf
        RACON3(long_reads, round3_paf, round2_asm)
        round3_asm = RACON3.out.assembly
        

        MAP4(long_reads, round3_asm)
        round4_paf = MAP4.out.paf
        RACON4(long_reads, round4_paf, round3_asm)
        assembly = RACON4.out.assembly
        ASSEMBLY_STATS(assembly)

    emit:
        assembly
        versions = ch_versions
        stats = ASSEMBLY_STATS.out.stats
        
}

