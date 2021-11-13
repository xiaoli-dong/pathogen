#!/usr/bin/env nextflow

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()
                    
include {UNICYCLERHYBRID} from '../../modules/local/unicyclerhybrid'          addParams( options: modules['unicycler_hybrid']) 
include {SPADESHYBRID} from '../../modules/local/spadeshybrid'                addParams( options: modules['spades_hybrid']) 
include {ASSEMBLY_STATS as STATS_SPADES} from '../../modules/local/assembly_stats' addParams(options: modules['spades_hybrid_assembly_stats'])   
include {ASSEMBLY_STATS as STATS_UNICYCLER} from '../../modules/local/assembly_stats' addParams(options: modules['unicycler_hybrid_assembly_stats'])   

workflow RUN_ASSEMBLE_HYBRID {   

    take:
        long_reads
        short_reads
    main:
        ch_versions = Channel.empty()
        //works the best, default can be set to unicycler
        //Unicycler works best when the short-read set is very good (deep and complete coverage) 
        //which yields a nice short-read assembly graph for scaffolding.
        if (params.hybrid_assembler  == 'unicycler' ) {
            UNICYCLERHYBRID(long_reads, short_reads)
           
            contigs = UNICYCLERHYBRID.out.scaffolds
            ch_versions = ch_versions.mix(UNICYCLERHYBRID.out.versions.first())
            STATS_SPADES(contigs)
            stats = STATS_SPADES.out.stats
            
        }

        //worst
        else if (params.hybrid_assembler  == 'spades'){
            SPADESHYBRID (long_reads, short_reads )   
            contigs = SPADESHYBRID.out.contigs    
            ch_versions = ch_versions.mix(SPADESHYBRID.out.versions.first())
            STATS_UNICYCLER(contigs)
            stats = STATS_UNICYCLER.out.stats
        } 

        
    emit:
        contigs
        versions = ch_versions
        stats
        
}
