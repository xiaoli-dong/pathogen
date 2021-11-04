#!/usr/bin/env nextflow

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()
                    
include {UNICYCLERHYBRID} from '../../modules/local/unicyclerhybrid'          addParams( options: modules['unicycler_hybrid']) 
include {SPADESHYBRID} from '../../modules/local/spadeshybrid'                addParams( options: modules['spades_hybrid']) 
include {SEQSTATS as SEQSTATS_SPADES} from '../../modules/local/seqstats' addParams( tool: '_spades_hybrid', options: modules['seqstats_assembly'])   
include {SEQSTATS as SEQSTATS_UNICYCLER} from '../../modules/local/seqstats' addParams( tool: '_unicycler_hybrid', options: modules['seqstats_assembly'])   

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
            SEQSTATS_SPADES(contigs)
            stats = SEQSTATS_SPADES.out.stats
            
        }

        //worst
        else if (params.hybrid_assembler  == 'spades'){
            SPADESHYBRID (long_reads, short_reads )   
            contigs = SPADESHYBRID.out.contigs    
            ch_versions = ch_versions.mix(SPADESHYBRID.out.versions.first())
            SEQSTATS_UNICYCLER(contigs)
            stats = SEQSTATS_UNICYCLER.out.stats
        } 

        
    emit:
        contigs
        versions = ch_versions
        stats
        
}
