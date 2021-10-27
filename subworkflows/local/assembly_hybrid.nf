#!/usr/bin/env nextflow

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()
                    
include {UNICYCLERHYBRID} from '../../modules/local/unicyclerhybrid'          addParams( options: modules['unicycler_hybrid']) 
include {SPADESHYBRID} from '../../modules/local/spadeshybrid'                addParams( options: modules['spades_hybrid']) 

workflow RUN_ASSEMBLE_HYBRID {   

    take:
        long_reads
        short_reads
    main:
        
        //works the best, default can be set to unicycler
        //Unicycler works best when the short-read set is very good (deep and complete coverage) 
        //which yields a nice short-read assembly graph for scaffolding.
        if (params.hybrid_assembler  == 'unicycler' ) {
            UNICYCLERHYBRID(long_reads, short_reads)
           
            contigs = UNICYCLERHYBRID.out.scaffolds
            versions = UNICYCLERHYBRID.out.versions
            
        }

        //worst
        else if (params.hybrid_assembler  == 'spades'){
            SPADESHYBRID (long_reads, short_reads )   
            contigs = SPADESHYBRID.out.contigs    
            versions = SPADESHYBRID.out.versions
        } 

        
    emit:
        contigs
        versions
        
}
