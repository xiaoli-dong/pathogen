#!/usr/bin/env nextflow

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()


include {UNICYCLER as UNICYCLER_LONG} from '../../modules/nf-core/modules/unicycler/main'        addParams( options: modules['unicycler_long'])                          
include {UNICYCLERHYBRID} from '../../modules/local/unicyclerhybrid'          addParams( options: modules['unicycler_hybrid']) 
include {SPADESHYBRID} from '../../modules/local/spadeshybrid'                addParams( options: modules['spades_hybrid']) 
include {FLYE} from '../../modules/local/flye'                                addParams( options: modules['flye']) 
include {MINIMAP2_ALIGN} from '../../modules/nf-core/modules/minimap2/align/main'       addParams( options: modules['minimap_align']) 
include {MINIASM} from '../../modules/local/miniasm'                          addParams( options: modules['miniasm']) 
include {RACON} from '../../modules/local/racon'                              addParams( options: modules['racon']) 
include {MEDAKA} from '../../modules/local/medaka'                            addParams( options: modules['medaka']) 
//illumina polisher
include {POLCA} from '../../modules/local/polca'                              addParams( options: modules['polca']) 
include {NEXTPOLISH} from '../../modules/local/nextpolish'                    addParams( options: modules['nextpolish']) 
include {PILON} from '../../modules/local/pilon'                              addParams( options: modules['pilon'])  


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
        else if ( hybrid_assembler  == 'spades'){
            SPADESHYBRID (long_reads, short_reads )   
            contigs = SPADESHYBRID.out.contigs    
            versions = SPADESHYBRID.out.versions
        } 

        
    emit:
        contigs
        versions
        
}
