#!/usr/bin/env nextflow

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//include { SHOVILL } from './../modules/shovill/main'                          addParams( options: modules['shovill']) 
include { SPADES } from '../../modules/nf-core/modules/spades/main'                              addParams( options: modules['spades']) 
include { SKESA } from '../../modules/local/skesa'                                addParams( options: modules['skesa']) 
include {UNICYCLER as UNICYCLER_SHORT} from '../../modules/nf-core/modules/unicycler/main'       addParams( options: modules['unicycler_short'])      
/* include {UNICYCLER as UNICYCLER_LONG} from '../../modules/nf-core/modules/unicycler/main'        addParams( options: modules['unicycler_long'])                          
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
include {PILON} from '../../modules/local/pilon'                              addParams( options: modules['pilon'])  */

workflow RUN_ASSEMBLE_SHORT {   

    take:
        reads
    main:
        //[shovill] Assembly failed - spades.fasta has zero contigs!
        //disable it
       /*  if ( params.assembler == 'shovill'){
            SHOVILL ( reads )   
            contigs = SHOVILL.out.contigs   
            version = SHOVILL.out.version 
        }  */
        if ( params.assembler == 'spades' ){
            SPADES ( reads)
            contigs = SPADES.out.contigs
            version = SPADES.out.version
        } 
        //default
        else if (params.assembler == 'skesa' ) {
            SKESA ( reads )
            contigs = SKESA.out.contigs
            version = SKESA.out.version
        }
        else if (params.assembler == 'unicycler' ) {
            UNICYCLER_SHORT(reads)
            contigs = UNICYCLER_SHORT.out.scaffolds
            version = UNICYCLER_SHORT.out.version
        }
        
    emit:
        contigs
        version
        
}

workflow RUN_ASSEMBLE_LONG {   

    take:
        long_reads
        short_reads
    main:
        //Flye to be the best-performing bacterial genome assembler in many metrics
        if ( params.assembler == 'flye'){
            FLYE ( long_reads )

            //4x long reads polish
            RACON(long_reads, FLYE.out.assembly)
            //1x consensus
            MEDAKA(long_reads, RACON.out.racon_assembly)

            if(params.illumina_polisher == "pilon"){
                PILON(short_reads, MEDAKA.out.assembly)
                contigs = PILON.out.assembly    
                version = PILON.out.version
            }
            else{
                // 1x polca
                POLCA(short_reads, MEDAKA.out.assembly)
                //2x nextpolish
                NEXTPOLISH(short_reads, POLCA.out.assembly)

                contigs = NEXTPOLISH.out.assembly    
                version = NEXTPOLISH.out.version
            }
        
            
        } 
        //failed testing with the current sample
        else if ( params.assembler == 'miniasm' ){
            MINIMAP2_ALIGN(long_reads)
            MINIASM ( long_reads, MINIMAP2_ALIGN.out.paf)
            RACON(long_reads, MINIASM.out.assembly)
            MEDAKA(long_reads, RACON.out.racon_assembly)
            
            if(params.illumina_polisher == "pilon"){
                PILON(short_reads, MEDAKA.out.assembly)
                contigs = PILON.out.assembly    
                version = PILON.out.version
            }
            else if (params.illumina_polisher == "polca+nextpolish"){
                // 1x polca
                POLCA(short_reads, MEDAKA.out.assembly)
                //2x nextpolish
                NEXTPOLISH(short_reads, POLCA.out.assembly)

                contigs = NEXTPOLISH.out.assembly    
                version = NEXTPOLISH.out.version
            }

            
        } 
        //If your long-read depth falls between those values, it might be worth trying both approaches.
        //If you have sparse long reads (~25× or less), flye is better

        else if (params.assembler == 'unicycler' ) {
            UNICYCLER_LONG(long_reads)
            contigs = UNICYCLER_LONG.out.scaffolds
            version = UNICYCLER_LONG.out.version   
        }
        // If you have lots of long reads (~100× depth or more), use Trycycler+polishing.
        // it need manual intervation, so skip it
        /* else if (params.assembler == 'trycycler' ) {
            UNICYCLER(reads)
            contigs = UNICYCLER.out.scaffolds
            version = UNICYCLER.out.version   
        } */
       
    emit:
        contigs
        version
}

workflow RUN_ASSEMBLE_HYBRID {   

    take:
        long_reads
        short_reads
    main:
        
        //works the best, default can be set to unicycler
        //Unicycler works best when the short-read set is very good (deep and complete coverage) 
        //which yields a nice short-read assembly graph for scaffolding.
        if (params.assembler == 'unicycler' ) {
            UNICYCLERHYBRID(long_reads, short_reads)
           
            contigs = UNICYCLERHYBRID.out.scaffolds
            version = UNICYCLERHYBRID.out.version
            
        }

        //worst
        else if ( params.assembler == 'spades'){
            SPADESHYBRID (long_reads, short_reads )   
            contigs = SPADESHYBRID.out.contigs    
            version = SPADESHYBRID.out.version 
        } 

        
    emit:
        contigs
        version
        
}
