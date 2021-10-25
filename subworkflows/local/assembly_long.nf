#!/usr/bin/env nextflow

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()
    
include {UNICYCLER as UNICYCLER_LONG} from '../../modules/local/unicycler'        addParams( options: modules['unicycler_long'])                          
include {FLYE} from '../../modules/local/flye'                                addParams( options: modules['flye']) 
include {MINIMAP2_ALIGN} from '../../modules/nf-core/modules/minimap2/align/main'       addParams( options: modules['minimap_align']) 
include {MINIASM} from '../../modules/local/miniasm'                          addParams( options: modules['miniasm']) 
include {RACON} from '../../modules/local/racon'                              addParams( options: modules['racon']) 
include {MEDAKA} from '../../modules/local/medaka'                            addParams( options: modules['medaka']) 
//illumina polisher
include {POLCA} from '../../modules/local/polca'                              addParams( options: modules['polca']) 
include {NEXTPOLISH} from '../../modules/local/nextpolish'                    addParams( options: modules['nextpolish']) 
include {PILON} from '../../modules/local/pilon'                              addParams( options: modules['pilon']) 

workflow RUN_ASSEMBLE_LONG {   

    take:
        long_reads
        short_reads
    main:
        //Flye to be the best-performing bacterial genome assembler in many metrics
        if ( params.long_reads_assembler == 'flye'){
            FLYE ( long_reads )

            //4x long reads polish
            RACON(long_reads, FLYE.out.assembly)
            //1x consensus
            MEDAKA(long_reads, RACON.out.racon_assembly)

            if(params.short_reads_polisher  == "pilon"){
                PILON(short_reads, MEDAKA.out.assembly)
                contigs = PILON.out.assembly    
                versions = PILON.out.versions
            }
            else if(params.short_reads_polisher  == "1xpolca+2xnextpolish"){
                // 1x polca
                POLCA(short_reads, MEDAKA.out.assembly)
                //2x nextpolish
                NEXTPOLISH(short_reads, POLCA.out.assembly)

                contigs = NEXTPOLISH.out.assembly    
                versions = NEXTPOLISH.out.versions
            }
            else{
                contigs = MEDAKA.out.assembly
                versions = FLYE.out.versions
            }
        
            
        } 
        //failed testing with the current sample
        else if ( params.long_reads_assembler == 'miniasm' ){
            MINIMAP2_ALIGN(long_reads)
            MINIASM ( long_reads, MINIMAP2_ALIGN.out.paf)
            RACON(long_reads, MINIASM.out.assembly)
            MEDAKA(long_reads, RACON.out.racon_assembly)
            
           if(params.short_reads_polisher  == "pilon"){
                PILON(short_reads, MEDAKA.out.assembly)
                contigs = PILON.out.assembly    
                versions = PILON.out.versions
            }
            else if(params.short_reads_polisher  == "1xpolca+2xnextpolish"){
                // 1x polca
                POLCA(short_reads, MEDAKA.out.assembly)
                //2x nextpolish
                NEXTPOLISH(short_reads, POLCA.out.assembly)

                contigs = NEXTPOLISH.out.assembly    
                versions = NEXTPOLISH.out.versions
            }
            else{
                contigs = MEDAKA.out.assembly
                versions = FLYE.out.versions
            }
        

            
        } 
        //If your long-read depth falls between those values, it might be worth trying both approaches.
        //If you have sparse long reads (~25× or less), flye is better

        else if (params.long_reads_assembler == 'unicycler' ) {
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
        versions
}