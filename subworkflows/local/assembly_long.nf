#!/usr/bin/env nextflow

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()
    
include {UNICYCLER as UNICYCLER_LONG} from '../../modules/local/unicycler'        addParams( options: modules['unicycler_long'])                          
include {FLYE} from '../../modules/local/flye'                                addParams( options: modules['flye']) 
include {MINIMAP2_ALIGN_LONG} from '../../modules/local/minimap_align_long'       addParams( options: modules['minimap_align_long']) 
include {MINIASM} from '../../modules/local/miniasm'                          addParams( options: modules['miniasm']) 


workflow RUN_ASSEMBLE_LONG {   

    take:
        long_reads
        short_reads
    main:
        //Flye to be the best-performing bacterial genome assembler in many metrics
        if ( params.long_reads_assembler == 'flye'){
            FLYE ( long_reads )
            contigs = FLYE.out.assembly
            versions = FLYE.out.versions
        } 
        //failed testing with the current sample
        else if ( params.long_reads_assembler == 'miniasm' ){
            MINIMAP2_ALIGN_LONG(long_reads)
            MINIASM ( long_reads, MINIMAP2_ALIGN.out.paf)
            contigs = MINIASM.out.assembly
            versions = MINIASM.out.versions   
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