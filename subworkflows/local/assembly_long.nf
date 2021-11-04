#!/usr/bin/env nextflow

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()
    
include {UNICYCLER as UNICYCLER_LONG} from '../../modules/local/unicycler'        addParams( options: modules['unicycler_long'])                          
include {FLYE} from '../../modules/local/flye'                                addParams( options: modules['flye']) 
include {MINIMAP2_ALIGN_LONG} from '../../modules/local/minimap_align_long'       addParams( options: modules['minimap_align_long']) 
include {MINIASM} from '../../modules/local/miniasm'                          addParams( options: modules['miniasm']) 
include {SEQSTATS as SEQSTATS_FLYE} from '../../modules/local/seqstats' addParams( tool: '_flye', options: modules['seqstats_assembly'])  
include {SEQSTATS as SEQSTATS_MINIASM} from '../../modules/local/seqstats' addParams( tool: '_miniasm', options: modules['seqstats_assembly'])      
include {SEQSTATS as SEQSTATS_UNICYCLER} from '../../modules/local/seqstats' addParams( tool: '_unicycler', options: modules['seqstats_assembly'])   

workflow RUN_ASSEMBLE_LONG {   

    take:
        long_reads
        short_reads
    main:
        ch_versions = Channel.empty()
        //Flye to be the best-performing bacterial genome assembler in many metrics
        if ( params.long_reads_assembler == 'flye'){
            FLYE ( long_reads )
            contigs = FLYE.out.assembly
            ch_versions = ch_versions.mix(FLYE.out.versions.first())
            SEQSTATS_FLYE(contigs)
            stats = SEQSTATS_FLYE.out.stats
        } 
        //failed testing with the current sample
        else if ( params.long_reads_assembler == 'miniasm' ){
            MINIMAP2_ALIGN_LONG(long_reads)
            MINIASM ( long_reads, MINIMAP2_ALIGN.out.paf)
            contigs = MINIASM.out.assembly
            ch_versions = ch_versions.mix(MINIASM.out.versions.first())
            SEQSTATS_MINIASM(contigs)
            stats = SEQSTATS_MINIASM.out.stats
            
        } 
        //If your long-read depth falls between those values, it might be worth trying both approaches.
        //If you have sparse long reads (~25× or less), flye is better

        else if (params.long_reads_assembler == 'unicycler' ) {
            UNICYCLER_LONG(long_reads)
            contigs = UNICYCLER_LONG.out.scaffolds
            ch_versions = ch_versions.mix(UNICYCLER_LONG.out.versions.first())
            SEQSTATS_UNICYCLER(contigs)
            stats = SEQSTATS_UNICYCLER.out.stats
            
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
        versions = ch_versions
        stats
}