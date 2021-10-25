#!/usr/bin/env nextflow

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//include { SHOVILL } from './../modules/shovill/main'                          addParams( options: modules['shovill']) 
include { SPADES } from '../../modules/nf-core/modules/spades/main'             addParams(spades_hmm: false , options: modules['spades']) 
include { SKESA } from '../../modules/local/skesa'                                addParams( options: modules['skesa']) 
include {UNICYCLER as UNICYCLER_SHORT} from '../../modules/local/unicycler'       addParams( options: modules['unicycler_short'])      

workflow RUN_ASSEMBLE_SHORT {   

    take:
        reads
    main:
        //[shovill] Assembly failed - spades.fasta has zero contigs!
        //disable it
       /*  if ( params.short_reads_assembler == 'shovill'){
            SHOVILL ( reads )   
            contigs = SHOVILL.out.contigs   
            version = SHOVILL.out.version 
        }  */
        if ( params.short_reads_assembler == 'spades' ){
            SPADES ( reads, [])
            contigs = SPADES.out.contigs
            versions = SPADES.out.versions
        } 
        //default
        else if (params.short_reads_assembler == 'skesa' ) {
            SKESA ( reads )
            contigs = SKESA.out.contigs
            versions = SKESA.out.versions
        }
        else if (params.short_reads_assembler == 'unicycler' ) {
            UNICYCLER_SHORT(reads)
            contigs = UNICYCLER_SHORT.out.scaffolds
            versions = UNICYCLER_SHORT.out.versions
        }
        
    emit:
        contigs
        versions
        
}
