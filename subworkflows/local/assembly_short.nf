#!/usr/bin/env nextflow

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//include { SHOVILL } from './../modules/shovill/main'                          addParams( options: modules['shovill']) 
include { SPADES } from '../../modules/nf-core/modules/spades/main'             addParams(spades_hmm: false , options: modules['spades']) 
include { SKESA } from '../../modules/local/skesa'                                addParams( options: modules['skesa']) 
include {UNICYCLER as UNICYCLER_SHORT} from '../../modules/local/unicycler'       addParams( options: modules['unicycler_short'])      
include {SEQSTATS as SEQSTATS_SPADES} from '../../modules/local/seqstats' addParams( tool: '_spades', options: modules['seqstats_assembly'])  
include {SEQSTATS as SEQSTATS_SKESA} from '../../modules/local/seqstats' addParams( tool: '_skesa', options: modules['seqstats_assembly'])      
include {SEQSTATS as SEQSTATS_UNICYCLER} from '../../modules/local/seqstats' addParams( tool: '_unicycler', options: modules['seqstats_assembly'])   
workflow RUN_ASSEMBLE_SHORT {   

    take:
        reads
    main:
        ch_versions = Channel.empty()

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
            ch_versions = ch_versions.mix(SPADES.out.versions.first())
            SEQSTATS_SPADES(contigs)
            stats = SEQSTATS_SPADES.out.stats
            ch_versions = ch_versions.mix(SEQSTATS.out.versions.first())

        } 
        //default
        else if (params.short_reads_assembler == 'skesa' ) {
            SKESA ( reads )
            contigs = SKESA.out.contigs
            ch_versions = ch_versions.mix(SKESA.out.versions.first())
            SEQSTATS_SKESA(contigs)
            stats = SEQSTATS_SKESA.out.stats
            ch_versions = ch_versions.mix(SEQSTATS.out.versions.first())
        }
        else if (params.short_reads_assembler == 'unicycler' ) {
            UNICYCLER_SHORT(reads)
            contigs = UNICYCLER_SHORT.out.scaffolds
            ch_versions = ch_versions.mix(UNICYCLER_SHORT.out.versions.first())
            SEQSTATS_UNICYCLER(contigs)
            stats = SEQSTATS_UNICYCLER.out.stats
            ch_versions = ch_versions.mix(SEQSTATS.out.versions.first())
        }
        
    emit:
        contigs
        stats
        versions = ch_versions
        
}
