#!/usr/bin/env nextflow

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//include { SHOVILL } from './../modules/shovill/main'                          addParams( options: modules['shovill']) 
include { SPADES } from '../../modules/nf-core/modules/spades/main'             addParams(spades_hmm: false , options: modules['spades']) 
include { SKESA } from '../../modules/local/skesa'                                addParams( options: modules['skesa']) 
include {UNICYCLER as UNICYCLER_SHORT} from '../../modules/local/unicycler'       addParams( options: modules['unicycler_short'])      
include {ASSEMBLY_STATS as STATS_SPADES} from '../../modules/local/assembly_stats' addParams(options: modules['spades_assembly_stats'])  
include {ASSEMBLY_STATS as STATS_SKESA} from '../../modules/local/assembly_stats' addParams(options: modules['skesa_assembly_stats'])      
include {ASSEMBLY_STATS as STATS_UNICYCLER} from '../../modules/local/assembly_stats' addParams(options: modules['unicycler_short_assembly_stats'])   
workflow RUN_ASSEMBLE_SHORT {   

    take:
        reads
    main:
        ch_versions = Channel.empty()

        if ( params.short_reads_assembler == 'spades' ){
            SPADES ( reads, [])
            contigs = SPADES.out.contigs
            ch_versions = ch_versions.mix(SPADES.out.versions.first())
            STATS_SPADES(contigs)
            stats = STATS_SPADES.out.stats
            ch_versions = ch_versions.mix(STATS_SPADES.out.versions.first())

        } 
        //default
        else if (params.short_reads_assembler == 'skesa' ) {
            SKESA ( reads )
            contigs = SKESA.out.contigs
            ch_versions = ch_versions.mix(SKESA.out.versions.first())
            STATS_SKESA(contigs)
            stats = STATS_SKESA.out.stats
            ch_versions = ch_versions.mix(STATS_SKESA.out.versions.first())
        }
        else if (params.short_reads_assembler == 'unicycler' ) {
            UNICYCLER_SHORT(reads)
            contigs = UNICYCLER_SHORT.out.scaffolds
            ch_versions = ch_versions.mix(UNICYCLER_SHORT.out.versions.first())
            STATS_UNICYCLER(contigs)
            stats = STATS_UNICYCLER.out.stats
            ch_versions = ch_versions.mix(STATS_UNICYCLER.out.versions.first())
        }
        
    emit:
        contigs
        stats
        versions = ch_versions
        
}
