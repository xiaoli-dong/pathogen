#!/usr/bin/env nextflow

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//nanopore
include {NANOPLOT as NANOPLOT_RAW} from '../../modules/nf-core/modules/nanoplot/main' addParams( options: modules['nanoplot']) 
include {NANOPLOT as NANOPLOT_QC}  from '../../modules/nf-core/modules/nanoplot/main' addParams( options: modules['nanoplot']) 
include {PORECHOP} from '../../modules/local//porechop'                              addParams( options: modules['porechop']) 
//include { SEQTK_FQCHK as NANOPORESTATS_RAW } from '../modules/seqtk/fqchk/main'    addParams( options: modules['seqtk_fqchk_longreads_raw'] )
//include { SEQTK_FQCHK as NANOPORESTATS_QC } from '../modules/seqtk/fqchk/main'     addParams( options: modules['seqtk_fqchk_longreads_qc'] )
include {SEQSTATS as SEQSTATS_RAW} from '../../modules/local/seqstats' addParams( tool: '', options: modules['seqstats_reads'])   
include {SEQSTATS as SEQSTATS_QC} from '../../modules/local/seqstats' addParams( tool: '', options: modules['seqstats_reads'])   


workflow RUN_NANOPORE_QC {

    take:
        reads
    main:

        ch_versions = Channel.empty()

        NANOPLOT_RAW(reads)
        ch_versions = ch_versions.mix(NANOPLOT_RAW.out.versions.first())
        SEQSTATS_RAW(reads)
        raw_stats = SEQSTATS_RAW.out.stats

        PORECHOP(reads)
        
        qc_reads = PORECHOP.out.reads
        ch_versions = ch_versions.mix(PORECHOP.out.versions.first())
        SEQSTATS_QC(reads)
        qc_stats = SEQSTATS_QC.out.stats


        NANOPLOT_QC(qc_reads)


    emit:
        raw_stats
        qc_reads
        qc_stats
        versions = ch_versions

}
