#!/usr/bin/env nextflow

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//nanopore
include {NANOPLOT as NANOPLOT_RAW} from '../../modules/nf-core/modules/nanoplot/main' addParams( options: modules['nanoplot']) 
include {NANOPLOT as NANOPLOT_QC}  from '../../modules/nf-core/modules/nanoplot/main' addParams( options: modules['nanoplot']) 
include {PORECHOP} from '../../modules/local//porechop'                              addParams( options: modules['porechop']) 
//include { SEQTK_FQCHK as NANOPORESTATS_RAW } from '../modules/seqtk/fqchk/main'    addParams( options: modules['seqtk_fqchk_longreads_raw'] )
//include { SEQTK_FQCHK as NANOPORESTATS_QC } from '../modules/seqtk/fqchk/main'     addParams( options: modules['seqtk_fqchk_longreads_qc'] )


workflow RUN_NANOPORE_QC {

    take:
        reads
    main:
        NANOPLOT_RAW(reads)

        //NANOPORESTATS_RAW(reads)
        //raw_stats = NANOPORESTATS_RAW.out.seqstats

        PORECHOP(reads)
        qc_reads = PORECHOP.out.reads

        NANOPLOT_QC(qc_reads)

        //NANOPORESTATS_QC (qc_reads)
        //qc_stats = NANOPORESTATS_QC.out.seqstats
        versions = PORECHOP.out.versions

    emit:
        //raw_stats
        qc_reads
        //qc_stats
        versions

}
