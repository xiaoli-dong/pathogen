#!/usr/bin/env nextflow

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include { SEQTK_FQCHK as ILLUMINASTATS_RAW} from '../modules/seqtk/fqchk/main'      addParams( options: modules['seqtk_fqchk_shortreads_raw'] )
include {SEQTK_FQCHK as ILLUMINASTATS_QC} from '../modules/seqtk/fqchk/main'        addParams( options: modules['seqtk_fqchk_shortreads_qc'] )
include {FASTQC as FASTQC_RAW} from  '../modules/fastqc/main'                       addParams( options: modules['fastqc']) 
include {FASTQC as FASTQC_QC} from  '../modules/fastqc/main'                        addParams( options: modules['fastqc']) 
include { BBMAP_BBDUK as BBDUK } from '../modules/bbmap/bbduk/main'                 addParams( options: modules['bbmap_bbduk']) 

//nanopore
include {NANOPLOT as NANOPLOT_RAW} from '../modules/nanopore/nanoplot'              addParams( options: modules['nanoplot']) 
include {NANOPLOT as NANOPLOT_QC} from '../modules/nanopore/nanoplot'               addParams( options: modules['nanoplot']) 
include {PORECHOP} from '../modules/nanopore/porechop'                              addParams( options: modules['porechop']) 
include { SEQTK_FQCHK as NANOPORESTATS_RAW } from '../modules/seqtk/fqchk/main'    addParams( options: modules['seqtk_fqchk_longreads_raw'] )
include { SEQTK_FQCHK as NANOPORESTATS_QC } from '../modules/seqtk/fqchk/main'     addParams( options: modules['seqtk_fqchk_longreads_qc'] )

workflow RUN_ILLUMINA_QC {   

    take:
        reads
    main:
        FASTQC_RAW(reads)
        ILLUMINASTATS_RAW ( reads )
        raw_stats = ILLUMINASTATS_RAW.out.seqstats
        BBDUK(reads, [])
        qc_reads = BBDUK.out.reads
        FASTQC_QC(qc_reads)
        ILLUMINASTATS_QC (qc_reads)
        qc_stats = ILLUMINASTATS_QC.out.seqstats
        version = BBDUK.out.version

    emit:
        raw_stats
        qc_reads
        qc_stats
        version
        
}


workflow RUN_NANOPORE_QC {

    take:
        reads
    main:
        NANOPLOT_RAW(reads)

        NANOPORESTATS_RAW(reads)
        raw_stats = NANOPORESTATS_RAW.out.seqstats

        PORECHOP(reads)
        qc_reads = PORECHOP.out.reads

        NANOPLOT_QC(qc_reads)

        NANOPORESTATS_QC (qc_reads)
        qc_stats = NANOPORESTATS_QC.out.seqstats
        version = PORECHOP.out.version

    emit:
        raw_stats
        qc_reads
        qc_stats
        version 

}
