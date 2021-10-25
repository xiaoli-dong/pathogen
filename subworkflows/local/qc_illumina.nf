#!/usr/bin/env nextflow

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//include { SEQTK_FQCHK as ILLUMINASTATS_RAW} from '../modules/seqtk/fqchk/main'            addParams( options: modules['seqtk_fqchk_shortreads_raw'] )
//include {SEQTK_FQCHK as ILLUMINASTATS_QC} from '../modules/seqtk/fqchk/main'              addParams( options: modules['seqtk_fqchk_shortreads_qc'] )
include {FASTQC as FASTQC_RAW} from     '../../modules/nf-core/modules/fastqc/main'         addParams( options: modules['fastqc_raw'] )
include {FASTQC as FASTQC_QC} from      '../../modules/nf-core/modules/fastqc/main'         addParams( options: modules['fastqc_qc']) 
include { BBMAP_BBDUK as BBDUK } from '../../modules/local/bbmap_bbduk'                       addParams( options: modules['bbmap_bbduk']) 

workflow RUN_ILLUMINA_QC {   

    take:
        reads
    main:
        FASTQC_RAW(reads)
        //ILLUMINASTATS_RAW ( reads )
        //raw_stats = ILLUMINASTATS_RAW.out.seqstats
        BBDUK(reads, [])
        qc_reads = BBDUK.out.reads
        FASTQC_QC(qc_reads)
        //ILLUMINASTATS_QC (qc_reads)
        //qc_stats = ILLUMINASTATS_QC.out.seqstats
        versions = BBDUK.out.versions 

    emit:
        /* raw_stats
        qc_reads
        qc_stats
        version */
        qc_reads
        versions
        
}

