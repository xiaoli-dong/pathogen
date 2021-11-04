#!/usr/bin/env nextflow

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//include { SEQTK_FQCHK as ILLUMINASTATS_RAW} from '../modules/seqtk/fqchk/main'            addParams( options: modules['seqtk_fqchk_shortreads_raw'] )
//include {SEQTK_FQCHK as ILLUMINASTATS_QC} from '../modules/seqtk/fqchk/main'              addParams( options: modules['seqtk_fqchk_shortreads_qc'] )
include {FASTQC as FASTQC_RAW} from     '../../modules/nf-core/modules/fastqc/main'         addParams( options: modules['fastqc_raw'] )
include {FASTQC as FASTQC_QC} from      '../../modules/nf-core/modules/fastqc/main'         addParams( options: modules['fastqc_qc']) 
include { BBMAP_BBDUK as BBDUK } from '../../modules/local/bbmap_bbduk'                       addParams( options: modules['bbmap_bbduk']) 
include {SEQSTATS as SEQSTATS_RAW} from '../../modules/local/seqstats' addParams( tool: '', options: modules['seqstats_reads'])   
include {SEQSTATS as SEQSTATS_QC} from '../../modules/local/seqstats' addParams( tool: '', options: modules['seqstats_reads'])   

workflow RUN_ILLUMINA_QC {   

    take:
        reads
    main:
        ch_versions = Channel.empty()


        FASTQC_RAW(reads)
        ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())
        SEQSTATS_RAW(reads)
        raw_stats = SEQSTATS_RAW.out.stats
        //ILLUMINASTATS_RAW ( reads )
        //raw_stats = ILLUMINASTATS_RAW.out.seqstats
        BBDUK(reads, [])
        ch_versions = ch_versions.mix(BBDUK.out.versions.first())
        qc_reads = BBDUK.out.reads
        SEQSTATS_QC(qc_reads)
        qc_stats = SEQSTATS_QC.out.stats

        FASTQC_QC(qc_reads)
        //ILLUMINASTATS_QC (qc_reads)
        //qc_stats = ILLUMINASTATS_QC.out.seqstats
        

    emit:
        qc_reads
        versions = ch_versions
        raw_stats
        qc_stats
        
}

