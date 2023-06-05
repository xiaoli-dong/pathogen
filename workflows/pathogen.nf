
// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )
include {MEDAKA} from '../modules/local/medaka' addParams( options: modules['medaka'])
include {POLCA} from '../modules/local/polca'                              addParams( options: modules['polca']) 
include {BRACKEN} from '../modules/local/bracken'        addParams( options: modules['kraken2_bracken'])   
include { BAKTA } from '../modules/local/bakta'                addParams( options: modules['bakta'])
include {GFF2FEATURES as BAKTA_FEATURES} from '../modules/local/gff2features'                addParams( options: modules['bakta_features'])
include {GFF2FEATURES as PROKKA_FEATURES} from '../modules/local/gff2features'                addParams( options: modules['prokka_features'])
include { MOBSUITE } from '../modules/local/mobsuite'          addParams( options: modules['mobsuite'])
include { ABRICATE as ABRICATE_VF} from '../modules/local/abricate' addParams( options: modules['abricate_vf'])
include { ABRICATE_SUMMARIZE as ABRICATE_SUMMARIZE_VF} from '../modules/local/abricate' addParams( options: modules['abricate_vf_summarize'])
include {AMRFINDERPLUS} from '../modules/local/amrfinderplus' addParams( options: modules['amrfinderplus'])


include {GET_SAMPLEIDS} from '../modules/local/get_sampleids' addParams( options: modules['csvtk_concat'])
include {
    CSVTK_CONCAT as CONCAT_STATS_SHORT_RAW;
    CSVTK_CONCAT as CONCAT_STATS_SHORT_QC;
    CSVTK_CONCAT as CONCAT_STATS_LONG_RAW;
    CSVTK_CONCAT as CONCAT_STATS_LONG_QC;
    CSVTK_CONCAT as CONCAT_STATS_ASM;
    CSVTK_CONCAT as CONCAT_STATS_BAKTA;
    CSVTK_CONCAT as CONCAT_STATS_PROKKA;
    CSVTK_CONCAT as CONCAT_MOBSUITE;
    CSVTK_CONCAT as CONCAT_AMR;

} from '../modules/local/csvtk_concat' addParams( header: '',  options: modules['csvtk_concat'])


include {
    CSVTK_CONCAT as CONCAT_MLST;

} from '../modules/local/csvtk_concat' addParams(header: '-H', options: modules['csvtk_concat'])

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/
//
// MODULE: Installed directly from nf-core/modules
//
//include { FASTQC  } from '../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
include {KRAKEN2_KRAKEN2} from '../modules/nf-core/modules/kraken2/kraken2/main'        addParams( options: modules['kraken2'])   
include { MLST } from '../modules/nf-core/modules/mlst/main'                  addParams( options: modules['mlst'])
include { PROKKA } from '../modules/nf-core/modules/prokka/main'              addParams( options: modules['prokka'])
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main' addParams( options: [publish_files : ['_versions.yml':'']] )



//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )
include {RUN_ILLUMINA_QC} from '../subworkflows/local/qc_illumina'
include {RUN_NANOPORE_QC} from '../subworkflows/local/qc_nanopore'
include {RUN_ASSEMBLE_SHORT} from '../subworkflows/local/assembly_short'
include {RUN_ASSEMBLE_LONG} from '../subworkflows/local/assembly_long'
include {RUN_ASSEMBLE_HYBRID} from '../subworkflows/local/assembly_hybrid'
include {
    RUN_RACON_POLISH;
} from '../subworkflows/local/long_reads_polisher'
include {
    RUN_PILON_POLISH;
    RUN_PILON_POLISH as RUN_PILON_POLISH2;
    RUN_PILON_POLISH as RUN_PILON_POLISH3;
    RUN_PILON_POLISH as RUN_PILON_POLISH4;
    RUN_NEXTPOLISH_POLISH;
    RUN_NEXTPOLISH_POLISH as RUN_NEXTPOLISH_POLISH2;
} from '../subworkflows/local/short_reads_polisher'

include {
    ARG
} from '../subworkflows/local/arg'

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowPathogen.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
//def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
def checkPathParamList = [ params.input]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

 
/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow PATHOGEN {

    ch_software_versions = Channel.empty()
    
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    
    reads = INPUT_CHECK.out.reads
    short_reads = INPUT_CHECK.out.shortreads
    long_reads = INPUT_CHECK.out.longreads

    
    INPUT_CHECK.out.ids.collect() | GET_SAMPLEIDS
    


    if(!params.skip_short_reads_qc){
        RUN_ILLUMINA_QC(short_reads)
        ch_software_versions = ch_software_versions.mix(RUN_ILLUMINA_QC.out.versions)
        short_reads = RUN_ILLUMINA_QC.out.qc_reads
        CONCAT_STATS_SHORT_RAW(RUN_ILLUMINA_QC.out.raw_stats.map { cfg, stats -> stats }.collect().map { files -> tuple("short_reads_raw_seqstats", files)} )
        CONCAT_STATS_SHORT_QC(RUN_ILLUMINA_QC.out.qc_stats.map { cfg, stats -> stats }.collect().map { files -> tuple("short_reads_qc_seqstats", files)} )
        
    }
    if(!params.skip_long_reads_qc){
        RUN_NANOPORE_QC(long_reads)
        long_reads = RUN_NANOPORE_QC.out.qc_reads
        ch_software_versions = ch_software_versions.mix(RUN_NANOPORE_QC.out.versions)
        CONCAT_STATS_LONG_RAW(RUN_NANOPORE_QC.out.raw_stats.map { cfg, stats -> stats }.collect().map { files -> tuple("long_reads_raw_seqstats", files)} )
        CONCAT_STATS_LONG_QC(RUN_NANOPORE_QC.out.qc_stats.map { cfg, stats -> stats }.collect().map { files -> tuple("long_reads_qc_seqstats", files)} )
        
    }

    //classify
    if(!params.skip_kraken2){
        Channel
            .value(file( "${params.kraken2_db}" ))
            .set { ch_kraken2_db_file }
        KRAKEN2_KRAKEN2 (short_reads, ch_kraken2_db_file)
        ch_software_versions = ch_software_versions.mix(KRAKEN2_KRAKEN2.out.versions)
        BRACKEN(KRAKEN2_KRAKEN2.out.txt, ch_kraken2_db_file)
        ch_software_versions = ch_software_versions.mix(BRACKEN.out.versions)
    }
    // assembly
    if(!params.skip_short_reads_assembly && params.assembly_type == 'short'){
    
        RUN_ASSEMBLE_SHORT ( short_reads)
        contigs = RUN_ASSEMBLE_SHORT.out.contigs
        ch_software_versions = ch_software_versions.mix(RUN_ASSEMBLE_SHORT.out.versions)
        stats = RUN_ASSEMBLE_SHORT.out.stats
          
    } 
    if(!params.skip_long_reads_assembly && params.assembly_type == 'long'){
    
        RUN_ASSEMBLE_LONG ( long_reads, short_reads)
        contigs = RUN_ASSEMBLE_LONG.out.contigs
        ch_software_versions = ch_software_versions.mix(RUN_ASSEMBLE_LONG.out.versions)
        stats = RUN_ASSEMBLE_LONG.out.stats
        if(!params.skip_racon){
            RUN_RACON_POLISH(long_reads, contigs)
            contigs = RUN_RACON_POLISH.out.assembly
            ch_software_versions = ch_software_versions.mix(RUN_RACON_POLISH.out.versions)
            stats = RUN_RACON_POLISH.out.stats
        }
        if(!params.skip_medaka){
            MEDAKA(long_reads,  contigs)
            contigs = MEDAKA.out.assembly
            ch_software_versions = ch_software_versions.mix(MEDAKA.out.versions)
            //stats = MEDAKA.out.stats
        }

        if(!params.skip_short_reads_polish  && params.short_reads_polisher  == "pilon"){
            // 4x iterations are recommended
            RUN_PILON_POLISH(short_reads, contigs)
            contigs = RUN_PILON_POLISH.out.assembly
            ch_software_versions = ch_software_versions.mix(RUN_PILON_POLISH.out.versions)

            RUN_PILON_POLISH2(short_reads, contigs)
            contigs = RUN_PILON_POLISH2.out.assembly

            RUN_PILON_POLISH3(short_reads, contigs)
            contigs = RUN_PILON_POLISH3.out.assembly

            RUN_PILON_POLISH4(short_reads, contigs)
            contigs = RUN_PILON_POLISH4.out.assembly

            stats = RUN_PILON_POLISH4.out.stats

        }
        if(!params.skip_short_reads_polish  && params.short_reads_polisher  == "1xpolca+2xnextpolish"){
            //when run as sigularity, the program is depend on bwa, it is not working 
            /* if(params.bwa_on_path && !params.skip_polca){
                POLCA(short_reads, contigs)
                contigs = POLCA.out.assembly
            } */

            // 2x iterations are recommended
            if(!params.skip_nextpolish && params.nextpolish_path != null ){
                RUN_NEXTPOLISH_POLISH(short_reads, contigs )
                contigs = RUN_NEXTPOLISH_POLISH.out.assembly 
                ch_software_versions = ch_software_versions.mix(RUN_NEXTPOLISH_POLISH.out.versions)

                RUN_NEXTPOLISH_POLISH2(short_reads, contigs )
                contigs = RUN_NEXTPOLISH_POLISH2.out.assembly 
                stats = RUN_NEXTPOLISH_POLISH2.out.stats

            }
        }

       
    
    }
    if(!params.skip_hybrid_reads_assembly && params.assembly_type == 'hybrid'){
        long_reads.view()
        short_reads.view()
        RUN_ASSEMBLE_HYBRID ( long_reads, short_reads)
        contigs =RUN_ASSEMBLE_HYBRID.out.contigs
        ch_software_versions = ch_software_versions.mix(RUN_ASSEMBLE_HYBRID.out.versions)

        stats = RUN_ASSEMBLE_HYBRID.out.stats
        
    }
    CONCAT_STATS_ASM(stats.map { cfg, stats -> stats }.collect().map { files -> tuple("assembly_stats", files)} )
        
     // analysis
     if(params.annotation_tool== "bakta"){
        Channel
        .value(file( "${params.bakta_db}" ))
        .set { ch_bakta_db_file }

        BAKTA(contigs, ch_bakta_db_file)
        BAKTA_FEATURES(BAKTA.out.gff)
        ch_software_versions = ch_software_versions.mix(BAKTA.out.versions)
        CONCAT_STATS_BAKTA(BAKTA_FEATURES.out.feature_count.map { cfg, stats -> stats }.collect().map { files -> tuple("bakta", files)} )
        gff = BAKTA.out.gff 
        ffn = BAKTA.out.ffn
        faa = BAKTA.out.faa
    }
     else if(params.annotation_tool == "prokka"){
        PROKKA(contigs, [], [])
        PROKKA_FEATURES(PROKKA.out.gff)
        ch_software_versions = ch_software_versions.mix(PROKKA.out.versions)
        CONCAT_STATS_PROKKA(PROKKA_FEATURES.out.feature_count.map { cfg, stats -> stats }.collect().map { files -> tuple("prokka", files)} )
        gff=PROKKA.out.gff
        ffn = PROKKA.out.ffn
        faa = PROKKA.out.faa
    }
    
    //ARG(contigs, ffn, faa)
    if(!params.skip_amr ){
        Channel
            .value(file( "${params.amrfinderplus_db}" ))
            .set { ch_amrfinderplus_db_file }

        AMRFINDERPLUS(contigs, ch_amrfinderplus_db_file)
        ch_software_versions = ch_software_versions.mix(AMRFINDERPLUS.out.versions)
        CONCAT_AMR(AMRFINDERPLUS.out.tsv.map { cfg, amr -> amr }.collect().map { files -> tuple("amrfinderplus", files)} )
    }
    
    if(!params.skip_mlst){
        MLST (contigs)
        ch_software_versions = ch_software_versions.mix(MLST.out.versions)
        CONCAT_MLST (MLST.out.tsv.map { cfg, mlst -> mlst }.collect().map { files -> tuple("mlst", files)} )
    }
    if(!params.skip_mobsuite){
        MOBSUITE (contigs )
        ch_software_versions = ch_software_versions.mix(MOBSUITE.out.versions)
        CONCAT_MOBSUITE (MOBSUITE.out.plasmid.map { cfg, plasmid -> plasmid }.collect().map { files -> tuple("mobsuite", files)} )
    } 
    //[[id:sample1, single_end:false, genome_size:2.8m], /data/deve/dsl2_workflows/pathogen/work/c5/0433a6c7c3da2c34511b6ad1c34a74/sample1_abricate.tsv]         
    if(!params.skip_virulome){
        //virulome
        ABRICATE_VF(contigs)
        ch_software_versions = ch_software_versions.mix(ABRICATE_VF.out.versions)
        ABRICATE_VF.out.report.collect{ it[1] } | ABRICATE_SUMMARIZE_VF
    }
    
    // MODULE: Pipeline reporting
    
    
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_software_versions.unique().collectFile()
        //ch_software_versions.collectFile()
    )
    
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
