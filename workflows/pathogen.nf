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
def checkPathParamList = [ params.input, params.multiqc_config]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
//include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )
include { INPUT_CHECK } from '../subworkflows/local/pathogen_input_check' addParams( options: [:] )
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

include {MEDAKA} from '../modules/local/medaka' addParams( options: modules['medaka'])
include {POLCA} from '../modules/local/polca'                              addParams( options: modules['polca']) 
include {KRAKEN2_KRAKEN2} from '../modules/nf-core/modules/kraken2/kraken2/main'        addParams( options: modules['kraken2'])   

//include { ABRITAMR} from '../modules/local/abritamr'          addParams( options: modules['abritamr'])      
include { MLST } from '../modules/nf-core/modules/mlst/main'                  addParams( options: modules['mlst'])
include { PROKKA } from '../modules/nf-core/modules/prokka/main'              addParams( options: modules['prokka'])
include { BAKTA } from '../modules/local/bakta'                addParams( options: modules['bakta'])
include { MOBSUITE } from '../modules/local/mobsuite'          addParams( options: modules['mobsuite'])
include {
    ARG
} from '../subworkflows/local/arg'

include { ABRICATE} from '../modules/local/abricate' addParams( options: modules['abricate_vf'])
include { ABRICATE_SUMMARIZE } from '../modules/local/abricate' 
/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/



def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
//include { FASTQC  } from '../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )

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
    short_reads = INPUT_CHECK.out.shortreads
    long_reads = INPUT_CHECK.out.longreads

    RUN_ILLUMINA_QC(short_reads)
    ch_software_versions = ch_software_versions.mix(RUN_ILLUMINA_QC.out.versions.first().ifEmpty(null))
    short_reads = RUN_ILLUMINA_QC.out.qc_reads
    
    
    RUN_NANOPORE_QC(long_reads)
    ch_software_versions = ch_software_versions.mix(RUN_NANOPORE_QC.out.versions.first().ifEmpty(null))
    long_reads = RUN_NANOPORE_QC.out.qc_reads

    //classify
    kraken2_db = Channel.fromPath( "${params.kraken2_db}")
    KRAKEN2_KRAKEN2 (short_reads, kraken2_db)

    // assembly
    if(params.assembly_type == 'short'){
    
        RUN_ASSEMBLE_SHORT ( short_reads)
        contigs = RUN_ASSEMBLE_SHORT.out.contigs
    } 
    else if(params.assembly_type == 'long'){
    
        RUN_ASSEMBLE_LONG ( long_reads, short_reads)
        contigs = RUN_ASSEMBLE_LONG.out.contigs

        if(!params.skip_racon){
            RUN_RACON_POLISH(long_reads, contigs)
            contigs = RUN_RACON_POLISH.out.assembly
        }
        if(!params.skip_medaka){
            MEDAKA(long_reads,  contigs)
            contigs = MEDAKA.out.assembly
        }

        if(!params.skip_short_reads_polish  && params.short_reads_polisher  == "pilon"){
            // 4x iterations are recommended
            RUN_PILON_POLISH(short_reads, contigs)
            contigs = RUN_PILON_POLISH.out.assembly

            RUN_PILON_POLISH2(short_reads, contigs)
            contigs = RUN_PILON_POLISH2.out.assembly

            RUN_PILON_POLISH3(short_reads, contigs)
            contigs = RUN_PILON_POLISH3.out.assembly

            RUN_PILON_POLISH4(short_reads, contigs)
            contigs = RUN_PILON_POLISH4.out.assembly

        }
        if(!params.skip_short_reads_polish  && params.short_reads_polisher  == "1xpolca+2xnextpolish"){
            //when run as sigularity, the program is depend on bwa, it is not working 
            /* if(params.bwa_on_path && !params.skip_polca){
                POLCA(short_reads, contigs)
                contigs = POLCA.out.assembly
            } */
run
            // 2x iterations are recommended
            if(!params.skip_nextpolish || params.nextpolish_path != null ){
                RUN_NEXTPOLISH_POLISH(short_reads, contigs )
                contigs = RUN_NEXTPOLISH_POLISH.out.assembly 

                RUN_NEXTPOLISH_POLISH2(short_reads, contigs )
                contigs = RUN_NEXTPOLISH_POLISH2.out.assembly 

            }
        }
    
    }
    else if(params.assembly_type == 'hybrid'){
        long_reads.view()
        short_reads.view()
        RUN_ASSEMBLE_HYBRID ( long_reads, short_reads)
        contigs =RUN_ASSEMBLE_HYBRID.out.contigs
        
    }
    
     // analysis
     if(params.annotation_tool== "bakta"){
        bakta_db = Channel.fromPath( "${params.bakta_db}")
        BAKTA(contigs, bakta_db)
     }
     else if(params.annotation_tool == "prokka"){
        PROKKA(contigs, [], [])
    }
    card_db = Channel.fromPath( "${params.card_db}")
    ARG(contigs, card_db)
    MLST (contigs)
    MOBSUITE (contigs )

    //virulome
    ABRICATE(contigs)
    ABRICATE.out.report.collect{ it[1] } | ABRICATE_SUMMARIZE
    //
    // MODULE: Run FastQC
    //
   /*  FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_software_versions = ch_software_versions.mix(FASTQC.out.versions.first().ifEmpty(null))
 */
    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
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
