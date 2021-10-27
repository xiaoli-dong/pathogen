#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/pathogen
========================================================================================
    Github : https://github.com/nf-core/pathogen
    Website: https://nf-co.re/pathogen
    Slack  : https://nfcore.slack.com/channels/pathogen
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

//params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { PATHOGEN } from './workflows/pathogen'

//
// WORKFLOW: Run main nf-core/pathogen analysis pipeline
//
workflow NFCORE_PATHOGEN {
    PATHOGEN ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_PATHOGEN ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
