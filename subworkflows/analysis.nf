#!/usr/bin/env nextflow

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include { ABRITAMR} from './../modules/resistome/main'          addParams( options: modules['abritamr'])      
include { MLST } from './../modules/mlst/main'                  addParams( options: modules['mlst'])
include { PROKKA } from './../modules/prokka/main'              addParams( options: modules['prokka'])
include { BAKTA } from './../modules/bakta/main'                addParams( options: modules['bakta'])
include { MOBSUITE } from './../modules/mobsuite/main'          addParams( options: modules['mobsuite'])

workflow RUN_ANALYSIS {   

    take:
        contigs
    main:
        ABRITAMR (contigs )
        MLST (contigs )
        if(params.annotator == "bakta"){
            BAKTA(contigs)
            gff = BAKTA.out.bakta_gff

        }
        else if(params.annotator == "prokka"){
            PROKKA(contigs)
            gff = PROKKA.out.gff
            //gff_summary = PROKKA.out.txt
        }

        MOBSUITE (contigs )

    emit:
        resistome = ABRITAMR.out.resistome
        virulome = ABRITAMR.out.abritamr_virulence
        mlst = MLST.out.mlst
        //gff = PROKKA.out.gff
        gff
        //prokka_txt = PROKKA.out.prokka_txt
        plasmid = MOBSUITE.out.plasmid
}
