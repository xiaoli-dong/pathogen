// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ABRICATE {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::abricate=0.8.13" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/abricate%3A1.0.1--ha8f3691_1"
    } else {
        container "quay.io/biocontainers/abricate:1.0.1--ha8f3691_1"
    }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.tsv'), emit: report
    path ("versions.yml"), emit: version

    script:
    def prefix  = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    //def datadir = params.abricate_datadir ? "--datadir ${params.abricate_datadir}" : ""
    """
    [ ! -f  ${prefix}.fasta ] && ln -s $fasta ${prefix}.fasta
    abricate $options.args ${prefix}.fasta --threads $task.cpus  > ${prefix}.tsv
    #extract_gene_fasta.py ${prefix}.tsv $fasta

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(abricate --version 2>& 1 | sed 's/^abricate //;')
    END_VERSIONS

    """
}

process ABRICATE_SUMMARIZE {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }
    
    conda (params.enable_conda ? "bioconda::abricate=0.8.13" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/abricate%3A1.0.1--ha8f3691_1"
    } else {
        container "quay.io/biocontainers/abricate:1.0.1--ha8f3691_1"
    }

    input:
    path('?.tsv')

    output:
    path('all.summary.tsv'), emit: summary
    path ("versions.yml"), emit: version
    script:
    """
    abricate --summary *.tsv > all.summary.tsv

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(abricate --version 2>& 1 | sed 's/^abricate //;')
    END_VERSIONS

    """
}