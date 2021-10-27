// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)


process RACON {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'racon=1.4.20-1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/racon:1.4.20--h9a82719_1"
    } else {
        container "quay.io/biocontainers/racon:1.4.20--h9a82719_1"
    }

    input:
    //tuple val(meta), val(reads), file(longreads), path('assembly.fasta'), path(paf)
    tuple val(meta), path(longreads)
    tuple val(meta), path(paf)
    tuple val(meta), path(assembly)


    output:
    tuple val(meta), path('*_racon.fasta') , emit: assembly
    path  'versions.yml'                     , emit: versions

    script:
    def round = params.racon_round ? params.racon_round : ""
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    racon -t ${task.cpus} ${longreads} ${paf} ${assembly} > ${prefix}_${round}_racon.fasta

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(racon --version 2>&1)
    END_VERSIONS
    
    """
}