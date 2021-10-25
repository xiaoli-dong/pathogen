// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)


process MINIASM {
    tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
    
    conda (params.enable_conda ? 'bioconda::miniasm=0.3_r179' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/miniasm:0.3_r179--h5bf99c6_2"
    } else {
        container "quay.io/biocontainers/miniasm:0.3_r179--h5bf99c6_2"
    }

    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(paf)

    output:
    tuple val(meta), path("*_miniasm.fasta") , emit: assembly
    tuple val(meta), path('*_miniasm.gfa') ,    emit: graph
    path  'versions.yml'                     , emit: versions

    script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    
    miniasm -f ${reads} ${paf} -p ug $options.args > ${prefix}_miniasm.gfa
    awk '/^S/{print ">"\$2"\\n"\$3}' "${prefix}_miniasm.gfa" | fold > ${prefix}_miniasm.fasta
    
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(miniasm -V 2>&1)
    END_VERSIONS

    
    """
}
