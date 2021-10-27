// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'
module_dir = moduleDir + "/bin"
params.options = [:]
options        = initOptions(params.options)

process PILON {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'pilon=1.24' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pilon%3A1.24--hdfd78af_0"
    } else {
        container "quay.io/biocontainers/pilon%3A1.24--hdfd78af_0"
    }

    input:
    tuple val(meta), path(sorted_bam)
    tuple val(meta), path(sorted_bamb_idex)
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*_pilon.fasta") , emit: assembly
    path  'versions.yml'                     , emit: version

    script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def maxmem = "-Xmx${task.memory.toGiga()}g"
    def round = params.pilon_round ? params.pilon_round : ""

    """
    pilon $options.args --genome ${assembly} --frags ${sorted_bam} --output round${round}_pilon 

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(pilon --version 2>&1 | sed 's/^Pilon version //; s/ .*\$//' )
    END_VERSIONS
    """
}
