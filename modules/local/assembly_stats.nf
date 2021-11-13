// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ASSEMBLY_STATS {
    tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'assembly-stats=1.0.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/assembly-stats%3A1.0.1--h7d875b9_4"
    } else {
        container "quay.io/biocontainers/assembly-stats:1.0.1--h7d875b9_4"
    }

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), file("*.tsv"), emit: stats
    path ("versions.yml"), emit: versions
    //refer to nullabour 
    script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
   // def tool = params.tool ? params.tool : ''
    def input = assembly.join('')
    
    """
    assembly-stats -t $input > ${prefix}_assembly_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(assembly-stats -v 2>&1 | sed -e 's/^Version: //;')
    END_VERSIONS

    """
}
