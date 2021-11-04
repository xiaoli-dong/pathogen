// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MEDAKA {
    tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'medaka=1.4.3-0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/medaka:1.4.3--py38h130def0_0"
    } else {
        container "quay.io/biocontainers/medaka:1.4.3--py38h130def0_0"
    }

    input:
    tuple val(meta), path(longreads)
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), file("*_medaka.fasta"), emit: assembly
    path ("versions.yml"), emit: versions

    script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    
    """
    medaka_consensus -i $longreads -d $assembly  -t $task.cpus $options.args
    mv medaka/consensus.fasta ${prefix}_medaka.fasta

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(medaka --version 2>&1 | sed -e 's/^medaka -v//;' | sed '/^[[:space:]]*\$/d')
    END_VERSIONS

    """
}
