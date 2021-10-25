// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SKESA {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
    
    conda (params.enable_conda ? 'bioconda::skesa=v2.4.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/skesa%3A2.4.0--he1c1bb9_0"
    } else {
        container "quay.io/biocontainers/skesa%3A2.4.0--he1c1bb9_0"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*_contigs.fasta'), emit: contigs
    path ("versions.yml"), emit: versions

    script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def input_reads = meta.single_end ? "-s $reads" : "${reads[0]},${reads[1]}"
    maxmem = task.memory.toGiga()

    """
    
    skesa --reads $input_reads --cores $task.cpus --memory $maxmem --contigs_out assembly.fasta
    mv assembly.fasta ${prefix}_contigs.fasta
    
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(skesa -v 2>&1 | sed -e 's/^skesa -v//;;s/^SKESA //;' | sed '/^[[:space:]]*\$/d')
    END_VERSIONS
    
    
    """
    
}
