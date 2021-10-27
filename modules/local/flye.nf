// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'


params.options = [:]
options        = initOptions(params.options)

process FLYE {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
    
    conda (params.enable_conda ? 'bioconda::flye=v2.9' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/flye%3A2.9--py39h39abbe0_0"
    } else {
        container "quay.io/biocontainers/flye%3A2.9--py39h39abbe0_0"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*_contigs.fasta') , emit: assembly
    tuple val(meta), path('*contig_info.txt') , emit: info
    tuple val(meta), path('*_graph.gfa') ,    emit: graph
    path  'versions.yml',             emit: versions

    script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    flye $options.args -t ${task.cpus}  --nano-raw $reads -o assembly &> flye_${prefix}.log
    mv assembly/assembly.fasta ${prefix}_contigs.fasta
    mv assembly/assembly_info.txt ${prefix}_contig_info.txt
    mv assembly/assembly_graph.gfa ${prefix}_graph.gfa

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(flye --version)
    END_VERSIONS

    
    """

}
