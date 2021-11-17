// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BAKTA {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    cache 'lenient'
    conda (params.enable_conda ? 'bioconda::bakta=1.2.2' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/bakta%3A1.2.2--pyhdfd78af_0'
    } else {
        container 'quay.io/biocontainers/bakta:1.2.2--pyhdfd78af_0'
    }

    scratch true
    
    input:
    tuple val(meta), path(contigs)
    path db

    output:
    tuple val(meta), path('*.gff3'), emit: gff
    tuple val(meta), path('*.json'), emit: json
    tuple val(meta), path('*.tsv'), emit: tsv
    tuple val(meta), path('*.gbff'), emit: gbff
    tuple val(meta), path('*.embl'), emit: embl
    tuple val(meta), path('*.ffn'), emit: ffn
    tuple val(meta), path("${meta.id}.faa"), emit: faa
    path ("versions.yml"), emit: versions
    
    
    script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    
    """
    
    bakta $options.args --db ${db} --output ./ --prefix ${prefix} --locus-tag ${prefix} --threads $task.cpus $contigs
    
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(bakta --version 2>& 1 | sed 's/^bakta //;')
    END_VERSIONS
    
    
    """
    
}
