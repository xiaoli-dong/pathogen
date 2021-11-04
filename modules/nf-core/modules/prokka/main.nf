include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PROKKA {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::prokka=1.14.6" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/prokka:1.14.6--pl526_0"
    } else {
        container "quay.io/biocontainers/prokka:1.14.6--pl526_0"
    }

    input:
    tuple val(meta), path(fasta)
    path proteins
    path prodigal_tf

    output:
    tuple val(meta), path("*.gff"), emit: gff
    tuple val(meta), path("*.gbk"), emit: gbk
    tuple val(meta), path("*.fna"), emit: fna
    tuple val(meta), path("*.faa"), emit: faa
    tuple val(meta), path("*.ffn"), emit: ffn
    tuple val(meta), path("*.sqn"), emit: sqn
    tuple val(meta), path("*.fsa"), emit: fsa
    tuple val(meta), path("*.tbl"), emit: tbl
    tuple val(meta), path("*.err"), emit: err
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*.txt"), emit: txt
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml" , emit: versions

    script:
    prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    def prodigal_opt = prodigal_tf ? "--prodigaltf ${prodigal_tf[0]}" : ""
    """
    prokka \\
        $options.args \\
        --cpus $task.cpus \\
        --prefix $prefix \\
        $proteins_opt \\
        $prodigal_tf \\
        $fasta
    mv ${prefix}/* .

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(prokka --version 2>&1) | sed 's/^.*prokka //')
    END_VERSIONS
    """
}
