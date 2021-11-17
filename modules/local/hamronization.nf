// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process HAMRONIZE_ABRICATE {
    tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    
    conda (params.enable_conda ? "bioconda::hamronization=1.0.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/hamronization%3A1.0.3--py_0"
    } else {
        container "quay.io/biocontainers/hamronization:1.0.3--py_0"
    }

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path('*.tsv'), emit: hamronized
    path "versions.yml"                    , emit: versions

    script:
    def prefix  = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    hamronize abricate ${report} --reference_database_version db_v_1 --analysis_software_version tool_v_1 --output ${prefix}_abricate_hamronized.tsv
    
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(hamronize -v | sed 's/hamronize //g' )
    END_VERSIONS

    """
}

process HAMRONIZE_AMRFINDERPLUS {
    tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    
    conda (params.enable_conda ? "bioconda::hamronization=1.0.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/hamronization%3A1.0.3--py_0"
    } else {
        container "quay.io/biocontainers/hamronization:1.0.3--py_0"
    }

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path('*.tsv'), emit: hamronized
    path "versions.yml"                    , emit: versions

    script:
    def prefix  = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    hamronize \\
        amrfinderplus \\
        ${report} \\
        --reference_database_version db_v_1 \\
        --analysis_software_version tool_v_1 \\
        --input_file_name aafasta \\
        --output ${prefix}_amrfinderplus_hamronized.tsv
    
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(hamronize -v | sed 's/hamronize //g' )
    END_VERSIONS

    """
}
process HAMRONIZE_SRAX {
    //tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:"") }

    conda (params.enable_conda ? "bioconda::hamronization=1.0.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/hamronization%3A1.0.3--py_0"
    } else {
        container "quay.io/biocontainers/hamronization:1.0.3--py_0"
    }
    

    input:
    path(report)

    output:
    path('*.tsv'), emit: hamronized
    path "versions.yml"                    , emit: versions

    script:
    //def prefix  = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    hamronize srax ${report} --input_file_name mysrax --reference_database_version db_v_1 --reference_database_id db_id_1 --analysis_software_version tool_v_1 --output srax_hamronized.tsv
    
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(hamronize -v | sed 's/hamronize //g' )
    END_VERSIONS

    """
}

process HAMRONIZE_RGI {
    tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::hamronization=1.0.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/hamronization%3A1.0.3--py_0"
    } else {
        container "quay.io/biocontainers/hamronization:1.0.3--py_0"
    }
    

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path('*.tsv'), emit: hamronized
    path "versions.yml"                    , emit: versions

    script:
    def prefix  = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    hamronize rgi ${report} --input_file_name ${prefix} --reference_database_version db_v_1 --analysis_software_version tool_v_1 --output ${prefix}_rgi_hamronized.tsv
    
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(hamronize -v | sed 's/hamronize //g' )
    END_VERSIONS

    """
}

process HAMRONIZE_SUMMARIZE {
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    conda (params.enable_conda ? "bioconda::hamronization=1.0.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/hamronization%3A1.0.3--py_0"
    } else {
        container "quay.io/biocontainers/hamronization:1.0.3--py_0"
    }
    

    input:
    path('?.tsv')

    output:
    path('amr_hamronized_summary.tsv'), emit: summary_tsv
    path('amr_hamronized_summary.html'), emit: summary_html
    path "versions.yml"                    , emit: versions

    script:
    """
    hamronize summarize --output amr_hamronized_summary.tsv --summary_type tsv *.tsv
    hamronize summarize --output amr_hamronized_summary.html --summary_type interactive *.tsv

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(hamronize -v | sed 's/hamronize //g' )
    END_VERSIONS

    """
}
