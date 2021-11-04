// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RGI {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "python=3.6 bioconda::rgi=5.1.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/rgi%3A5.1.1--py_0"
    } else {
        container "quay.io/biocontainers/rgi:5.1.1--py_0"
    }

    input:
    tuple val(meta), path(fasta)
    path card_db

    output:
    tuple val(meta), path('*_rgi.json'), emit: json
    tuple val(meta), path('*_rgi.txt'), emit: txt
    path "versions.yml"                    , emit: versions

    script:
    def prefix  = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    # Place card_rgi source in a read/write location for container
    #mkdir card_temp && cp -r /opt/conda/lib/python3.6/site-packages/app/ card_temp
    mkdir card_temp && cp -r /usr/local/lib/python3.6/site-packages/app/ card_temp
    export PYTHONPATH="\$(pwd)/card_temp/:\$PATH"

    rgi load --card_json ${card_db} --local
    rgi main -i $fasta -o ${prefix}_rgi -n $task.cpus $options.args
    
    #clean up work dir, if it exists
    [[ -d card_temp ]] && rm -r card_temp

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(rgi main --version | sed 's/rgi //g')
    END_VERSIONS 
    """
}

process RGI_HEATMAP {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "python=3.6 bioconda::rgi=5.1.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/rgi%3A5.1.1--py_0"
    } else {
        container "quay.io/biocontainers/rgi:5.1.1--py_0"
    }

    input:
    //path json
    path('?.json')

    output:
    path('*.png'), emit: heatmap
    path('rgi_heatmap.eps'), emit: eps
    path('rgi_heatmap.csv'),   emit: csv
    path "versions.yml"                    , emit: versions

    script:
    """
    mkdir dir
    cp ${json.join(' ')} dir
    rgi heatmap -i dir -o rgi_heatmap

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(rgi main --version | sed 's/rgi //g')
    END_VERSIONS 
    """
}
