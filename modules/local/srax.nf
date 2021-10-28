// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SRAX {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    //conda (params.enable_conda ? "lgpdevtools::srax" : null)
    conda (params.enable_conda ? "bioconda::srax=1.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/srax%3A1.5--pl5262ha8f3691_1"
    } else {
        container "quay.io/biocontainers/srax:1.5--pl5262ha8f3691_1"
    }
    input:
        path fasta
    output:
        path('srax_results'), emit: result
        path ("versions.yml"), emit: version

    script:
    //def prefix  = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    mkdir srax_in
    cp ${fasta.join(' ')} srax_in
    sraX $options.args -i srax_in -o srax_out
    mv srax_out srax_results

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(sraX --version | sed '/^[[:space:]]*\$/d; /sraX -.*\$/d; /Copyright.*\$/d; s/.*version: sraXv//')
    END_VERSIONS

    """
}
