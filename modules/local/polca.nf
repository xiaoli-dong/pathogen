// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

//has dependency on bwa
process POLCA {

	tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'masurca=3.4.2' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/masurca%3A3.4.2--pl5262h86ccdc5_1"
        
        
    } else {
        container "quay.io/biocontainers/masurca%3A3.4.2--pl5262h86ccdc5_1"
       
    }

    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), file("*_polca.fasta"), emit: assembly
    tuple val(meta), file("*.report"), emit: report
    path ("versions.yml"), emit: versions

	when:
	//!params.skip_illumina
	script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
	"""
	polca.sh -a ${assembly} -r '${reads[0]} ${reads[1]}' -t $task.cpus 
    mv ${assembly}.report ${prefix}_polca.report
    mv ${assembly}.PolcaCorrected.fa ${prefix}_polca.fasta
	
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(masurca --version 2>&1 | sed 's/^version //;')
    END_VERSIONS
    
    
	"""
}
