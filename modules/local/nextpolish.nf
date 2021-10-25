// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'
params.options = [:]
options        = initOptions(params.options)


process NEXTPOLISH {

	tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    
    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(assembly)
	output:
    tuple val(meta), file("*_nextpolish.fasta"), emit: assembly
    path ("versions.yml"), emit: versions

	when:
	!params.skip_illumina
	script:
	def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

	"""
	nextpolish_sr.sh ${assembly} ${reads[0]} ${reads[1]} $task.cpus $options.args
	
	cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
    	${getSoftwareName(task.process)}: \$(nextPolish -v 2>&1 | sed -e 's/^nextPolish //;' | sed '/^[[:space:]]*\$/d')
    END_VERSIONS

	"""
}
