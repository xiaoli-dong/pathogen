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
    tuple val(meta), path(sorted_bam)
    tuple val(meta), path(sorted_bamb_idex)
    tuple val(meta), path(contigs)

	output:
    tuple val(meta), file("*_nextpolish.fasta"), emit: assembly
    path ("versions.yml"), emit: versions

	// when:
	// !params.skip_illumina
	script:
	def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def algorithm = params.algorithm ? params.algorithm : "1"

	"""
    python ${params.nextpolish_path}/lib/nextpolish1.py $options.args -g ${contigs} -t ${algorithm} -p ${task.cpus} -s ${sorted_bam} >{$prefix}_t${algorithm}_nextpolish.fasta
       
	cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
    	${getSoftwareName(task.process)}: \$(${params.nextpolish_path}/nextPolish -v 2>&1 | sed -e 's/^nextPolish //;' | sed '/^[[:space:]]*\$/d')
    END_VERSIONS

	"""
}
