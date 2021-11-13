
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'
params.options = [:]
options        = initOptions(params.options)

process GFF2FEATURES{

    tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
    
    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path('*.tsv'), emit: feature_count

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.suffix ? "${meta.id}.${options.suffix}" : "${meta.id}"
        """
        extract_info_from_gff.pl  -g ${gff} -n ${meta.id} > ${prefix}.feature_count.tsv
        """
   
}
