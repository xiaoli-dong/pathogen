// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

module_dir = moduleDir + "/bin"

params.options = [:]
def options    = initOptions(params.options)


process ABRITAMR {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
    
    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("*summary_virulence.txt"), emit: abritamr_virulence
    //tuple val(meta),path('*summary_matches.txt'), emit: abritamr_matches
    //tuple val(meta),path('*summary_partials.txt'), emit: abritamr_partials
    tuple val(meta),path('*_resistome.txt'), emit: resistome
    path ("versions.yml"), emit: version
    
    script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    //def organism = params.species_options.any { it.contains(params.species) } ? "-sp $params.species": ""
    """
    source /data/software/miniconda3/etc/profile.d/conda.sh
    conda activate ncbi-amrfinderplus
    abritamr run -c $contigs -px ${meta.id} -j $task.cpus $options.args

    cp ${meta.id}/summary_virulence.txt ${prefix}_summary_virulence.txt
    cp ${meta.id}/summary_matches.txt ${prefix}_summary_matches.txt
    cp ${meta.id}/summary_partials.txt ${prefix}_summary_partials.txt

    
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(abritamr -v 2>&1 | sed 's/abritamr //')
    END_VERSIONS
    
    conda deactivate
    $module_dir/combine.py ${meta.id} ${prefix}_summary_matches.txt ${prefix}_summary_partials.txt
    mv resistome.txt ${prefix}_resistome.txt
    """
}
