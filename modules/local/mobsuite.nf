// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'
module_dir = moduleDir + "/bin"
params.options = [:]
def options    = initOptions(params.options)

process MOBSUITE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    cache 'lenient'
    conda (params.enable_conda ? 'bioconda::shovill=1.1.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0'
    } else {
        container 'quay.io/biocontainers/fastp:0.20.1--h8b12597_0'
    }

    input:
    tuple val(meta), path(asm)

    output:
    tuple val(meta), path('*_contig_report.txt'), emit: contig_report
    //tuple val(meta), path('*_mobtyper_results.txt'), emit: mobs
    tuple val(meta), path('*_plasmid.txt'), emit: plasmid
    path ("versions.yml"), emit: version
    
    script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    maxmem = task.memory.toGiga()

    """
    source /data/software/miniconda3/etc/profile.d/conda.sh
    conda activate mob_suite
    mob_recon -i $asm -s ${prefix} -n $task.cpus -o mob

    if [ ! -f mob/mobtyper_results.txt ];then
        touch mob/mobtyper_results.txt
    fi
    
    mv mob/contig_report.txt ${prefix}_contig_report.txt
    mv mob/mobtyper_results.txt ${prefix}_mobtyper_results.txt
    
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(mob_recon -V | sed 's/^mob_recon //;')
    END_VERSIONS

    conda deactivate
    $module_dir/wrangle_mobsuite.py ${prefix}_mobtyper_results.txt ${prefix}
    mv plasmid.txt ${prefix}_plasmid.txt
    """
    
}
