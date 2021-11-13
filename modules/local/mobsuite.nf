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
    conda (params.enable_conda ? 'bioconda::mob_suite=3.0.3' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/mob_suite%3A3.0.3--pyhdfd78af_0'
    } else {
        container 'quay.io/biocontainers/mob_suite%3A3.0.3--pyhdfd78af_0'
    }

    input:
    tuple val(meta), path(asm)

    output:
    tuple val(meta), path('*report.txt'), emit: contig_report
    tuple val(meta), path('*mobtyper_results.txt'), emit: mobs
    tuple val(meta), path('*plasmid.txt'), emit: plasmid
    tuple val(meta), path('*.fasta'), emit: fasta
    
    path ("versions.yml"), emit: versions
    
    script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    maxmem = task.memory.toGiga()

    """
    mob_recon -i $asm -s ${prefix} -n $task.cpus -o mob

    if [ ! -f mob/mobtyper_results.txt ];then
        touch mob/mobtyper_results.txt
    fi
    
    mv mob/contig_report.txt ${prefix}.contig_report.txt
    mv mob/mobtyper_results.txt ${prefix}.mobtyper_results.txt
    mv mob/chromosome.fasta ${prefix}.chromosome.fasta
    wrangle_mobsuite.py ${prefix}.mobtyper_results.txt ${prefix}
    mv plasmid.txt ${prefix}.plasmid.txt

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(mob_recon -V | sed 's/^mob_recon //;')
    END_VERSIONS

    
    """
    
}
