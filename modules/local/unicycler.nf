// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process UNICYCLER {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::unicycler=0.4.8' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/unicycler:0.4.8--py38h8162308_3"
    } else {
        container "quay.io/biocontainers/unicycler:0.4.8--py38h8162308_3"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*_contigs.fasta'), emit: scaffolds
    tuple val(meta), path('*_graph.gfa'), emit: gfa
    tuple val(meta), path('*.log')         , emit: log
    path('versions.yml')                   , emit: versions

    script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def input_reads = ""

    if(params.assembly_type == 'short'){
        input_reads = meta.single_end ? "-s $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    }
    else if(params.assembly_type == 'long'){
        input_reads = "-l $reads"
    }
    
    """
    unicycler --threads $task.cpus $options.args $input_reads --out .
    mv assembly.fasta ${prefix}_contigs.fasta
    mv assembly.gfa ${prefix}_graph.gfa
    mv unicycler.log ${prefix}_unicycler.log

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(unicycler --version 2>&1 | sed 's/^.*Unicycler v//; s/ .*\$//')
    END_VERSIONS
    

    """
}
