// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)



process SPADESHYBRID {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::spades=3.15.3' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/spades:3.15.3--h95f258a_0"
    } else {
        container "quay.io/biocontainers/spades:3.15.3--h95f258a_0"
    }

    input:
    tuple val(meta), path(long_reads)
    tuple val(meta), path(short_reads)

   output:
    tuple val(meta), path('*scaffolds.fasta'), emit: assembly
    tuple val(meta), path('*contigs.fasta'), emit: contigs
    tuple val(meta), path('*.gfa'),  emit: gfa
    tuple val(meta), path('*.log'), emit: log
    //path("versions.yml"), emit: versions
    path  "versions.yml"                       , emit: versions

    script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def input_reads = meta.single_end ? "-s $short_reads --nanopore ${long_reads}" : "-1 ${short_reads[0]} -2 ${short_reads[1]} --nanopore ${long_reads}"

    
    maxmem = task.memory.toGiga()
    // if ( params.spadeshybrid_fix_cpus == -1 || task.cpus == params.spadeshybrid_fix_cpus )
        """
        spades.py $options.args --threads $task.cpus --memory $maxmem $input_reads -o spades 

        mv spades/scaffolds.fasta ${prefix}_scaffolds.fasta
        mv spades/contigs.fasta ${prefix}_contigs.fasta
        mv spades/assembly_graph_with_scaffolds.gfa ${prefix}_graph.gfa
        mv spades/spades.log ${prefix}.log

       cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(spades.py --version 2>&1 | sed 's/SPAdes genome assembler //' )
        END_VERSIONS

        """
    // else
    //     error "ERROR: '--spadeshybrid_fix_cpus' was specified, but not succesfully applied. Likely this is caused by changed process properties in a custom config file."
}
