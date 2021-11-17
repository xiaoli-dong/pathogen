// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
def options    = initOptions(params.options)
Channel
  .fromPath(params.amrfinderplus_db)
  .map { it -> [it.simpleName, it]}
  .set { bam_files }

process AMRFINDERPLUS {

 tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    cache 'lenient'
    conda (params.enable_conda ? 'bioconda::ncbi-amrfinderplus=3.10.18' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus%3A3.10.18--h17dc2d4_0'
    } else {
        container 'quay.io/biocontainers/ncbi-amrfinderplus:3.10.18--h17dc2d4_0'
    }

  input:
  tuple val(meta), path(fasta)
  path(db)

  output:
  //tuple val(meta), file("AMRFinder_resistance-only.tsv")
  tuple val(meta), path("*_amrfinderplus.tsv"), emit: tsv
  tuple val(meta), path("*_amrfinderplus.ffn"), emit: ffn
  path ("versions.yml"), emit: versions

//awk -F '\t' '{ if (\$3 != "") { print } }' AMRFinder_complete.tsv | grep -v "VIRULENCE" > AMRFinder_resistance-only.tsv ;
  script:
  def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    
  """
    amrfinder ${options.args} --nucleotide $fasta -o ${prefix}_amrfinderplus.tsv --name ${prefix} --threads $task.cpus --nucleotide_output ${prefix}_amrfinderplus.ffn --database $db

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(amrfinder --version 2>&1)
    END_VERSIONS

  """
}
