process filter_contigs_length {
    tag "$name"

    label 'process_medium'
    label 'process_ignore'

    publishDir "${params.outdir}/fasta_filter/${name}", mode: 'copy'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
        tuple val(name), path(fasta)
    output:
        tuple val(name), path("*.size_filtered.fa")
    script:
        filter_out = name+".size_filtered.fa"
        """
        filter_contigs_length.py -min ${params.minlen} -p ${task.cpus} $fasta -o $filter_out
        """
}