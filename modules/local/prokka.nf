// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process prokka {
    tag "$name"

    label 'process_high'
    label 'process_ignore'

    publishDir "${params.outdir}/prokka", mode: 'copy'

    conda (params.enable_conda ? "bioconda::prokka=1.14.6" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/prokka:1.14.6--pl526_0"
    } else {
        container "quay.io/biocontainers/prokka:1.14.6--pl526_0"
    }

    input:
        tuple val(name), path(contigs)
    output:
        path("${name}")
        path "*.version.txt", emit: version
    script:
        def software = getSoftwareName(task.process)
        """
        prokka --metagenome --cpus ${task.cpus} --outdir $name --prefix $name $contigs
        echo \$(prokka --version 2>&1) | sed 's/^.*prokka //' > ${software}.version.txt
        """
}
