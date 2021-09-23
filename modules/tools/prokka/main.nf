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
    script:
        """
        prokka --metagenome --cpus ${task.cpus} --outdir $name --prefix $name $contigs
        """
}
