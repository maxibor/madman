process quast {
    tag "$name"

    label 'process_medium' 
    label 'process_ignore'

    publishDir "${params.outdir}/quast/${outdir}", mode: 'copy'

    conda (params.enable_conda ? "bioconda::quast=5.0.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/quast:5.0.2--py37pl526hb5aa323_2"
    } else {
        container "quay.io/biocontainers/quast:5.0.2--py37pl526hb5aa323_2"
    }

    input:
        tuple val(name), path(contigs)
        val(step)
    output:
        path("*_quast_${step}")
    script:
        outdir = name+"_quast_${step}"
        """
        quast -o $outdir -t ${task.cpus} $contigs
        """
}