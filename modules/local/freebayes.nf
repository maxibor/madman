// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process freebayes {
    tag "$name"

    label 'process_low'
    label 'process_ignore'

    publishDir "${params.outdir}/consensus_called_contigs/${name}", mode: 'copy'

    conda (params.enable_conda ? "bioconda::freebayes=1.3.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/freebayes:1.3.5--py39hd2e4403_0"
    } else {
        container "quay.io/biocontainers/freebayes:1.3.5--py39hd2e4403_0"
    }

    input:
        tuple val(name), path(contigs), path(bam)
    output:
        tuple val(name), path("*.vcf"), emit: vcf
        path  '*.version.txt', emit: version
    script:
        def software = getSoftwareName(task.process)
        """
        freebayes -j -f $contigs -p 1 $bam > calls.vcf
        echo \$(freebayes --version 2>&1) | sed 's/^.*version:\sv// > ${software}.version.txt
        """
}
