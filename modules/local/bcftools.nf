// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process bcftools {
    tag "$name"

    label 'process_low'
    label 'process_ignore'

    publishDir "${params.outdir}/consensus_called_contigs/${name}", mode: 'copy'

    conda (params.enable_conda ? 'bioconda::bcftools=1.13' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/bcftools:1.13--h3a49de5_0'
    } else {
        container 'quay.io/biocontainers/bcftools:1.13--h3a49de5_0'
    }

    input:
        tuple val(name), path(contigs), path(vcf)
    output:
        tuple val(name), path("*.fa"), emit: contigs
        path  "*.version.txt", emit: version
    script:
        def software = getSoftwareName(task.process)
        """
        bcftools view -i '%QUAL>=${params.min_variant_qual}' -Oz -o calls.vcf.gz calls.vcf
        bcftools index calls.vcf.gz

        cat $contigs | bcftools consensus calls.vcf.gz > ${name}_consensus_recalled.fa
        echo \$(bcftools --version 2>&1) | sed 's/^.*bcftools //; s/ .*\$//' > ${software}.version.txt
        """
}
