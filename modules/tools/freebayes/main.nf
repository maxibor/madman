process freebayes_consensus_calling {
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
        path("*.fa")
    script:
        """
        freebayes -j -f $contigs -p 1 $bam > calls.vcf
        bcftools view -i '%QUAL>=${params.min_variant_qual}' -Oz -o calls.vcf.gz calls.vcf
        bcftools index calls.vcf.gz

        cat $contigs | bcftools consensus calls.vcf.gz > ${name}_consensus_recalled.fa
        """    
}