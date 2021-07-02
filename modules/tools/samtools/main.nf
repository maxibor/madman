process consensus_calling {
    tag "$name"

    label 'process_low'
    label 'process_ignore'

    publishDir "${params.outdir}/consensus_called_contigs/${name}", mode: 'copy'

    input:
        tuple val(name), path(contigs), path(bam)
    output:
        path("*.fa")
    script:
        """
        # call variants
        bcftools mpileup -Ou -f $contigs $bam | bcftools call -mv -Oz -o calls.vcf.gz
        bcftools index calls.vcf.gz

        cat $contigs | bcftools consensus calls.vcf.gz > ${name}_consensus_recalled.fa
        """    
}