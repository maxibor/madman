process freebayes_consensus_calling {
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
        freebayes -j -f $contigs -p 1 $bam > calls.vcf
        bcftools view -i '%QUAL>=${params.min_variant_qual}' -Oz -o calls.vcf.gz calls.vcf
        bcftools index calls.vcf.gz

        cat $contigs | bcftools consensus calls.vcf.gz > ${name}_consensus_recalled.fa
        """    
}