process multiqc {
    
    label 'process_low'

    publishDir "${params.results}", mode: 'copy'

    input:
        path('adapterRemoval/*')
        path('fastqc/*')
        path('quast_pre/*')
        path('quast_post/*')
        path('prokka/*')
        path('fastp/*')
        path('dmgprof_pre/*')
        path('dmgprof_post/*')
        path(multiqc_conf)
    output:
        path('multiqc_report.html')
    script:
        """
        multiqc -c $multiqc_conf .
        """
}