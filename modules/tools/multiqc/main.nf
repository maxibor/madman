process multiqc {
    tag "$custom_runName"

    label 'process_low'

    publishDir "${params.outdir}", mode: 'copy'

    input:
        path('adapterRemoval/*')
        path('fastqc/*')
        path('fastp/*')
        path('quast_pre/*')
        path('quast_post/*')
        path('dmgprof_pre/*')
        path('dmgprof_post/*')
        path('prokka/*')
        path('software_versions/*')
        path(multiqc_conf)
        val(custom_runName)
        file workflow_summary name "workflow_summary_mqc.yaml"
    output:
        path('*multiqc_report.html')
    script:
        """
        multiqc -f -c $multiqc_conf --title $custom_runName ./
        """
}