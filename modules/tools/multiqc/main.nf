include { getSoftwareName } from "$baseDir/modules/tools/nf_core_utils/main.nf"

process multiqc {
    tag "$custom_runName"

    label 'process_low'
    label 'process_mandatory'

    publishDir "${params.outdir}", mode: 'copy'

    conda (params.enable_conda ? "bioconda::multiqc=1.9" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/multiqc:1.9--pyh9f0ad1d_0"
    } else {
        container "quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    }

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
        path(workflow_summary) name "workflow_summary_mqc.yaml"
    output:
        path('*multiqc_report.html')
        path  "*.version.txt", emit: version
    script:
        def software = getSoftwareName(task.process)
        """
        multiqc -f -c $multiqc_conf --title $custom_runName ./
        multiqc --version | sed -e "s/multiqc, version //g" > ${software}.version.txt
        """
}