process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    conda (params.enable_conda ? "conda-forge::markown=2.6.9" : null)
    container "datafolklabs/markdown:latest"

    input:
    path output_docs

    output:
    path "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }
    input:
    path versions

    output:
    path "software_versions.tsv"     , emit: tsv
    path 'software_versions_mqc.yaml', emit: yaml

    script: // This script is bundled with the pipeline, in nf-core/madman/bin/
    """
    echo $workflow.manifest.version > pipeline.version.txt
    echo $workflow.nextflow.version > nextflow.version.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

def getSoftwareName(task_process) {
    return task_process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()
}
