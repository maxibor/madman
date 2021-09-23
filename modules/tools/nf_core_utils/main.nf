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

    output:
    path 'software_versions_mqc.yaml', emit: yaml
    path "software_versions.csv", emit: csv

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    AdapterRemoval --version  &> v_adapterremoval.txt 2>&1 || true
    fastp --version &> v_fastp.txt 2>&1 || true 
    megahit --version > v_megahit.txt
    spades.py --version > v_spades.txt
    quast --version > v_quast.txt
    pydamage --version > v_pydamage.txt
    damageprofiler --version > v_damageprofiler.txt
    prokka --version &> v_prokka.txt 2>&1 || true
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}
