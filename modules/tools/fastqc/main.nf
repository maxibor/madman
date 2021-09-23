process fastqc {
    tag "$name"

    label 'process_low'
    label 'process_mandatory'

    conda (params.enable_conda ? "bioconda::fastqc=0.11.9" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1"
    }
    
    input:
        tuple val(name), path(reads)
    output:
        path '*_fastqc.{zip,html}'
    script:
        """
        fastqc -t ${task.cpus} -q $reads
        """
}