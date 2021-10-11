include { getSoftwareName } from "$baseDir/modules/tools/nf_core_utils/main.nf"

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
        path '*_fastqc.{zip,html}', emit: logs
        path  "*.version.txt", emit: version
    script:
        def software = getSoftwareName(task.process)
        """
        fastqc -t ${task.cpus} -q $reads
        fastqc --version | sed -e "s/FastQC v//g"  > ${software}.version.txt
        """
}