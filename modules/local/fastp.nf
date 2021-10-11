// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process fastp {
    tag "$name"

    label 'process_medium'
    label 'process_mantory'

    conda (params.enable_conda ? "bioconda::fastp=0.20.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0"
    } else {
        container "quay.io/biocontainers/fastp:0.20.1--h8b12597_0"
    }

    input:
        tuple val(name), path(reads)

    output:
        tuple val(name), path("*.fq.gz"), emit: trimmed_reads
        path "*.json", emit: settings
        path  '*.version.txt', emit: version

    script:
        def software = getSoftwareName(task.process)
        if (!params.single_end) {
            out1 = name+".pair1.trimmed.fq.gz"
            out2 = name+".pair2.trimmed.fq.gz"
            """
            fastp --in1 ${reads[0]} --in2 ${reads[1]} --out1 $out1 --out2 $out2 -A -g --poly_g_min_len "${params.complexity_filter_poly_g_min}" -Q -L -w ${task.cpus} --json ${name}.json
            """
        } else {
            se_out = name+".trimmed.fq.gz"
            """
            fastp --in1 ${reads[0]} --out1 $se_out -A -g --poly_g_min_len "${params.complexity_filter_poly_g_min}" -Q -L -w ${task.cpus} --json ${name}.json
            echo \$(fastp --version 2>&1) | sed -e "s/fastp //g" > ${software}.version.txt
            """
        }
}
