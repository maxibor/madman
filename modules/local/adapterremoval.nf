// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process adapterremoval {
    tag "$name"

    label 'process_medium'
    label 'process_mandatory'

    conda (params.enable_conda ? "bioconda::adapterremoval=2.3.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/adapterremoval:2.3.2--h33c0355_1"
    } else {
        container "quay.io/biocontainers/adapterremoval:2.3.2--h33c0355_1"
    }

    input:
        tuple val(name), path(reads)
        path(adapter_list)

    output:
        tuple val(name), path('*.trimmed.fastq'), emit: trimmed_reads
        path "*.settings", emit: settings
        path "*.version.txt", emit: version

    script:
        def software = getSoftwareName(task.process)
        out1 = name+".pair1.trimmed.fastq"
        out2 = name+".pair2.trimmed.fastq"
        se_out = name+".trimmed.fastq"
        settings = name+".settings"
        if (! params.single_end){
            """
            AdapterRemoval --basename $name \
                           --file1 ${reads[0]} \
                           --file2 ${reads[1]} \
                           --trimns \
                           --trimqualities \
                           --minquality 20 \
                           --minlength 30 \
                           --output1 $out1 \
                           --output2 $out2 \
                           --threads ${task.cpus} \
                           --qualitybase ${params.phred} \
                           --adapter-list ${adapter_list} \
                           --settings $settings
            AdapterRemoval --version 2>&1 | sed -e "s/AdapterRemoval ver. //g" > ${software}.version.txt
            """
        } else {
            """
            AdapterRemoval --basename $name \
                           --file1 ${reads[0]} \
                           --trimns \
                           --trimqualities \
                           --minquality 20 \
                           --minlength 30 \
                           --output1 $se_out \
                           --threads ${task.cpus} \
                           --qualitybase ${params.phred} \
                           --adapter-list ${adapter_list} \
                           --settings $settings
            AdapterRemoval --version 2>&1 | sed -e "s/AdapterRemoval ver. //g" > ${software}.version.txt
            """
        }

}
