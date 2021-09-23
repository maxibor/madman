process megahit {
    tag "MEGAHIT - $name"

    label 'process_high'
    label 'process_mandatory'

    publishDir "${params.outdir}/assembly/megahit/${name}", mode: 'copy'

    conda (params.enable_conda ? "bioconda::megahit=1.2.9" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/megahit:1.2.9--h8b12597_0"
    } else {
        container "quay.io/biocontainers/megahit:1.2.9--h8b12597_0"
    }


    input:
        tuple val(name), path(reads)

    output:
        tuple val(name), path("megahit_out/*_megahit.contigs.fa"), emit: contigs
        tuple val(name), path("megahit_out/*.log"), emit: log
    script:
        mem = task.memory.toBytes()
        """
        megahit -1 ${reads[0]} -2 ${reads[1]} -t ${task.cpus} -m $mem --out-prefix $name
        mv megahit_out/${name}.contigs.fa megahit_out/${name}_megahit.contigs.fa
        """       
}