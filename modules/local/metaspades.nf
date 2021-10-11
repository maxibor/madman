// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process metaspades {
    tag "$name"

    label 'process_very_bigmem'

    publishDir "${params.outdir}/assembly/metaspades/${name}", mode: 'copy'

    conda (params.enable_conda ? "bioconda::spades=3.9.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/spades:3.9.1--h9ee0642_1"
    } else {
        container "quay.io/biocontainers/spades:3.9.1--h9ee0642_1"
    }

    input:
        tuple val(name), path(reads)
    output:
        tuple val(name), path("*_metaspades.contigs.fa"), emit: contigs
        tuple val(name), path("${name}/*"), emit: metaspades_logs
        path  '*.version.txt', emit: version

    script:
        mem = task.memory.toGiga()
        def software = getSoftwareName(task.process)
        """
        spades.py --meta \
                  -1 ${reads[0]} \
                  -2 ${reads[1]} \
                  -t ${task.cpus} \
                  -m $mem \
                  --phred-offset ${params.phred} \
                  -o $name
        cp $name/contigs.fasta ${name}_metaspades.contigs.fa
        echo \$(spades.py --version 2>&1) | sed 's/^.*SPAdes genome assembler v//; s/ .*\$//' > ${software}.version.txt
        """
}
