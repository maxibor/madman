// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process filter_contigs_damage {
    tag "$name"

    label 'process_low'
    label 'process_ignore'

    publishDir "${params.outdir}/fasta_filter/${name}", mode: 'copy'

    conda (params.enable_conda ? "conda-forge::pandas=1.1.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pandas:1.1.5"
    } else {
        container "quay.io/biocontainers/pandas:1.1.5"
    }


    input:
        tuple val(name), path(pydamage_csv), path(contigs)
    output:
        tuple val(name), path("*.ancient_filtered.fa"), emit: fasta
        tuple val(name), path("*_ancient_contigs.txt"), emit: txt
    script:
        outfile = name + ".ancient_filtered.fa"
        ancient_contigs = name+"_ancient_contigs.txt"
        """
        filter_contigs_damage.py -d ${params.mindamage} --acc ${params.minaccuracy} -o $outfile $contigs $pydamage_csv
        awk -F "," '{if ((\$11 <= 0.05) && (\$8 >= ${params.mindamage}) && (\$1 >= ${params.minaccuracy})) { print \$1 }}' $pydamage_csv > $ancient_contigs
        """
}

