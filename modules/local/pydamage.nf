// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process pydamage {
    tag "$name"

    label 'process_high'
    label 'process_ignore'

    publishDir "${params.outdir}/pydamage/$name", mode: 'copy'

    conda (params.enable_conda ? "bioconda::pydamage=0.60" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pydamage:0.62--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/pydamage:0.62--pyhdfd78af_0"
    }

    input:
        tuple val(name), path(bam)
    output:
        tuple val(name), path("*pydamage_results.csv"), emit: csv
        path "${name}/plots", optional: true, emit: plot
         path "*.version.txt", emit: version
    script:
        def software = getSoftwareName(task.process)
        output = name
        if (params.pydamage_plot) {
            plot = "--plot"
        } else {
            plot = ""
        }
        """
        samtools index $bam
        pydamage analyze --force -p ${task.cpus} -w ${params.wlen} $plot $bam
        pydamage filter pydamage_results/pydamage_results.csv
        mv pydamage_results ${name}
        mv ${name}/pydamage_results.csv ${name}.pydamage_results.csv
        mv ${name}/pydamage_filtered_results.csv ${name}.filtered_pydamage_results.csv
        echo \$(pydamage --version 2>&1) | sed -e 's/pydamage, version //g' > ${software}.version.txt
        """
}
