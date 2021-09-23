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
        tuple val(name), path("*.pydamage_results.csv"), emit: csv
        path "${name}/plots", optional: true, emit: plot
    script:
        output = name
        if (params.pydamage_plot) {
            plot = "--plot"
        } else {
            plot = ""
        }
        """
        samtools index $bam
        pydamage -o $output analyze --force -p ${task.cpus} -w ${params.wlen} $plot $bam
        mv ${name}/pydamage_results.csv ${name}.pydamage_results.csv
        """
}
