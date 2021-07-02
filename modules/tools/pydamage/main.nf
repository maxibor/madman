process pydamage {
    tag "$name"
    
    label 'process_high'
    label 'process_ignore'

    publishDir "${params.outdir}/pydamage/$name", mode: 'copy'

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
