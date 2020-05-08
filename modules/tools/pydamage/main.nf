process pydamage {
    tag "$name"
    label 'process_medium'

    publishDir "${params.results}/pydamage/$name", mode: 'copy'

    input:
        tuple val(name), path(bam)
    output:
        tuple val(name), path("*.pydamage_results.csv"), emit: csv
        path "${name}/plots", emit: plot
    script:
        output = name
        // if (params.pydamage_plot) {
        //     plot = "--plot"
        // } else {
        //     plot = ""
        // }
        """
        samtools index $bam
        pydamage --force -p ${task.cpus} -m ${params.minread} -c ${params.coverage} -w ${params.wlen} --plot -o $output $bam
        mv ${name}/pydamage_results.csv ${name}.pydamage_results.csv
        """
}
