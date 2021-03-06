process prokka {
    tag "$name"

    label 'process_high'
    label 'process_ignore'

    publishDir "${params.outdir}/prokka", mode: 'copy'

    input:
        tuple val(name), path(contigs)
    output:
        path("${name}")
    script:
        """
        prokka --metagenome --cpus ${task.cpus} --outdir $name --prefix $name $contigs
        """
}
