process filter_contigs_length {
    tag "$name"

    label 'intenso'

    publishDir "${params.outdir}/fasta_filter/${name}", mode: 'copy'

    input:
        tuple val(name), path(fasta)
    output:
        tuple val(name), path("*.size_filtered.fa")
    script:
        filter_out = name+".size_filtered.fa"
        """
        filter_contigs_length.py -min ${params.minlen} -p ${task.cpus} $fasta -o $filter_out
        """
}