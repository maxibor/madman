process quast {
    tag "$name"

    label 'process_medium' 

    publishDir "${params.outdir}/quast/${outdir}", mode: 'copy'

    input:
        tuple val(name), path(contigs)
        val(step)
    output:
        path("*_quast_${step}")
    script:
        outdir = name+"_quast_${step}"
        """
        quast -o $outdir -t ${task.cpus} $contigs
        """
}