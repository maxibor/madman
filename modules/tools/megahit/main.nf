process megahit {
    tag "MEGAHIT - $name"

    label 'process_high'
    label 'process_mandatory'

    publishDir "${params.outdir}/assembly/megahit/${name}", mode: 'copy'

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