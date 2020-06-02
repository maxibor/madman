process megahit {
    tag "MEGAHIT - $name"

    label 'process_bigmem'
    label 'process_mandatory'

    publishDir "${params.outdir}/assembly/megahit/${name}", mode: 'copy'

    input:
        tuple val(name), path(reads)

    output:
        tuple val(name), path("megahit_out/*.megahit_contigs.fa"), emit: contigs
        tuple val(name), path("megahit_out/*.log"), emit: log
    script:
        mem = task.memory.toBytes()
        name2 = name+"_megahit"
        """
        megahit -1 ${reads[0]} -2 ${reads[1]} -t ${task.cpus} -m $mem --out-prefix $name
        mv megahit_out/${name}.contigs.fa megahit_out/${name}.megahit_contigs.fa
        """       
}