process metaspades {
    tag "$name"

    label 'process_very_bigmem'

    publishDir "${params.outdir}/assembly/metaspades/${name}", mode: 'copy'

    input:
        tuple val(name), path(reads)
    output:
        tuple val(name), path("*_metaspades.contigs.fa"), emit: contigs
        tuple val(name), path("${name}/*"), emit: metaspades_logs
    script:
        mem = task.memory.toGiga()
        """
        spades.py --meta \
                  -1 ${reads[0]} \
                  -2 ${reads[1]} \
                  -t ${task.cpus} \
                  -m $mem \
                  --phred-offset ${params.phred} \
                  -o $name
        cp $name/contigs.fasta ${name}_metaspades.contigs.fa
        """ 
}

process biospades {
    tag "$name"

    label 'process_very_bigmem'
    
    publishDir "${params.outdir}/assembly/biospades/${name}", mode: 'copy'

    input:
        tuple val(name), path(reads)
    output:
        tuple val(name), path("*_biospades.contigs.fa"), emit: contigs
        tuple val(name), path("${name}/*"), emit: biospades_logs
    script:
        mem = task.memory.toGiga()
        """
        spades.py --meta \
                  --bio \
                  -1 ${reads[0]} \
                  -2 ${reads[1]} \
                  -t ${task.cpus} \
                  -m $mem \
                  --phred-offset ${params.phred} \
                  -o $name
        cp $name/gene_clusters.fasta ${name}_biospades.contigs.fa
        """ 
}