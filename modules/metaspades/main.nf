process metaspades {
    tag "metaSPAdes - $name"

    label 'process_bigmem'

    publishDir "${params.results}/assembly/metaspades/${name}", mode: 'copy'

    input:
        tuple val(name), path(reads)
    output:
        tuple val(name2), path("*.metaspades_contigs.fa"), emit: contigs
        set val(name2), path("${name}/*"), emit: all
    script:
        mem = task.memory.toGiga()
        name2 = name+"_metaspades"
        """
        spades.py --meta \
                  -1 ${reads[0]} \
                  -2 ${reads[1]} \
                  -t ${task.cpus} \
                  -m $mem \
                  --phred-offset ${params.phred} \
                  -o $name
        cp $name/contigs.fasta ${name}.metaspades_contigs.fa
        """ 
}

process biospades {
    tag "bioSPAdes - $name"

    label 'process_bigmem'

    publishDir "${params.results}/assembly/metaspades/${name}", mode: 'copy'

    input:
        tuple val(name), path(reads)
    output:
        tuple val(name2), path("*.biospades_contigs.fa"), emit: contigs
        set val(name2), path("${name}/*"), emit: all
    script:
        mem = task.memory.toGiga()
        name2 = name+"_biospades"
        """
        spades.py --meta \
                  --bio \
                  -1 ${reads[0]} \
                  -2 ${reads[1]} \
                  -t ${task.cpus} \
                  -m $mem \
                  --phred-offset ${params.phred} \
                  -o $name
        cp $name/gene_clusters.fasta ${name}.biospades_contigs.fa
        """ 
}