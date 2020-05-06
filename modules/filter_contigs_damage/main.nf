process filter_contigs_damage {
    tag "$name"

    label 'process_low'

    publishDir "${params.results}/fasta_filter/${name}", mode: 'copy'

    input:
        tuple val(name), path(pydamage_csv), path(contigs)
    output:
        tuple val(name), path("*.ancient_filtered.fa"), emit: fasta
        tuple val(name), path("*_ancient_contigs.txt"), emit: txt
    script:
        outfile = name + ".ancient_filtered.fa"
        ancient_contigs = name+"_ancient_contigs.txt"
        """
        filter_contigs_damage.py -d ${params.mindamage} -o $outfile $contigs $pydamage_csv
        awk -F "," '{if ((\$11 <= 0.05) && (\$8 >= ${params.mindamage})) { print \$1 }}' $pydamage_csv > $ancient_contigs
        """
}

