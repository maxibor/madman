process damageprofiler {
    tag "$name"

    label 'process_medium'
    label 'process_ignore'

    publishDir "${params.outdir}/damageProfiler/${name}_${step}", mode: 'copy'

    input:
        tuple val(name), path(contig), path(bam)
        val(step)
    output:
        path("*.dmgprof_${step}.json")
    script:
        outfile = name+".dmgprof_${step}.json"
        maxmem = task.memory.toGiga()
        """
        damageprofiler -Xmx${maxmem}g -i $bam -r $contig -o tmp
        mv tmp/${name}.sorted/dmgprof.json $outfile
        """
}