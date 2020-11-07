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
        file("*.pdf")
    script:
        outfile = name+".dmgprof_${step}.json"
        ref = ${contig.baseName()}
        maxmem = task.memory.toGiga()
        """
        damageprofiler -Xmx${maxmem}g -i $bam -r $contig -s ${ref} -o tmp
        mv tmp/${name}.sorted/dmgprof.json $outfile
        """
}
