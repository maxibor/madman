process damageprofiler {
    tag "$name"

    label 'process_high'

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
        """
        damageprofiler -i $bam -r $contig -o tmp -s ${ref}
        mv tmp/${name}.sorted/dmgprof.json $outfile
        """
}
