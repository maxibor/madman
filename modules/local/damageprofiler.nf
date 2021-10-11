// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process damageprofiler {
    tag "$name"

    label 'process_medium'
    label 'process_ignore'

    publishDir "${params.outdir}/damageProfiler/${name}_${step}", mode: 'copy'

    conda (params.enable_conda ? "bioconda::damageprofiler=1.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/damageprofiler:1.1--hdfd78af_2"
    } else {
        container "quay.io/biocontainers/damageprofiler:1.1--hdfd78af_2"
    }

    input:
        tuple val(name), path(contig), path(bam)
        val(step)
    output:
        path("*.dmgprof_${step}.json")
        path  "*.version.txt", emit: version
    script:
        def software = getSoftwareName(task.process)
        outfile = name+".dmgprof_${step}.json"
        maxmem = task.memory.toGiga()
        """
        damageprofiler -Xmx${maxmem}g -i $bam -r $contig -o tmp
        mv tmp/${name}.sorted/dmgprof.json $outfile
        echo \$(damageprofiler -v) | sed 's/^DamageProfiler v//' > ${software}.version.txt
        """
}
