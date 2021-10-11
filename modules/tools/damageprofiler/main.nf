include { getSoftwareName } from "$baseDir/modules/tools/nf_core_utils/main.nf"

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
        path("*.dmgprof_${step}.json"), emit: json
        path  "*.version.txt", emit: version
    script:
        outfile = name+".dmgprof_${step}.json"
        maxmem = task.memory.toGiga()
        def software = getSoftwareName(task.process)
        """
        damageprofiler -Xmx${maxmem}g -i $bam -r $contig -o tmp
        mv tmp/${name}.sorted/dmgprof.json $outfile
        echo \$(damageprofiler -v) | sed 's/^DamageProfiler v//' > ${software}.version.txt
        """
}