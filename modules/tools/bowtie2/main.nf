process align_reads_to_contigs {
    tag "$name"

    label 'process_high'
    label 'process_ignore'

    publishDir "${params.outdir}/alignment/${name}", mode: 'copy'

    conda (params.enable_conda ? "bioconda::bowtie2=2.4.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bowtie2:2.4.4--py39hbb4e92a_0"
    } else {
        container "quay.io/biocontainers/bowtie2:2.4.4--py39hbb4e92a_0"
    }

    input:
        tuple val(name), path(contigs), path(reads)
    output:
        tuple val(name), file("*.sorted.bam")
    script:
        outfile = name+".sorted.bam"
        if (!params.single_end) {
            """
            bowtie2-build --threads ${task.cpus} $contigs $name
            bowtie2 -x $name -1 ${reads[0]} -2 ${reads[1]} --very-sensitive -N 1 --threads ${task.cpus} | samtools view -S -b -F 4 - | samtools sort - > $outfile
            """
        } else {
            """
            bowtie2-build --threads ${task.cpus} $contigs $name
            bowtie2 -x $name -U $reads --very-sensitive -N 1 --threads ${task.cpus} | samtools view -S -b -F 4 - | samtools sort - > $outfile
            """
        }
}