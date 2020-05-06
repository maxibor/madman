process fastqc {
    tag "FastQC - $name"

    label 'process_low'

    input:
        tuple val(name), path(reads)
    output:
        path '*_fastqc.{zip,html}'
    script:
        """
        fastqc -t ${task.cpus} -q $reads
        """
}