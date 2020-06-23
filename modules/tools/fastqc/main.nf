process fastqc {
    tag "$name"

    label 'process_low'
    label 'process_mandatory'
    
    input:
        tuple val(name), path(reads)
    output:
        path '*_fastqc.{zip,html}'
    script:
        """
        fastqc -t ${task.cpus} -q $reads
        """
}