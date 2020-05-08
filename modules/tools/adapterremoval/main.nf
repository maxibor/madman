process adapterremoval {
    tag "$name"

    label 'process_medium'

    input:
        tuple name, path(reads)
        path adapter_list

    output:
        tuple val(name), path('*.trimmed.fastq'), emit: trimmed_reads
        path "*.settings", emit: settings

    script:
        out1 = name+".pair1.trimmed.fastq"
        out2 = name+".pair2.trimmed.fastq"
        se_out = name+".trimmed.fastq"
        settings = name+".settings"
        if (! params.single_end){
            """
            AdapterRemoval --basename $name \
                           --file1 ${reads[0]} \
                           --file2 ${reads[1]} \
                           --trimns \
                           --trimqualities \
                           --minquality 20 \
                           --minlength 30 \
                           --output1 $out1 \
                           --output2 $out2 \
                           --threads ${task.cpus} \
                           --qualitybase ${params.phred} \
                           --adapter-list ${adapter_list} \
                           --settings $settings
            """
        } else {
            """
            AdapterRemoval --basename $name \
                           --file1 ${reads[0]} \
                           --trimns \
                           --trimqualities \
                           --minquality 20 \
                           --minlength 30 \
                           --output1 $se_out \
                           --threads ${task.cpus} \
                           --qualitybase ${params.phred} \
                           --adapter-list ${adapter_list} \
                           --settings $settings
            """
        }    
}