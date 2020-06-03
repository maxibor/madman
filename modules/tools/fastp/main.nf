process fastp {
    tag "$name"

    label 'process_high'
    label 'process_mantory'
    input:
        tuple val(name), path(reads)

    output:
        tuple val(name), path("*.fq.gz"), emit: trimmed_reads
        path "*.json", emit: settings

    script:
        if (!params.single_end) {
            out1 = name+".pair1.trimmed.fq.gz"
            out2 = name+".pair2.trimmed.fq.gz"
            """
            fastp --in1 ${reads[0]} --in2 ${reads[1]} --out1 $out1 --out2 $out2 -A -g --poly_g_min_len "${params.complexity_filter_poly_g_min}" -Q -L -w ${task.cpus} --json ${name}.json 
            """
        } else {
            se_out = name+".trimmed.fq.gz"
            """
            fastp --in1 ${reads[0]} --out1 $se_out -A -g --poly_g_min_len "${params.complexity_filter_poly_g_min}" -Q -L -w ${task.cpus} --json ${name}.json 
            """
        }
}