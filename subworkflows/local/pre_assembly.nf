include { adapterremoval } from "$baseDir/modules/local/adapterremoval.nf" params(params)
include { fastp } from "$baseDir/modules/local/fastp.nf" params(params)

workflow PRE_ASSEMBLY {
    take:
        reads
        adapter_list

    main:
        adapterremoval(reads, adapter_list)
        fastp(adapterremoval.out.trimmed_reads)
    emit:
        trimmed_reads = fastp.out.trimmed_reads
        fastp_logs = fastp.out.settings
        adapater_removal_logs = adapterremoval.out.settings
}
