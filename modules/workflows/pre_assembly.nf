include adapterremoval from "$baseDir/modules/tools/adapterremoval/main.nf" params(params)
include fastp from "$baseDir/modules/tools/fastp/main.nf" params(params)
include fastqc from "$baseDir/modules/tools/fastqc/main.nf" params(params)

workflow PRE_ASSEMBLY {
    take:
        reads
        adapter_list

    main:
        adapterremoval(reads, adapter_list)
        fastp(adapterremoval.out.trimmed_reads)
        fastqc(adapterremoval.out.trimmed_reads)
    emit:
        trimmed_reads = fastp.out.trimmed_reads
        fastqc_logs = fastqc.out
        fastp_logs = fastp.out.settings
        adapater_removal_logs = adapterremoval.out.settings
}