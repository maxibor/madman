#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

def helpMessage() {
    log.info"""
     MADMAN: Metagenomic Assembly of Ancient DaMaged reads with Nextflow
     Homepage: https://github.com/maxibor/madman
     Author: Maxime Borry <borry@shh.mpg.de>
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run maxibor/madman --reads '/path/to/paired_end_reads_*.{1,2}.fastq.gz'
    Mandatory arguments:
      --reads                           Path to input data (must be surrounded with quotes)

    Settings:
      --phred                           Specifies the fastq quality encoding (33 | 64). Default: ${params.phred}
      --single_end                      To specify if reads are single-end.
      --adapter_list                    List of sequencing adapters to trim. Default: ${params.adapter_list}
      --complexity_filter_poly_g_min    Length of poly-g min for clipping to be performed. Default: ${params.complexity_filter_poly_g_min}
      --minlen                          Minimum contig length to retain. Default:  ${params.minlen}
      --minread                         Minimum number of reads aligned to contig to consider contig. Default: ${params.minread}
      --coverage                        Minimum coverage to consider contig. Default: ${params.coverage}
      --wlen                            Window length from 5' end to consider for damage estimation. Default: ${params.wlen}
      --mindamage                       Mimimum amount of CtoT damage on the 5' end of the read. Default: ${params.mindamage}
      --assembly_tool                   Choose de novo assembly tool, seperated by ',' (megahit | metaspades). Default: ${params.assembly_tool}
      

    Options:
      --results                         The output directory where the results will be saved. Default: ${params.results}
      --help  --h                       Shows this help page
    """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

Channel
    .fromFilePairs( params.reads, size: params.single_end ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\n" }
	.set {ch_reads}

ch_multiqc_config = file(params.multiqc_config, checkIfExists: true)

ch_adapter_list = file(params.adapter_list, checkIfExists: true)

log.info "================================================================"
def summary = [:]
summary['Reads'] = params.reads
summary['phred'] = params.phred
summary['single_end'] = params.single_end
summary['Run Megahit'] = params.megahit
summary['Run MetaSpades'] = params.metaspades
summary['Run Biosynthetic Spades'] = params.biospades
summary['minlen'] = params.minlen
summary['minread'] = params.minread
summary['results'] = params.results
summary['Config Profile'] = workflow.profile
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
log.info summary.collect { k,v -> "${k.padRight(25)}: $v" }.join("\n")
log.info "----------------------------------------------------------------"

include megahit from "$baseDir/modules/tools/megahit/main.nf" params(params)
include {metaspades ; biospades} from "$baseDir/modules/tools/metaspades/main.nf" params(params)
include multiqc from "$baseDir/modules/tools/multiqc/main.nf" params(params)
include PRE_ASSEMBLY from "$baseDir/modules/workflows/pre_assembly.nf" params(params)
include {POST_ASSEMBLY as POST_ASSEMBLY_MEGAHIT ; POST_ASSEMBLY as POST_ASSEMBLY_BIOSPADES; POST_ASSEMBLY as POST_ASSEMBLY_METASPADES} from "$baseDir/modules/workflows/post_assembly.nf" params(params)


workflow {
    PRE_ASSEMBLY(ch_reads, ch_adapter_list)
    quast_pre_ch = Channel.empty()
    quast_post_ch = Channel.empty()
    damageprofiler_pre_ch = Channel.empty()
    damageprofiler_post_ch = Channel.empty()
    prokka_ch = Channel.empty()

    if (params.megahit) {
        megahit(PRE_ASSEMBLY.out.trimmed_reads)
        POST_ASSEMBLY_MEGAHIT(megahit.out.contigs, PRE_ASSEMBLY.out.trimmed_reads, "megahit")
        quast_pre_ch.mix(POST_ASSEMBLY_MEGAHIT.out.quast_pre)
        quast_post_ch.mix(POST_ASSEMBLY_MEGAHIT.out.quast_post)
        damageprofiler_pre_ch.mix(POST_ASSEMBLY_MEGAHIT.out.damageprofiler_pre)
        damageprofiler_post_ch.mix(POST_ASSEMBLY_MEGAHIT.out.damageprofiler_post)
        prokka_ch.mix(POST_ASSEMBLY_MEGAHIT.out.prokka)
    }

    if (params.biospades) {
        biospades(PRE_ASSEMBLY.out.trimmed_reads)
        POST_ASSEMBLY_BIOSPADES(biospades.out.contigs, PRE_ASSEMBLY.out.trimmed_reads, "biospades")
        quast_pre_ch.mix(POST_ASSEMBLY_BIOSPADES.out.quast_pre)
        quast_post_ch.mix(POST_ASSEMBLY_BIOSPADES.out.quast_post)
        damageprofiler_pre_ch.mix(POST_ASSEMBLY_BIOSPADES.out.damageprofiler_pre)
        damageprofiler_post_ch.mix(POST_ASSEMBLY_BIOSPADES.out.damageprofiler_post)
        prokka_ch.mix(POST_ASSEMBLY_BIOSPADES.out.prokka)
    } 

    if (params.metaspades) {
        metaspades(PRE_ASSEMBLY.out.trimmed_reads)
        POST_ASSEMBLY_METASPADES(metaspades.out.contigs, PRE_ASSEMBLY.out.trimmed_reads, "metaspades")
        quast_pre_ch.mix(POST_ASSEMBLY_METASPADES.out.quast_pre)
        quast_post_ch.mix(POST_ASSEMBLY_METASPADES.out.quast_post)
        damageprofiler_pre_ch.mix(POST_ASSEMBLY_METASPADES.out.damageprofiler_pre)
        damageprofiler_post_ch.mix(POST_ASSEMBLY_METASPADES.out.damageprofiler_post)
        prokka_ch.mix(POST_ASSEMBLY_METASPADES.out.prokka)
    }
    
    multiqc(PRE_ASSEMBLY.out.adapater_removal_logs.collect().ifEmpty([]), 
            PRE_ASSEMBLY.out.fastqc_logs.collect().ifEmpty([]), 
            PRE_ASSEMBLY.out.fastp_logs.collect().ifEmpty([]),
            quast_pre_ch.collect().ifEmpty([]),
            quast_post_ch.collect().ifEmpty([]),
            damageprofiler_pre_ch.collect().ifEmpty([]),
            damageprofiler_post_ch.collect().ifEmpty([]),
            prokka_ch.collect().ifEmpty([]),
            ch_multiqc_config)
}