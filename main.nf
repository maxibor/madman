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
summary['assembly tool'] = params.assembly_tool 
summary['minlen'] = params.minlen
summary['minread'] = params.minread
summary['results'] = params.results
summary['Config Profile'] = workflow.profile
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
log.info summary.collect { k,v -> "${k.padRight(25)}: $v" }.join("\n")
log.info "----------------------------------------------------------------"


include adapterremoval from "$baseDir/modules/adapterremoval/main.nf" params(params)
include align_reads_to_contigs from "$baseDir/modules/bowtie2/main.nf" params(params)
include {damageprofiler as damageprofiler_pre; damageprofiler as damageprofiler_post} from "$baseDir/modules/damageprofiler/main.nf" params(params)
include fastp from "$baseDir/modules/fastp/main.nf" params(params)
include fastqc from "$baseDir/modules/fastqc/main.nf" params(params)
include filter_contigs_length from "$baseDir/modules/filter_contigs_length/main.nf" params(params)
include filter_contigs_damage from "$baseDir/modules/filter_contigs_damage/main.nf" params(params)
include megahit from "$baseDir/modules/megahit/main.nf" params(params)
include {metaspades; biospades} from "$baseDir/modules/metaspades/main.nf" params(params)
include prokka from "$baseDir/modules/prokka/main.nf" params(params)
include pydamage from "$baseDir/modules/pydamage/main.nf" params(params)
include {quast as quast_pre ; quast as quast_post} from "$baseDir/modules/quast/main.nf" params(params)
include multiqc from "$baseDir/modules/multiqc/main.nf" params(params)


workflow {
    adapterremoval(ch_reads, ch_adapter_list)
    fastp(adapterremoval.out.trimmed_reads)
    fastqc(adapterremoval.out.trimmed_reads)
    megahit(fastp.out.trimmed_reads)
    quast_pre(megahit.out.contigs, "pre")
    filter_contigs_length(megahit.out.contigs)
    align_reads_to_contigs(filter_contigs_length.out.join(fastp.out.trimmed_reads))
    damageprofiler_pre(megahit.out.contigs.join(align_reads_to_contigs.out), "pre")
    pydamage(align_reads_to_contigs.out)
    filter_contigs_damage(pydamage.out.csv.join(megahit.out.contigs))
    quast_post(filter_contigs_damage.out.fasta, "post")
    damageprofiler_post(filter_contigs_damage.out.fasta.join(align_reads_to_contigs.out),"post")
    prokka(filter_contigs_damage.out.fasta)
    multiqc(adapterremoval.out.settings, 
            fastqc.out.collect().ifEmpty([]), 
            quast_pre.out.collect().ifEmpty([]),
            quast_post.out.collect().ifEmpty([]),
            prokka.out.collect().ifEmpty([]),
            fastp.out.settings.collect().ifEmpty([]),
            damageprofiler_pre.out.collect().ifEmpty([]),
            damageprofiler_post.out.collect().ifEmpty([]),
            ch_multiqc_config)

}

// // Channel duplication and tuple index renaming to execute megahit and/or metaspades

// ch_trimmed_reads_mapping_pre.into {ch_trimmed_reads_mapping_megahit_pre; ch_trimmed_reads_mapping_metaspades_pre}

// ch_trimmed_reads_mapping_megahit_pre
//     .map {it -> [it[0]+"_megahit",[file(it[1][0]),file(it[1][1])]]}
//     .set {ch_trimmed_reads_mapping_megahit}
// ch_trimmed_reads_mapping_metaspades_pre
//     .map {it -> [it[0]+"_metaspades",[file(it[1][0]),file(it[1][1])]]}
//     .set {ch_trimmed_reads_mapping_metaspades}

// ch_trimmed_reads_mapping_megahit.mix(ch_trimmed_reads_mapping_metaspades).set{ch_trimmed_reads_mapping}


// if (params.assembly_tool.toString().contains('megahit') && params.assembly_tool.toString().contains('metaspades')) {
//     ch_trimmed_reads_assembly.into {ch_trimmed_reads_assembly_megahit; ch_trimmed_reads_assembly_metaspades}

// } else if (params.assembly_tool == 'metaspades') {
//     ch_trimmed_reads_assembly
//         .set {ch_trimmed_reads_assembly_metaspades}
//   ch_trimmed_reads_assembly_megahit = Channel.empty()
// } else if (params.assembly_tool == 'megahit') {
//     ch_trimmed_reads_assembly
//         .set {ch_trimmed_reads_assembly_megahit}
//     ch_trimmed_reads_assembly_metaspades = Channel.empty()
// }


// ch_megahit_filter.mix(ch_metaspades_filter).set{ch_contigs_filter}
// ch_megahit_quast.mix(ch_metaspades_quast).set{ch_contigs_quast}
// ch_megahit_log.mix(ch_metaspades_log).set{ch_contigs_log}
