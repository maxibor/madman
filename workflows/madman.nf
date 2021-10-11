/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMadman.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )
include { PRE_ASSEMBLY } from "$baseDir/$baseDir/subworkflows/local/pre_assembly.nf" params(params)

if (params.modern){
    include {POST_ASSEMBLY_MODERN as POST_ASSEMBLY_MEGAHIT ;
             POST_ASSEMBLY_MODERN as POST_ASSEMBLY_BIOSPADES;
             POST_ASSEMBLY_MODERN as POST_ASSEMBLY_METASPADES
            } from "$baseDir/subworkflows/local/post_assembly.nf" params(params)
} else {
    include {POST_ASSEMBLY_ANCIENT as POST_ASSEMBLY_MEGAHIT ;
             POST_ASSEMBLY_ANCIENT as POST_ASSEMBLY_BIOSPADES;
             POST_ASSEMBLY_ANCIENT as POST_ASSEMBLY_METASPADES
            } from "$baseDir/subworkflows/local/post_assembly.nf" params(params)
}
/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC  } from '../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Setting up necessaru Channels
ch_adapter_list = file(params.adapter_list, checkIfExists: true)


// Info required for completion email and summary
def multiqc_report = []

workflow MADMAN {

    ch_software_versions = Channel.empty()
    PRE_ASSEMBLY(ch_reads, ch_adapter_list)
    quast_pre_ch = Channel.empty()
    quast_post_ch = Channel.empty()
    damageprofiler_pre_ch = Channel.empty()
    damageprofiler_post_ch = Channel.empty()
    prokka_ch = Channel.empty()
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )

    if (params.megahit) {
        megahit(PRE_ASSEMBLY.out.trimmed_reads)
        POST_ASSEMBLY_MEGAHIT(megahit.out.contigs, PRE_ASSEMBLY.out.trimmed_reads, "megahit")
        quast_pre_ch.mix(POST_ASSEMBLY_MEGAHIT.out.quast_pre).set{quast_pre_ch}
        if (! params.modern) {
            quast_post_ch.mix(POST_ASSEMBLY_MEGAHIT.out.quast_post).set{quast_post_ch}
            damageprofiler_pre_ch.mix(POST_ASSEMBLY_MEGAHIT.out.damageprofiler_pre).set{damageprofiler_pre_ch}
            damageprofiler_post_ch.mix(POST_ASSEMBLY_MEGAHIT.out.damageprofiler_post).set{damageprofiler_post_ch}
        }
        prokka_ch.mix(POST_ASSEMBLY_MEGAHIT.out.prokka).set{prokka_ch}
    }

    if (params.metaspades) {
        metaspades(PRE_ASSEMBLY.out.trimmed_reads)
        POST_ASSEMBLY_METASPADES(metaspades.out.contigs, PRE_ASSEMBLY.out.trimmed_reads, "metaspades")
        quast_pre_ch.mix(POST_ASSEMBLY_METASPADES.out.quast_pre).set{quast_pre_ch}
        if (! params.modern){
            quast_post_ch.mix(POST_ASSEMBLY_METASPADES.out.quast_post).set{quast_post_ch}
            damageprofiler_pre_ch.mix(POST_ASSEMBLY_METASPADES.out.damageprofiler_pre).set{damageprofiler_pre_ch}
            damageprofiler_post_ch.mix(POST_ASSEMBLY_METASPADES.out.damageprofiler_post).set{damageprofiler_post_ch}
        }
        prokka_ch.mix(POST_ASSEMBLY_METASPADES.out.prokka).set{prokka_ch}
    }

    if (params.biospades) {
        biospades(PRE_ASSEMBLY.out.trimmed_reads)
        POST_ASSEMBLY_BIOSPADES(biospades.out.contigs, PRE_ASSEMBLY.out.trimmed_reads, "biospades")
        quast_pre_ch.mix(POST_ASSEMBLY_BIOSPADES.out.quast_pre).set{quast_pre_ch}
        if (! params.modern){
            quast_post_ch.mix(POST_ASSEMBLY_BIOSPADES.out.quast_post).set{quast_post_ch}
            damageprofiler_pre_ch.mix(POST_ASSEMBLY_BIOSPADES.out.damageprofiler_pre).set{damageprofiler_pre_ch}
            damageprofiler_post_ch.mix(POST_ASSEMBLY_BIOSPADES.out.damageprofiler_post).set{damageprofiler_post_ch}
        }
        prokka_ch.mix(POST_ASSEMBLY_BIOSPADES.out.prokka).set{prokka_ch}
    }



    ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))

    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMadman.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
