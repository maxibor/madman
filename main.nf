#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/madman
========================================================================================
 nf-core/madman Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/madman
----------------------------------------------------------------------------------------
*/


nextflow.preview.dsl = 2

def helpMessage() {
    log.info nfcoreHeader()
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
      --results                         The output directory where the results will be saved. Default: ${params.outdir}
      --help  --h                       Shows this help page
    """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

Channel
    .fromFilePairs( params.reads, size: params.single_end ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\n" }
	.set {ch_reads}

ch_multiqc_config = file(params.multiqc_config, checkIfExists: true)
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)
ch_adapter_list = file(params.adapter_list, checkIfExists: true)

runName = workflow.runName

log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Reads'] = params.reads
summary['phred'] = params.phred
summary['single_end'] = params.single_end
summary['Run Megahit'] = params.megahit
summary['Run MetaSpades'] = params.metaspades
summary['Run Biosynthetic Spades'] = params.biospades
summary['Make Pydamage plots'] = params.pydamage_plot
summary['minlen'] = params.minlen
summary['minread'] = params.minread
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(25)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-madman-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/madman Workflow Summary'
    section_href: 'https://github.com/nf-core/madman'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

include megahit from "$baseDir/modules/tools/megahit/main.nf" params(params)
include {metaspades ; biospades} from "$baseDir/modules/tools/metaspades/main.nf" params(params)
include multiqc from "$baseDir/modules/tools/multiqc/main.nf" params(params)
include PRE_ASSEMBLY from "$baseDir/modules/workflows/pre_assembly.nf" params(params)
include {POST_ASSEMBLY as POST_ASSEMBLY_MEGAHIT ; POST_ASSEMBLY as POST_ASSEMBLY_BIOSPADES; POST_ASSEMBLY as POST_ASSEMBLY_METASPADES} from "$baseDir/modules/workflows/post_assembly.nf" params(params)
include {output_documentation ; get_software_versions} from "$baseDir/modules/tools/nf_core_utils/main.nf" params(params)

workflow {
    get_software_versions()
    PRE_ASSEMBLY(ch_reads, ch_adapter_list)
    quast_pre_ch = Channel.empty()
    quast_post_ch = Channel.empty()
    damageprofiler_pre_ch = Channel.empty()
    damageprofiler_post_ch = Channel.empty()
    prokka_ch = Channel.empty()

    if (params.megahit) {
        megahit(PRE_ASSEMBLY.out.trimmed_reads)
        POST_ASSEMBLY_MEGAHIT(megahit.out.contigs, PRE_ASSEMBLY.out.trimmed_reads, "megahit")
        quast_pre_ch.mix(POST_ASSEMBLY_MEGAHIT.out.quast_pre).set{quast_pre_ch}
        quast_post_ch.mix(POST_ASSEMBLY_MEGAHIT.out.quast_post).set{quast_post_ch}
        damageprofiler_pre_ch.mix(POST_ASSEMBLY_MEGAHIT.out.damageprofiler_pre).set{damageprofiler_pre_ch}
        damageprofiler_post_ch.mix(POST_ASSEMBLY_MEGAHIT.out.damageprofiler_post).set{damageprofiler_post_ch}
        prokka_ch.mix(POST_ASSEMBLY_MEGAHIT.out.prokka).set{prokka_ch}
    }

    if (params.biospades) {
        biospades(PRE_ASSEMBLY.out.trimmed_reads)
        POST_ASSEMBLY_BIOSPADES(biospades.out.contigs, PRE_ASSEMBLY.out.trimmed_reads, "biospades")
        quast_pre_ch.mix(POST_ASSEMBLY_BIOSPADES.out.quast_pre).set{quast_pre_ch}
        quast_post_ch.mix(POST_ASSEMBLY_BIOSPADES.out.quast_post).set{quast_post_ch}
        damageprofiler_pre_ch.mix(POST_ASSEMBLY_BIOSPADES.out.damageprofiler_pre).set{damageprofiler_pre_ch}
        damageprofiler_post_ch.mix(POST_ASSEMBLY_BIOSPADES.out.damageprofiler_post).set{damageprofiler_post_ch}
        prokka_ch.mix(POST_ASSEMBLY_BIOSPADES.out.prokka).set{prokka_ch}
    } 

    if (params.metaspades) {
        metaspades(PRE_ASSEMBLY.out.trimmed_reads)
        POST_ASSEMBLY_METASPADES(metaspades.out.contigs, PRE_ASSEMBLY.out.trimmed_reads, "metaspades")
        quast_pre_ch.mix(POST_ASSEMBLY_METASPADES.out.quast_pre).set{quast_pre_ch}
        quast_post_ch.mix(POST_ASSEMBLY_METASPADES.out.quast_post).set{quast_post_ch}
        damageprofiler_pre_ch.mix(POST_ASSEMBLY_METASPADES.out.damageprofiler_pre).set{damageprofiler_pre_ch}
        damageprofiler_post_ch.mix(POST_ASSEMBLY_METASPADES.out.damageprofiler_post).set{damageprofiler_post_ch}
        prokka_ch.mix(POST_ASSEMBLY_METASPADES.out.prokka).set{prokka_ch}
    }
    
    multiqc(PRE_ASSEMBLY.out.adapater_removal_logs.collect().ifEmpty([]), 
            PRE_ASSEMBLY.out.fastqc_logs.collect().ifEmpty([]), 
            PRE_ASSEMBLY.out.fastp_logs.collect().ifEmpty([]),
            quast_pre_ch.collect().ifEmpty([]),
            quast_post_ch.collect().ifEmpty([]),
            damageprofiler_pre_ch.collect().ifEmpty([]),
            damageprofiler_post_ch.collect().ifEmpty([]),
            prokka_ch.collect().ifEmpty([]),
            get_software_versions.out.yaml,
            ch_multiqc_config,
            runName,
            ch_workflow_summary)

    output_documentation(ch_output_docs)
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/madman] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/madman] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // TODO nf-core: If not using MultiQC, strip out this code (including params.max_multiqc_email_size)
    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/madman] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/madman] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/madman] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            [ 'mail', '-s', subject, email_address ].execute() << email_txt
            log.info "[nf-core/madman] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/madman]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/madman]${c_red} Pipeline completed with errors${c_reset}-"
    }

}


def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/madman v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
