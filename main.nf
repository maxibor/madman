#!/usr/bin/env nextflow

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
      --paired_end                      Specifies if reads are paired-end (true | false). Default: ${params.paired_end}
      --complexity_filter_poly_g_min    Length of poly-g min for clipping to be performed. Default: ${params.complexity_filter_poly_g_min}
      --minlen                          Minimum contig length to retain. Default:  ${params.minlen}
      --minread                         Minimum number of reads aligned to contig to consider contig. Default: ${params.minread}
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
    .fromFilePairs( params.reads, size: params.paired_end ? 2 : 1 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\n" }
	.set {ch_reads_to_trim}

ch_multiqc_config = file(params.multiqc_config, checkIfExists: true)

log.info "================================================================"
def summary = [:]
summary['Reads'] = params.reads
summary['phred'] = params.phred
summary['paired_end'] = params.paired_end
summary['assembly tool'] = params.assembly_tool 
summary['minlen'] = params.minlen
summary['minread'] = params.minread
summary['results'] = params.results
summary['Config Profile'] = workflow.profile
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
log.info summary.collect { k,v -> "${k.padRight(25)}: $v" }.join("\n")
log.info "----------------------------------------------------------------"


process AdapterRemoval {
    tag "$name"

    label 'expresso'

    input:
        set val(name), file(reads) from ch_reads_to_trim

    output:
        set val(name), file('*.trimmed.fastq') into ch_trimmed_reads_fastqc, ch_trimmed_reads_fastp
        file("*.settings") into ch_adapter_removal_results

    script:
        out1 = name+".pair1.trimmed.fastq"
        out2 = name+".pair2.trimmed.fastq"
        se_out = name+".trimmed.fastq"
        settings = name+".settings"
        if (params.paired_end){
            """
            AdapterRemoval --basename $name --file1 ${reads[0]} --file2 ${reads[1]} --trimns --trimqualities --minquality 20 --minlength 30 --output1 $out1 --output2 $out2 --threads ${task.cpus} --qualitybase ${params.phred} --settings $settings
            """
        } else {
            """
            AdapterRemoval --basename $name --file1 ${reads[0]} --trimns --trimqualities --minquality 20 --minlength 30 --output1 $se_out --threads ${task.cpus} --qualitybase ${params.phred} --settings $settings
            """
        }    
}

process fastp {
    tag "$name"

    label 'expresso'

    input:
        set val(name), file(reads) from ch_trimmed_reads_fastp

    output:
        set val(name), file("*.fq.gz") into ch_trimmed_reads_assembly, ch_trimmed_reads_mapping_pre
        file("*.json") into ch_fastp_for_multiqc

    script:
    if (params.paired_end) {
        out1 = name+".pair1.trimmed.fq.gz"
        out2 = name+".pair2.trimmed.fq.gz"
        se_out = name+".trimmed.fq.gz"
    """
    fastp --in1 ${reads[0]} --in2 ${reads[1]} --out1 $out1 --out2 $out2 -A -g --poly_g_min_len "${params.complexity_filter_poly_g_min}" -Q -L -w ${task.cpus} --json ${name}_fastp.json 
    """
    } else {
    """
    fastp --in1 ${reads[0]} --out1 $se_out -A -g --poly_g_min_len "${params.complexity_filter_poly_g_min}" -Q -L -w ${task.cpus} --json ${name}_fastp.json 
    """
    }
}


process fastqc {
    tag "$name"

    label 'intenso'

    input:
        set val(name), file(reads) from ch_trimmed_reads_fastqc
    output:
        file '*_fastqc.{zip,html}' into ch_fastqc_results
    script:
        """
        fastqc -t ${task.cpus} -q $reads
        """
}

// Channel duplication and tuple index renaming to execute megahit and/or metaspades

ch_trimmed_reads_mapping_pre.into {ch_trimmed_reads_mapping_megahit_pre; ch_trimmed_reads_mapping_metaspades_pre}

ch_trimmed_reads_mapping_megahit_pre
    .map {it -> [it[0]+"_megahit",[file(it[1][0]),file(it[1][1])]]}
    .set {ch_trimmed_reads_mapping_megahit}
ch_trimmed_reads_mapping_metaspades_pre
    .map {it -> [it[0]+"_metaspades",[file(it[1][0]),file(it[1][1])]]}
    .set {ch_trimmed_reads_mapping_metaspades}

ch_trimmed_reads_mapping_megahit.mix(ch_trimmed_reads_mapping_metaspades).set{ch_trimmed_reads_mapping}


if (params.assembly_tool.toString().contains('megahit') && params.assembly_tool.toString().contains('metaspades')) {
    ch_trimmed_reads_assembly.into {ch_trimmed_reads_assembly_megahit; ch_trimmed_reads_assembly_metaspades}

} else if (params.assembly_tool == 'metaspades') {
    ch_trimmed_reads_assembly
        .set {ch_trimmed_reads_assembly_metaspades}
  ch_trimmed_reads_assembly_megahit = Channel.empty()
} else if (params.assembly_tool == 'megahit') {
    ch_trimmed_reads_assembly
        .set {ch_trimmed_reads_assembly_megahit}
    ch_trimmed_reads_assembly_metaspades = Channel.empty()
}

process megahit {
    tag "$name"

    label 'bigmem'

    errorStrategy 'ignore'

    when:
        params.assembly_tool.toString().contains('megahit')

    publishDir "${params.results}/assembly/megahit/${name}", mode: 'copy'

    input:
        set val(name), file(reads) from ch_trimmed_reads_assembly_megahit

    output:
        set val(name2), file("megahit_out/*.megahit_contigs.fa") into ch_megahit_filter, ch_megahit_quast
        set val(name2), file("megahit_out/*.log") into ch_megahit_log
    script:
        mem = task.memory.toBytes()
        name2 = name+"_megahit"
        """
        megahit -1 ${reads[0]} -2 ${reads[1]} -t ${task.cpus} -m $mem --out-prefix $name
        mv megahit_out/${name}.contigs.fa megahit_out/${name}.megahit_contigs.fa
        """       
}

process metaspades {
    tag "$name"

    label 'bigmem'

    errorStrategy 'ignore'

    when: 
        params.assembly_tool.toString().contains('metaspades')

    publishDir "${params.results}/assembly/metaspades", mode: 'copy'

    input:
        set val(name), file(reads) from ch_trimmed_reads_assembly_metaspades

    output:
        set val(name2), file("*.metaspades_contigs.fa") into ch_metaspades_filter, ch_metaspades_quast
        set val(name2), file("$name") into ch_metaspades_log
    script:
        mem = task.memory.toGiga()
        name2 = name+"_metaspades"
        """
        spades.py --meta \
                  -1 ${reads[0]} \
                  -2 ${reads[1]} \
                  -t ${task.cpus} \
                  -m $mem \
                  --phred-offset ${params.phred} \
                  -o $name
        cp $name/contigs.fasta ${name}.metaspades_contigs.fa
        """ 
}

ch_megahit_filter.mix(ch_metaspades_filter).set{ch_contigs_filter}
ch_megahit_quast.mix(ch_metaspades_quast).set{ch_contigs_quast}
ch_megahit_log.mix(ch_metaspades_log).set{ch_contigs_log}

process quast_pre {
    tag "$name"

    label 'intenso'

    errorStrategy 'ignore'

    publishDir "${params.results}/quast/${name}", mode: 'copy'

    input:
        set val(name), file(contigs) from ch_contigs_quast
    output:
        file("*_quast_pre") into ch_quast_results_pre
    script:
        outdir = name+"_quast_pre"
        """
        quast -o $outdir -t ${task.cpus} $contigs
        """
}

process filter_contigs_size {
    tag "$name"

    label 'intenso'

    errorStrategy 'ignore'

    publishDir "${params.results}/fasta_filter/${name}", mode: 'copy'

    input:
        set val(name), file(fasta) from ch_contigs_filter
    output:
        set val(name), file("*.size_filtered.fa") into ch_contigs_filter_size, ch_contigs_filter_ancient, ch_contigs_dp
    script:
        filter_out = name+".size_filtered.fa"
        """
        filter_contigs_length.py -min ${params.minlen} -p ${task.cpus} $fasta -o $filter_out
        """
}

process align_reads_to_contigs {
    tag "$name"

    label 'intenso'

    errorStrategy 'ignore'

    publishDir "${params.results}/alignment/${name}", mode: 'copy'

    input:
        set val(name), file(contigs), file(reads) from ch_contigs_filter_size.join(ch_trimmed_reads_mapping)
    output:
        set val(name), file("*.sorted.bam") into (ch_alignment_to_dp_pre, ch_alignment_to_dp_post, ch_alignment_to_pydamage)
    script:
        outfile = name+".sorted.bam"
        if (params.paired_end) {
            """
            bowtie2-build --threads ${task.cpus} $contigs $name
            bowtie2 -x $name -1 ${reads[0]} -2 ${reads[1]} --threads ${task.cpus} | samtools view -S -b -F 4 - | samtools sort - > $outfile
            """
        } else {
            """
            bowtie2-build --threads ${task.cpus} $contigs $name
            bowtie2 -x $name -U $reads --threads ${task.cpus} | samtools view -S -b -F 4 - | samtools sort - > $outfile
            """
        }
}

process damageProfiler_pre {
    tag "$name"

    label 'expresso'

    errorStrategy 'ignore'

    publishDir "${params.results}/damageProfiler/${name}", mode: 'copy'

    input:
        set val(name), file(contig), file(bam) from ch_contigs_dp.join(ch_alignment_to_dp_pre)
    output:
        file("*_pre_filtering.dmgprof.json") into dmgProf_pre
    script:
        outfile = name+"_pre_filtering.dmgprof.json"
        """
        damageprofiler -i $bam -r $contig -o tmp
        mv tmp/${name}.sorted/dmgprof.json $outfile
        """
}

process pydamage {
    tag "$name"

    label 'intenso'

    publishDir "${params.results}/pydamage/", mode: 'copy'

    errorStrategy 'ignore'

    input:
        set val(name), file(bam) from ch_alignment_to_pydamage
    output:
        set val(name), file("*.csv") into ch_pydamage_stats
    script:
        output = name+".pydamage"
        """
        samtools index $bam
        pydamage -p ${task.cpus} -m ${params.minread} -o $output $bam
        """
}

process filter_contigs_damage {
    tag "$name"

    label 'ristretto'

    errorStrategy 'ignore'

    echo true

    publishDir "${params.results}/fasta_filter/${name}", mode: 'copy'

    input:
        set val(name), file(pydamage_csv), file(contigs) from ch_pydamage_stats.join(ch_contigs_filter_ancient)
    output:
        set val(name), file("*.ancient_filtered.fa") into ch_filtered_ancient_contigs_prokka, ch_filtered_ancient_contigs_quast, ch_ancient_contigs_dp
        set val(name), file("*_ancient_contigs.txt") into ch_ancient_contigs_list
    script:
        outfile = name + ".ancient_filtered.fa"
        ancient_contigs = name+"_ancient_contigs.txt"
        """
        filter_contigs_damage.py -d ${params.mindamage} -o $outfile $contigs $pydamage_csv
        awk -F "," '{if ((\$7 <= 0.05) && (\$5 >= ${params.mindamage})) { print \$1 }}' $pydamage_csv > $ancient_contigs
        """
}

process quast_post {
    tag "$name"

    label 'intenso'

    errorStrategy 'ignore'

    publishDir "${params.results}/quast/${name}", mode: 'copy'

    input:
        set val(name), file(contigs) from ch_filtered_ancient_contigs_quast
    output:
        file("*_quast_post") into ch_quast_results_post
    script:
        outdir = name+"_quast_post"
        """
        quast -o $outdir -t ${task.cpus} $contigs
        """
}

process damageProfiler_post {
    tag "$name"

    label 'expresso'

    errorStrategy 'ignore'

    publishDir "${params.results}/damageProfiler/${name}", mode: 'copy'

    input:
        set val(name), 
            file(contig), 
            file(bam), 
            file(acontigs) from ch_ancient_contigs_dp 
                                    .join(ch_alignment_to_dp_post)
                                    .join(ch_ancient_contigs_list)
    output:
        file("*_post_filtering.dmgprof.json") into dmgProf_post
    script:
        outfile = name+"_post_filtering.dmgprof.json"
        """
        damageprofiler -i ${bam} -r $contig -o tmp
        mv tmp/${name}.sorted/dmgprof.json $outfile
        """
}


// samtools view -@ ${task.cpus} $bam | grep -f $acontigs - > ${name}.ancient_contigs.bam -


process prokka {
    tag "$name"

    label 'intenso'

    errorStrategy 'ignore'

    publishDir "${params.results}/prokka", mode: 'copy'

    input:
        set val(name), file(contigs) from ch_filtered_ancient_contigs_prokka
    output:
        file("${name}") into ch_prokka_results
    script:
        """
        prokka --metagenome --cpus ${task.cpus} --outdir $name --prefix $name $contigs
        """
}

process multiqc {
    label 'ristretto'

    publishDir "${params.results}", mode: 'copy'

    input:
        file ('adapterRemoval/*') from ch_adapter_removal_results.collect().ifEmpty([])
        file ('fastqc/*') from ch_fastqc_results.collect().ifEmpty([])
        file ('quast_pre/*') from ch_quast_results_pre.collect().ifEmpty([])
        file ('quast_post/*') from ch_quast_results_post.collect().ifEmpty([])
        file ('prokka/*') from ch_prokka_results.collect().ifEmpty([])
        file ('fastp/*') from ch_fastp_for_multiqc.collect().ifEmpty([])
        file ('dmgprof_pre/*') from dmgProf_pre.collect().ifEmpty([])
        file ('dmgprof_post/*') from dmgProf_post.collect().ifEmpty([])
        file(multiqc_conf) from ch_multiqc_config
    output:
        file 'multiqc_report.html' into multiqc_report
    script:
        """
        multiqc -c $multiqc_conf .
        """
}


// process PMDtools {
//     tag "$name"

//     label 'intenso'

//     publishDir "${params.results}/pmdtools/", mode: 'copy', pattern: '*.fastq'

//     input:
//         set val(name), file(bam) from alignment_to_pmd
//     output:
//         file("*.fastq") into pmd_assemble
//     script:
//         fwd_out = name+"_pmd_R1.fastq"
//         rev_out = name+"_pmd_R2.fastq"
//         """
//         samtools view -h -F 4 $bam |\\
//         pmdtools -t ${params.pmdscore} --header |\\
//         samtools view -Sbh -@ ${task.cpus} - |\\
//         samtools fastq -1 $fwd_out -2 $rev_out -
//         """
// }

// process damageProfiler {
//     tag "$name"

//     label 'expresso'

//     errorStrategy 'ignore'

//     publishDir "${params.results}/damageProfiler/${name}", mode: 'copy'

//     input:
//         set val(name), file(bam), file(contig) from pmd_alignment_to_dp.join(pmd_contig_filtered_dp)
//     output:
//         file("*dmgprof.json") into dmgProf
//     script:
//         """
//         damageprofiler -i $bam -r $contig -o tmp
//         mv tmp/${name}.sorted/dmgprof.json ${name}.dmgprof.json
//         """
// }




// process miniKraken {
//     tag "$name"

//     conda 'bioconda::kraken2'

//     label 'intenso'

//     input:
//         set val(name), file(reads) from trimmed_reads

//     output:
//         set val(name), file('*.kraken.out') into kraken_out
//         set val(name), file('*.kreport') into kraken_report

//     script:
//         out = name+".kraken.out"
//         kreport = name+".kreport"
//         if (params.paired_end){
//             """
//             kraken2 --db ${params.krakendb} --threads ${task.cpus} --output $out --report $kreport --paired ${reads[0]} ${reads[1]}
//             """    
//         } else {
//             """
//             kraken2 --db ${params.krakendb} --threads ${task.cpus} --output $out --report $kreport ${reads[0]}
//             """
//         }
        
// }

// process kraken_parse {
//     tag "$name"

//     conda 'python=3.6'

//     label 'ristretto'

//     input:
//         set val(name), file(kraken_r) from kraken_report

//     output:
//         set val(name), file('*.kraken_parsed.csv') into kraken_parsed

//     script:
//         out = name+".kraken_parsed.csv"
//         """
//         kraken_parse.py -c ${params.minhit} -o $out $kraken_r
//         """    
// }

// process kraken_merge {

//     conda 'python=3.6 pandas numpy'

//     label 'ristretto'

//     publishDir "${params.results}", mode: 'copy'

//     input:
//         file(csv_count) from kraken_parsed.collect()

//     output:
//         file('kraken_otu_table.csv') into kraken_merged

//     script:
//         out = "kraken_otu_table.csv"
//         """
//         merge_kraken_res.py -o $out
//         """    
// }