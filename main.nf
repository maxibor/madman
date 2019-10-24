#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
     megahit-nf: simple Megahit assembler Nextflow pipeline
     Homepage: https://github.com/maxibor/megahit-nf
     Author: Maxime Borry <borry@shh.mpg.de>
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run maxibor/megahit-nf --reads '/path/to/paired_end_reads_*.{1,2}.fastq.gz'
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)

    Settings:
      --phred                       Specifies the fastq quality encoding (33 | 64). Defaults to ${params.phred}
      --pairedEnd                   Specifies if reads are paired-end (true | false). Default = ${params.pairedEnd}
      --minlen                      Minimum contig length to retain. Default =  ${params.minlen}
      --pmd_minlen                     Minimum pmd contig length to retain. Default =  ${params.pmd_minlen}
      --pmdscore                    PMDTools score threshold. Default = ${params.pmdscore}

    Options:
      --results                     The output directory where the results will be saved. Defaults to ${params.results}
      --help  --h                   Shows this help page
    """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

Channel
    .fromFilePairs( params.reads, size: params.pairedEnd ? 2 : 1 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\n" }
	.set {reads_to_trim}

log.info "================================================================"
def summary = [:]
summary['Reads'] = params.reads
summary['minlen'] = params.minlen
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(25)}: $v" }.join("\n")
log.info "----------------------------------------------------------------"


process AdapterRemoval {
    tag "$name"

    label 'intenso'

    input:
        set val(name), file(reads) from reads_to_trim

    output:
        set val(name), file('*.trimmed.fastq') into trimmed_reads_assembly, trimmed_reads_mapping, trimmed_reads_fastqc, trimmed_reads_pmd_mapping
        file("*.settings") into adapter_removal_results
        val(name) into samp_name

    script:
        out1 = name+".pair1.trimmed.fastq"
        out2 = name+".pair2.trimmed.fastq"
        se_out = name+".trimmed.fastq"
        settings = name+".settings"
        if (params.pairedEnd){
            """
            AdapterRemoval --basename $name --file1 ${reads[0]} --file2 ${reads[1]} --trimns --trimqualities --minquality 20 --minlength 30 --output1 $out1 --output2 $out2 --threads ${task.cpus} --qualitybase ${params.phred} --settings $settings
            """
        } else {
            """
            AdapterRemoval --basename $name --file1 ${reads[0]} --trimns --trimqualities --minquality 20 --minlength 30 --output1 $se_out --threads ${task.cpus} --qualitybase ${params.phred} --settings $settings
            """
        }    
}

process fastqc {
    tag "$fastqc"

    label 'intenso'

    input:
        set val(name), file(reads) from trimmed_reads_fastqc
    output:
        file '*_fastqc.{zip,html}' into fastqc_results
    script:
        """
        fastqc -t ${task.cpus} -q $reads
        """
}

process megahit {
    tag "$name"

    label 'bigmem'

    publishDir "${params.results}/assembly/${name}", mode: 'copy'

    input:
        set val(name), file(reads) from trimmed_reads_assembly

    output:
        set val(name), file("megahit_out/*.contigs.fa") into contigs_filter
        set val(name), file("megahit_out/*.log") into megahit_log
    script:
        mem = task.memory.toBytes()
        """
        megahit -1 ${reads[0]} -2 ${reads[1]} -t ${task.cpus} -m $mem --out-prefix $name
        """       
}

process filter_fasta {
    tag "$name"

    label 'intenso'

    publishDir "${params.results}/fasta_filter/${name}", mode: 'copy'

    input:
        set val(name), file(fasta) from contigs_filter
    output:
        set val(name), file("*.filtered.fa") into filtered_fa
    script:
        filter_out = name+".filtered.fa"
        """
        filter_fasta_length.py -min ${params.minlen} -p ${task.cpus} $fasta -o $filter_out
        """
}

process split_fasta {
    tag "$name"

    label 'intenso'

    publishDir "${params.results}/fasta_filter/${name}", mode: 'copy'

    input:
        set val(name), file(fasta) from filtered_fa
    output:
        file("*.split.fa") into split_contigs
    script:
        """
        split_fasta.py -p ${task.cpus} $fasta
        """
}

split_contigs
    .flatten()
    .map { it -> tuple(it.baseName, it) }  
    .into {filtered_contigs_bt; filtered_contigs_dp}

trimmed_reads_mapping
    .map {it -> it}
    .set {trimmed_reads_mapping_ch}

process align_reads_to_contigs {
    tag "$name"

    label 'intenso'

    errorStrategy 'ignore'

    publishDir "${params.results}/alignment/${name}", mode: 'copy'

    input:
        set val(name), file(contigs), val(name2), file(reads) from filtered_contigs_bt.combine(trimmed_reads_mapping_ch)
    output:
        set val(name), file("*.sorted.bam") into (alignment_to_dp, alignment_to_pmd)
    script:
        outfile = name+".sorted.bam"
        if (params.pairedEnd) {
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

process PMDtools {
    tag "$name"

    label 'intenso'

    publishDir "${params.results}/pmdtools/", mode: 'copy', pattern: '*.fastq'

    input:
        set val(name), file(bam) from alignment_to_pmd
    output:
        file("*.fastq") into pmd_assemble
    script:
        fwd_out = name+"_pmd_R1.fastq"
        rev_out = name+"_pmd_R2.fastq"
        """
        samtools view -h -F 4 $bam |\\
        pmdtools -t ${params.pmdscore} --header |\\
        samtools view -Sbh -@ ${task.cpus} - |\\
        samtools fastq -1 $fwd_out -2 $rev_out -
        """
}

process megahit_pmd {
    tag "$name"

    label 'bigmem'

    publishDir "${params.results}/pmd_assembly", mode: 'copy'

    input:
        file(reads) from pmd_assemble.collect()
        val(name) from samp_name

    output:
        set val(name), file("megahit_out/*.contigs.fa") into pmd_contigs_filter, pmd_contigs_quast
        file("megahit_out/*.log") into pmd_megahit_log
    script:
        mem = task.memory.toBytes()
        """
        cat *_pmd_R1.fastq > forward.fq
        cat *_pmd_R2.fastq > reverse.fq
        megahit -1 forward.fq -2 reverse.fq -t ${task.cpus} -m $mem --out-prefix $name
        """       
}

process quast {
    tag "$name"

    label 'intenso'

    publishDir "${params.results}/quast", mode: 'copy'

    input:
        set val(name), file(contigs) from pmd_contigs_quast
    output:
        file("quast_result/report.tsv") into quast_multiqc
        file("quast_result/*") into quast_results
    script:
        """
        quast -o quast_result -t ${task.cpus} $contigs
        """
}

process filter_fasta_pmd {
    tag "$name"

    label 'intenso'

    publishDir "${params.results}/fasta_filter/${name}", mode: 'copy'

    input:
        set val(name), file(fasta) from pmd_contigs_filter
    output:
        set val(name), file("*.filtered.fa") into pmd_filtered_fa, pmd_contig_filtered_dp, pmd_contigs_prokka
    script:
        filter_out = name+".filtered.fa"
        """
        filter_fasta_length.py -min ${params.pmd_minlen} -p ${task.cpus} $fasta -o $filter_out
        """
}

process align_reads_to_pmd_contigs {
    tag "$name"

    label 'intenso'

    errorStrategy 'ignore'

    publishDir "${params.results}/pmd_alignment/${name}", mode: 'copy'

    input:
        set val(name), file(contigs), val(name2), file(reads) from pmd_filtered_fa.combine(trimmed_reads_pmd_mapping)
    output:
        set val(name), file("*.sorted.bam") into pmd_alignment_to_dp
        file ("*.samtools.stats") into pmd_align_stats
    script:
        outfile = name+".sorted.bam"
        outstats = name + ".samtools.stats"
        if (params.pairedEnd) {
            """
            bowtie2-build --threads ${task.cpus} $contigs $name
            bowtie2 -x $name -1 ${reads[0]} -2 ${reads[1]} --threads ${task.cpus} | samtools view -S -b -F 4 - | samtools sort - > $outfile
            samtools stats $outfile > $outstats
            """
        } else {
            """
            bowtie2-build --threads ${task.cpus} $contigs $name
            bowtie2 -x $name -U $reads --threads ${task.cpus} | samtools view -S -b -F 4 - | samtools sort - > $outfile
            samtools stats $outfile > $outstats
            """
        }
}

process prokka {
    tag "$name"

    label 'intenso'

    publishDir "${params.results}/prokka", mode: 'copy'

    input:
        set val(name), file(contigs) from pmd_contigs_prokka
    output:
        file("prokka_result/*") into prokka_result
    script:
        """
        prokka --metagenome --cpus ${task.cpus} --outdir prokka_result $contigs
        """
}

process damageProfiler {
    tag "$name"

    label 'expresso'

    errorStrategy 'ignore'

    publishDir "${params.results}/damageProfiler/${name}", mode: 'copy'

    input:
        set val(name), file(bam), file(contig) from pmd_alignment_to_dp.join(pmd_contig_filtered_dp)
    output:
        file("*dmgprof.json") into dmgProf
    script:
        """
        damageprofiler -i $bam -r $contig -o tmp
        mv tmp/${name}.sorted/dmgprof.json ${name}.dmgprof.json
        """
}



process multiqc {
    label 'ristretto'

    publishDir "${params.results}/multiqc", mode: 'copy'

    input:
        file ('damageProfiler/*') from dmgProf.collect()
        file ('adapterRemoval/*') from adapter_removal_results.collect()
        file ('fastqc/*') from fastqc_results.collect()
        file ('samtools/*') from pmd_align_stats.collect()
        file ('quast/*') from quast_multiqc.collect()
        file ('prokka/*') from prokka_result.collect()
    output:
        file 'multiqc_report.html' into multiqc_report
    script:
        """
        multiqc -f -d fastqc adapterRemoval damageProfiler samtools quast prokka
        """
}




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
//         if (params.pairedEnd){
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