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
        set val(name), file('*.trimmed.fastq') into trimmed_reads_assembly, trimmed_reads_mapping
        file("*.settings") into adapter_removal_results

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

// process split_fasta {
//     tag "$name"

//     label 'intenso'

//     publishDir "${params.results}/fasta_filter/${name}", mode: 'copy'

//     input:
//         set val(name), file(fasta) from filtered_fa
//     output:
//         file("*.split.fa") into split_contigs
//     script:
//         filter_out = name+".filtered.fa"
//         """
//         split_fasta.py -p ${task.cpus} ${name}.filtered.fa
//         """
// }

// split_contigs
//     .flatten()
//     .map { it -> tuple(it.baseName, it) }
//     .into {filtered_contigs_bt; filtered_contigs_dp}

// trimmed_reads_mapping
//     .map {it -> it}
//     .set {trimmed_reads_mapping_ch}

    
// process bowtie_index_contigs{
//     tag "$name"

//     label 'intenso'

//     input:
//         set val(name), file(contig) from filtered_contigs_bt
//     output:
//         set val(name), file("*.bt2") into bt_index
//     script:
//         """
//         bowtie2-build --threads ${task.cpus} $contig $name
//         """
// }




process align_reads_to_contigs{
    tag "$name"

    label 'intenso'

    errorStrategy 'ignore'

    publishDir "${params.results}/alignment/${name}", mode: 'copy'

    input:
        // set val(read_name), file(reads) from trimmed_reads_mapping
        set val(name), file(contigs), val(name2), file(reads) from filtered_fa.combine(trimmed_reads_mapping)
    output:
        set val(name), file("*.sorted.bam") into alignment_to_dp, alignment_to_pmd
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

// process damageProfiler {
//     tag "$name"

//     label 'expresso'

//     errorStrategy 'ignore'

//     publishDir "${params.results}/damageProfiler/${name}", mode: 'copy'

//     input:
//         set val(name), file(bam), file(contig) from alignment_to_dp.join(filtered_contigs_dp)
//     output:
//         file("*dmgprof.json") into dmgProf
//     script:
//         """
//         damageprofiler -i $bam -r $contig -o tmp
//         mv tmp/${name}.sorted/dmgprof.json ${name}.dmgprof.json
//         """
// }



process PMDtools {
        tag "$name"

        label 'intenso'

        publishDir "${params.results}/pmdtools/", mode: 'copy', pattern: '*.fastq'

        input:
            set val(name), file(bam) from alignment_to_pmd
        output:
            set val(name), file("*.fastq") into pmd_assemble, pmd_map
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

    publishDir "${params.results}/pmd_assembly/${name}", mode: 'copy'

    input:
        set val(name), file(reads) from pmd_assemble

    output:
        set val(name), file("megahit_out/*.contigs.fa") into pmd_contigs_filter
        set val(name), file("megahit_out/*.log") into pmd_megahit_log
    script:
        mem = task.memory.toBytes()
        """
        cat *_pmd_R1.fastq > forward.fq
        cat *_pmd_R2.fastq > reverse.fq
        megahit -1 forward.fq -2 reverse.fq -t ${task.cpus} -m $mem --out-prefix $name
        """       
}

process multiqc {
    label 'ristretto'

    publishDir "${params.results}/multiqc", mode: 'copy'

    input:
        // file ('DamageProfiler/*') from dmgProf.collect()
        file ('AdapterRemoval/*') from adapter_removal_results.collect()
    output:
        file 'multiqc_report.html' into multiqc_report
    script:
        """
        multiqc -f -d AdapterRemoval
        """   
        // """
        // multiqc -f -d AdapterRemoval DamageProfiler
        // """   
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