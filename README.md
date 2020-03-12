![megahit_nf_CI](https://github.com/maxibor/megahit-nf/workflows/megahit_nf_CI/badge.svg)

# megahit-nf

Nextflow assembly pipeline for ancient DNA

## Dependancies

- [Nextflow](https://www.nextflow.io/) : `conda install -c bioconda nextflow`

## Usage

```
nextflow run maxibor/megahit-nf -profile docker --reads "/path/to/paired_end_reads_*.{1,2}.fastq.gz"
```

## Input

### --reads

Use this to specify the location of your input FastQ files. For example:

`--reads 'path/to/data/sample_*_{1,2}.fastq'`

**Please note the following requirements:**

- The path must be enclosed in quotes
- The path must have at least one * wildcard character
- When using the pipeline with paired end data, the path must use {1,2} notation to specify read pairs.

## Arguments


### `--paired_end`

Specifies if reads are paired-end (true | false). Default: `true`

### `--phred`

Fastq quality encoding. Default: `33`

### `--minlen`

Minimum length of contig to keep it. Default: `500`

### `--minread`

Minimum number of reads aligned back to a contig to keep it. Default: `150`

### `--mindamage`

 Mimimum amount of CtoT damage on the 5' end of the aligned reads to keep the contig. Default=`0.2`

### `--results`

Path to directory to save results. Default: `results`

## Output

The output is generated in the `results` directory

### `multiqc_report.html`

Report summarizing the results

### `execution_report.html`

Report summarizing the execution of the pipeline

### `assembly`

Contains the output of megahit assembler:

- contigs
- logfile

### `alignment`

Contains the bam file of the reads aligned back to the contigs

### `fasta_filter`

Contains the contigs filter:

- by size: `*.size_filtered.fa`
- by damage (pydamage) and size `*.ancient_filtered.fa`

### `pydamage`

Contains the [pydamage](https://github.com/maxibor/pydamage) ancient DNA estimation report file

### `quast`

Contains the [Quast](https://github.com/ablab/quast) assembly statistics report files

### `prokka`

Contains the [prokka](https://github.com/tseemann/prokka) annotation report files


## Help

```
$ nextflow run maxibor/megahit-nf --help
N E X T F L O W  ~  version 19.10.0
Launching `maxibor/megahit-nf` [maniac_hypatia] - revision: cf6cdbd49c
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
  --phred                       Specifies the fastq quality encoding (33 | 64). Default: 33
  --paired_end                   Specifies if reads are paired-end (true | false). Default: null
  --minlen                      Minimum contig length to retain. Default:  1000
  --minread                     Minimum number of reads aligned to contig to consider contig. Default: 200

Options:
  --results                     The output directory where the results will be saved. Default: ./results
```
