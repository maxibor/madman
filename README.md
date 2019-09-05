# megahit-nf

Simple megahit2 Nextflow pipeline

## Dependancies

- [conda](https://conda.io/en/latest/) 
- [Nextflow](https://www.nextflow.io/) : `conda install -c bioconda nextflow`

## Usage

```
nextflow run maxibor/megahit-nf --reads "/path/to/paired_end_reads_*.{1,2}.fastq.gz"
```

### Input

#### --reads

Use this to specify the location of your input FastQ files. For example:

`--reads 'path/to/data/sample_*_{1,2}.fastq'`

**Please note the following requirements:**

- The path must be enclosed in quotes
- The path must have at least one * wildcard character
- When using the pipeline with paired end data, the path must use {1,2} notation to specify read pairs.


### Output


## Help

```
$ nextflow run maxibor/megahit-nf --help
N E X T F L O W  ~  version 19.04.0
Launching `maxibor/megahit-nf` [compassionate_morse] - revision: f534a6a703
megahit-nf: simple Megahit assembler Nextflow pipeline
 Homepage: https://github.com/maxibor/megahit-nf
 Author: Maxime Borry <borry@shh.mpg.de>
=========================================
Usage:
The typical command for running the pipeline is as follows:
nextflow run maxibor/megahit-nf --reads "/path/to/paired_end_reads_*.{1,2}.fastq.gz"
Mandatory arguments:
  --reads                       Path to input data (must be surrounded with quotes)

Settings:
  --phred                       Specifies the fastq quality encoding (33 | 64). Defaults to 33
  --pairedEnd                   Specified if reads are paired-end (true | false). Default = true
  --minlen                      Minimum contig length to retain. Default =  1000

Options:
  --results                     The output directory where the results will be saved. Defaults to ./results
  --help  --h                   Shows this help page
```
