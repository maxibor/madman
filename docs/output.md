# nf-core/madman: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

* [nf-core/madman: Output](#nf-coremadman-output)
  * [**Output directory: `results/fastqc`**](#output-directory-resultsfastqc)
  * [**Output directory: `results/multiqc`**](#output-directory-resultsmultiqc)
  * [**Output directory: `assembly`**](#output-directory-assembly)
  * [**Output directory: `alignment`**](#output-directory-alignment)
  * [**Output directory: `fasta_filter`**](#output-directory-fastafilter)
  * [**Output directory: `pydamage`**](#output-directory-pydamage)
  * [**Output directory: `quast`**](#output-directory-quast)
  * [**Output directory: `prokka`**](#output-directory-prokka)

## **Output directory: `results/fastqc`**

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the `trim_galore` directory.

* `sample_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files
* `zips/sample_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images

## **Output directory: `results/multiqc`**

[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see [http://multiqc.info](http://multiqc.info)

## **Output directory: `assembly`**

Contains the output of de novo assembler:

* contigs (in FASTA format)
* logfile

## **Output directory: `alignment`**

Contains BAM files of the reads aligned back to the contigs

## **Output directory: `fasta_filter`**

Contains the filtered contigs:

* by size: `*.size_filtered.fa`
* by damage (pydamage) and size `*.ancient_filtered.fa`

## **Output directory: `pydamage`**

Contains the [pydamage](https://github.com/maxibor/pydamage) ancient DNA estimation report file

## **Output directory: `quast`**

Contains the [Quast](https://github.com/ablab/quast) assembly statistics report files.

This can help you assess the quality of the assemblies (contig lengths, N50s, chimeric contigs etc.).

## **Output directory: `prokka`**

Contains the [prokka](https://github.com/tseemann/prokka) annotation report files.

This allows assessment of assembly completeness, such as defined in with the [MIMAG criteria](https://doi.org/10.1038/nbt.3893)
