################################################################################
# Snakemake pipeline scaffold for MAG binning and quality checks
#
# Alex Huebner, 08/10/21
################################################################################


rule all:
    input: 
        expand("dastool/{sample}/{sample}_DASTool_{files}.txt", sample=SAMPLES, files=['labels', 'summary']),
        expand("gunc/{sample}/GUNC_checkM.merged.tsv", sample=SAMPLES),
        expand("cmseq/{bin}.{analysis}.txt", bin=BINS, analysis=['breadth_depth', 'polymut'])


#### Determine depth ###########################################################

rule index:
    output:
        "alignment/{sample}_megahit/{sample}_megahit.sorted.bam.bai"
    message: "Index BAM file of {wildcards.sample}"
    conda: "scaffold_binning.yaml"
    params:
        bam = "alignment/{sample}_megahit/{sample}_megahit.sorted.bam"  # output of modules/tools/bowtie2/main.nf
    shell:
        """
        samtools index {params.bam}
        """

rule depth_metabat2:
    input:
        "alignment/{sample}_megahit/{sample}_megahit.sorted.bam.bai"
    output:
        "metabat2/{sample}/{sample}.depth.txt"
    message: "Calculate the mean and variance of the coverage per contig: {wildcards.sample}"
    conda: "scaffold_binning.yaml"
    params:
        min_contig_length = 1000 # defined by variable params.minlen
        bam = "alignment/{sample}_megahit/{sample}_megahit.sorted.bam"  # output of modules/tools/bowtie2/main.nf
    shell:
        """
        jgi_summarize_bam_contig_depths \
            --outputDepth {output} \
            --minContigDepth 1 \
            --minContigLength {params.min_contig_length} \
            {input.bam}
        """

rule depth_maxbin2:
    input:
        "metabat2/{sample}/{sample}.depth.txt"
    output:
        "maxbin2/{sample}/{sample}.depth.txt"
    message: "Prepare depth for MaxBin2: {wildcards.sample}"
    conda: "scaffold_binning.yaml"
    shell:
        """
        bioawk -t '{{if (NR > 1){{print $1, $3}}}}' {input} > {output}
        """

################################################################################

#### Binning ###################################################################

rule metabat2:
    input:
        "metabat2/{sample}/{sample}.depth.txt"
    output:
        "metabat2/{sample}/{sample}_contigBin"
    message: "Bin contigs of {wildcards.sample} with metaBAT2"
    conda: "scaffold_binning.yaml"
    params:
        fasta = "assembly/{assembler}/{sample}/{assembler}_out/{sample}_{assembler}.contigs.fa" # output of modules/tools/megahit/main.nf or modules/tools/spades/main.nf
    threads: 16
    shell:
        """
        metabat2 -i {params.fasta} \
                 -o {output} \
                 -a {input} \
                 --minContig 1500 \
                 --minClsSize 100000 \
                 --numThreads {threads} \
                 --verbose && touch {output}
        """

rule maxbin2:
    input:
        "maxbin2/{sample}/{sample}.depth.txt"
    output:
        "maxbin2/{sample}/{sample}.summary"
    message: "Cluster contigs using MaxBin2 for sample {wildcards.sample}"
    conda: "scaffold_binning.yaml"
    params:
        fasta = "assembly/{assembler}/{sample}/{assembler}_out/{sample}_{assembler}.contigs.fa", # output of modules/tools/megahit/main.nf or modules/tools/spades/main.nf
        outprefix = "tmp/maxbin2/{sample}/{sample}",
        outdir = "maxbin2/{sample}"
    threads: 8
    shell:
        """
        mkdir $(dirname {params.outprefix})
        run_MaxBin.pl -contig {params.fasta} \
                      -abund {input} \
                      -thread {threads} \
                      -verbose \
                      -out {params.outprefix}
        cp -r {params.outprefix}*.fasta {params.outprefix}.{{marker,summary}} {params.outdir}/
        rm -r $(dirname {params.outprefix})
        """

################################################################################

#### Bin refinement ############################################################

rule prepare_metabat2:
    input:
        "metabat2/{sample}/{sample}_contigBin"
    output:
        "dastool/{sample}.metaBAT2.scaffolds2bin"
    message: "Generate table of contigs and their respective bins for metaBAT2: {wildcards.sample}"
    conda: "scaffold_binning.yaml"
    params:
        dir = "metabat2/{sample}"
    shell:
        """
        bin/prepare_metabat2_dastool.py -d {params.dir} -o {output}
        """

rule prepare_maxbin2:
    input:
        "maxbin2/{sample}/{sample}.summary"
    output:
        "dastool/{sample}.maxbin2.scaffolds2bin"
    message: "Generate table of contigs and their respective bins for maxBin2: {wildcards.sample}"
    conda: "scaffold_binning.yaml"
    params:
        dir = "maxbin2/{sample}"
    shell:
        """
        bin/prepare_maxbin2_dastool.py -d {params.dir} -o {output}
        """

rule dastool:
    input:
        metabat2 = "dastool/{sample}.metaBAT2.scaffolds2bin",
        maxbin2 = "dastool/{sample}.maxbin2.scaffolds2bin"
    output:
        "dastool/{sample}/{sample}_DASTool_summary.txt"
    message: "Redefine bins using DASTool for sample {wildcards.sample}"
    conda: "scaffold_binning.yaml"
    params:
        dir = "dastool/{sample}/{sample}",
        fasta = "assembly/{assembler}/{sample}/{assembler}_out/{sample}_{assembler}.contigs.fa", # output of modules/tools/megahit/main.nf or modules/tools/spades/main.nf
    threads: 8
    shell:
        """
        DAS_Tool -i {input.metabat2},{input.maxbin2} \
                 -c {params.fasta} \
                 -l metaBAT2,MaxBin2 \
                 -o {params.dir} \
                 -t {threads}
        touch {output}
        """

rule generate_bins:
    input:
        "dastool/{sample}/{sample}_DASTool_summary.txt"
    output:
        "dastool/{sample}/{sample}_DASTool_labels.txt"
    message: "Write all bins of sample {wildcards.sample} to file"
    conda: "scaffold_binning.yaml"
    params:
        fasta = "assembly/{assembler}/{sample}/{assembler}_out/{sample}_{assembler}.contigs.fa", # output of modules/tools/megahit/main.nf or modules/tools/spades/main.nf
        outdir = "bins",
        scaffolds2bin = "dastool/{sample}/{sample}_DASTool_scaffolds2bin.txt"
    shell:
        """
        bin/dastool_generate_bins.py \
            --fasta {params.fasta} \
            --summary {input} \
            --contiglist {params.scaffolds2bin} \
            --samplename {wildcards.sample} \
            --outputdir {params.outdir} \
            --map {output}
        """

################################################################################

#### CheckM ####################################################################

rule checkM_prepareDatabase:
    output:
        "checkM/.dmanifest"
    message: "Download the checkM databases and extract the tarballs"
    params:
        outdir = "checkM",
        url = "https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz"
    shell:
        """
        cd {params.outdir}
        wget {params.url}
        tar xvf checkm_data_2015_01_16.tar.gz
        """

rule checkM_setRoot:
    input:
        "checkM/.dmanifest"
    output:
        touch("checkM/setRoot.done")
    message: "Specify the database folder in checkM"
    conda: "scaffold_binning.yaml"
    params:
        dbdir = "checkM" 
    shell:
        """
        echo {params.dbdir} | checkm data setRoot {params.dbdir}
        """

rule unzip_compressed_fastas:
    output:
        touch("tmp/checkM/{sample}/{sample}.done")
    message: "De-compress the FastA files of bins belonging to sample {wildcards.sample}"
    params:
        tmpdir = "tmp/checkM/{sample}",
        dastool = "bins/{sample}"
    shell:
        """
        for bin in $(ls {params.dastool}_*.fasta.gz); do
            gunzip -c ${{bin}} > {params.tmpdir}/$(basename ${{bin}} fasta.gz)fa
        done
        """

rule checkM:
    input:
        db = "checkM/setRoot.done",
        fasta = "tmp/checkM/{sample}/{sample}.done" 
    output:
        "checkM/{sample}/storage/marker_gene_stats.tsv"
    message: "Run checkM using lineage-specific workflow on sample {wildcards.sample}"
    conda: "scaffold_binning.yaml"
    params: 
        bindir = "tmp/checkM/{sample}",
        outdir = "checkM/{sample}"
    log: "checkM/{sample}/{sample}.checkM.log"
    threads: 8
    shell:
        """
        rm -r {params.outputfd} # Remove output folder automatically created by Snakemake
        checkm lineage_wf \
                -x fa \
                -t {threads} \
                {params.bindir} {params.outdir} > {log} && \
        """

rule checkM_extendedreport_table:
    input:
        "checkM/{sample}/storage/marker_gene_stats.tsv"
    output:
        "checkM/{sample}/{sample}.checkM.txt"
    message: "Generate extended checkM report for sample {wildcards.sample}"
    conda: "scaffold_binning.yaml"
    params:
        outdir = "checkM/{sample}"
    shell:
        """
        checkm qa -o 2 --tab_table -f {output} \
                {params.outdir}/lineage.ms \
                {params.outdir}
        """

################################################################################

#### GUNC ######################################################################

rule install_gunc_database:
    output:
        "gunc/db/gunc_db_progenomes2.1.dmnd"
    message: "Install GUNC database"
    conda: "scaffold_binning.yaml"
    params:
        dir = "gunc/db"
    shell:
        """
        gunc download_db {params.dir}
        """

rule gunc_run:
    input:
        fastadir = "tmp/checkM/{sample}/{sample}.done",
        db = "gunc/db/gunc_db_progenomes2.1.dmnd"
    output:
        "gunc/{sample}/GUNC.progenomes_2.1.maxCSS_level.tsv"
    message: "Run GUNC on bins: {wildcards.sample}"
    conda: "scaffold_binning.yaml"
    params:
        fadir = "tmp/checkM/{sample}",
        outdir = "gunc/{sample}",
        tmpdir = "gunc/{sample}/tmp"
    threads: 8
    shell:
        """
        mkdir -p {params.tmpdir}
        gunc run \
            -r {input.db} \
            --input_dir {params.fadir} \
            -t {threads} \
            -o {params.outdir} \
            --temp_dir {params.tmpdir} \
            --detailed_output
        rmdir {params.tmpdir}
        """

rule gunc_merge:
    input:
        gunc = "gunc/{sample}/GUNC.progenomes_2.1.maxCSS_level.tsv",
        checkm = "checkM/{sample}/{sample}.checkM.txt"
    output:
        "gunc/{sample}/GUNC_checkM.merged.tsv"
    message: "Merge GUNC and CheckM output: {wildcards.sample}"
    conda: "scaffold_binning.yaml"
    params:
        dir = "gunc/{sample}"
    shell:
        """
        gunc merge_checkm \
            -g {input.gunc} \
            -c {input.checkm} \
            -o {params.dir}
        """

################################################################################

#### CMSeq #####################################################################

rule generate_contiglist_per_bin:
    input:
        "dastool/{sample}/{sample}_DASTool_labels.txt"
    output:
        temp("cmseq/{bin}/{bin}.contiglist.txt"),
    message: "Generate list of contig names belonging to bin {wildcards.bin} for subsetting the BAM file"
    conda: "scaffold_binning.yaml"
    params:
        fasta = "dastool/bins/{bin}.fasta.gz"
    shell:
        """
        bioawk -c fastx '{{print $name ":1-" length($seq)}}' {params.fasta} > {output}
        """

rule subset_bam:
    input:
        "cmseq/{bin}/{bin}.contiglist.txt",
    output:
        temp("cmseq/{bin}/{bin}.subset.bam"),
    message: "Subset BAM file for contigs of bin {wildcards.bin}"
    conda: "scaffold_binning.yaml"
    params:
        bam = "alignment/{sample}_megahit/{sample}_megahit.sorted.bam"  # output of modules/tools/bowtie2/main.nf
    shell:
        """
        cat {input} | xargs samtools view -hb {params.bam} > {output}
        """

rule clean_header:
    input:
        contigs = "cmseq/{bin}/{bin}.contiglist.txt",
        bam = "cmseq/{bin}/{bin}.subset.bam"
    output:
        bam = temp("cmseq/{bin}/{bin}.bam"),
        bai = temp("cmseq/{bin}/{bin}.bam.bai")
    message: "Remove contigs not belonging to bin to allow for faster iteration: {wildcards.bin}"
    conda: "scaffold_binning.yaml"
    params:
        fasta = lambda wildcards: f"assembly/{assembler}/{wildcards.bin.split('_')[0]}/{assembler}_out/{wildcards.bin.split('_')[0]}_{assembler}.contigs.fa" # output of modules/tools/megahit/main.nf or modules/tools/spades/main.nf
    shell:
        """
        bin/cmseq_fill_header.py \
            --bam {input.bam} \
            --fasta {params.fasta} \
            --contiglist {input.contigs} \
            --output {output.bam}
        """

rule breadth_depth:
    input: 
        bam = "cmseq/{bin}/{bin}.bam",
        bai = "cmseq/{bin}/{bin}.bam.bai"
    output:
        "cmseq/{bin}.breadth_depth.txt"
    message: "Investigate the breadth and depth of sample {wildcards.bin}"
    conda: "scaffold_binning.yaml"
    shell:
        """
        breadth_depth.py \
            -f \
            --minqual 30 \
            --mincov 1 \
            {input.bam} > {output}
        """

rule polymut:
    # Used by Pasolli et al. (2019); returns three-column table: number of
    # non-syn. mutations, number of syn. mutations, total number of positions
    input:
        bam = "cmseq/{bin}/{bin}.bam",
        bai = "cmseq/{bin}/{bin}.bam.bai",
    output:
        "cmseq/{bin}.polymut.txt"
    message: "Estimate the polymorphic rate of sample {wildcards.bin}"
    conda: "scaffold_binning.yaml"
    params:
        gff = "prokka/{sample}_{assember}/{sample}_{assember}.gff" # from modules/tools/prokka/main.nf
    shell:
        """
        polymut.py \
            -f \
            --gff_file {params.gff} \
            --minqual 30 \
            --mincov 3 \
            --dominant_frq_thrsh 0.8 \
            {input.bam} > {output}
        """

################################################################################
