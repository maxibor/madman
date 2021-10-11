include { align_reads_to_contigs } from "$baseDir/modules/tools/bowtie2/main.nf" params(params)
include { damageprofiler as damageprofiler_pre; damageprofiler as damageprofiler_post } from "$baseDir/modules/tools/damageprofiler/main.nf" params(params)
include { filter_contigs_length } from "$baseDir/modules/tools/filter_contigs_length/main.nf" params(params)
include { filter_contigs_damage } from "$baseDir/modules/tools/filter_contigs_damage/main.nf" params(params)
include { megahit } from "$baseDir/modules/tools/megahit/main.nf" params(params)
include  {metaspades; biospades } from "$baseDir/modules/tools/spades/main.nf" params(params)
include { prokka } from "$baseDir/modules/tools/prokka/main.nf" params(params)
include { pydamage } from "$baseDir/modules/tools/pydamage/main.nf" params(params)
include { quast as quast_pre ; quast as quast_post } from "$baseDir/modules/tools/quast/main.nf" params(params)
include { freebayes } from "$baseDir/modules/tools/freebayes/main.nf" params(params)
include { bcftools } from "$baseDir/modules/tools/bcftools/main.nf" params(params)

workflow POST_ASSEMBLY_MODERN {
    versions = Channel.empty()
    take:
        contigs
        trimmed_reads
        assembler_name

    main:
        trimmed_reads
            .map {it -> ["${it[0]}_${assembler_name}", [file(it[1][0]),file(it[1][1])]]}
            .set { trimmed_reads }
        contigs
            .map {it -> ["${it[0]}_${assembler_name}", file(it[1])]}
            .set { contigs }

        quast_pre(contigs, "pre")
        filter_contigs_length(contigs)
        align_reads_to_contigs(filter_contigs_length.out.join(trimmed_reads))
        freebayes(filter_contigs_length.out.join(align_reads_to_contigs.out.bam))
        bcftools(filter_contigs_length.out.join(freebayes.out.vcf))
        prokka(bcftools.out.contigs)
        versions.mix(quast_pre.out.version, align_reads_to_contigs.out.version, bcftools.out.version, prokka.out.version)

    emit:
        quast_pre = quast_pre.out
        prokka = prokka.out
        versions


}

workflow POST_ASSEMBLY_ANCIENT {
    versions = Channel.empty()
    take:
        contigs
        trimmed_reads
        assembler_name

    main:
        trimmed_reads
            .map {it -> ["${it[0]}_${assembler_name}", [file(it[1][0]),file(it[1][1])]]}
            .set { trimmed_reads }
        contigs
            .map {it -> ["${it[0]}_${assembler_name}", file(it[1])]}
            .set { contigs }

        quast_pre(contigs, "pre")
        filter_contigs_length(contigs)
        align_reads_to_contigs(filter_contigs_length.out.join(trimmed_reads))
        freebayes(filter_contigs_length.out.join(align_reads_to_contigs.out.bam))
        bcftools(filter_contigs_length.out.join(freebayes.out.vcf))
        damageprofiler_pre(bcftools.out.contigs.join(align_reads_to_contigs.out.bam), "pre")
        pydamage(align_reads_to_contigs.out.bam)
        filter_contigs_damage(pydamage.out.filtered_csv.join(bcftools.out.contigs))
        quast_post(filter_contigs_damage.out.fasta, "post")
        damageprofiler_post(filter_contigs_damage.out.fasta.join(align_reads_to_contigs.out.bam),"post")
        prokka(bcftools.out.contigs)
        versions.mix(quast_pre.out.version, align_reads_to_contigs.out.version, freebayes.out.version, bcftools.out.version, damageprofiler_pre.out.version, prokka.out.version)

    emit:
        quast_pre = quast_pre.out.result
        quast_post = quast_post.out.result
        damageprofiler_pre = damageprofiler_pre.out.json
        damageprofiler_post = damageprofiler_post.out.json
        prokka = prokka.out.result
        versions

}
