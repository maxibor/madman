include { align_reads_to_contigs } from "$baseDir/modules/local/bowtie2.nf" params(params)
include { damageprofiler as damageprofiler_pre; damageprofiler as damageprofiler_post } from "$baseDir/modules/local/damageprofiler.nf" params(params)
include { filter_contigs_length } from "$baseDir/modules/local/filter_contigs_length.nf" params(params)
include { filter_contigs_damage } from "$baseDir/modules/local/filter_contigs_damage.nf" params(params)
include { megahit } from "$baseDir/modules/local/megahit.nf" params(params)
include  { metaspades } from "$baseDir/modules/local/metaspades.nf" params(params)
include  { biospades } from "$baseDir/modules/local/biospades.nf" params(params)
include { prokka } from "$baseDir/modules/local/prokka.nf" params(params)
include { pydamage } from "$baseDir/modules/local/pydamage.nf" params(params)
include { quast as quast_pre ; quast as quast_post } from "$baseDir/modules/local/quast.nf" params(params)
include { freebayes } from "$baseDir/modules/local/freebayes.nf" params(params)
include { bcftools } from "$baseDir/modules/local/bcftools.nf" params(params)

workflow POST_ASSEMBLY_MODERN {
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
        freebayes(filter_contigs_length.out.join(align_reads_to_contigs.out))
        bcftools(filter_contigs_length.out.join(freebayes.out.vcf))
        prokka(bcftools.out.contigs)

    emit:
        quast_pre = quast_pre.out
        prokka = prokka.out

}

workflow POST_ASSEMBLY_ANCIENT {
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
        freebayes(filter_contigs_length.out.join(align_reads_to_contigs.out))
        bcftools(filter_contigs_length.out.join(freebayes.out.vcf))
        damageprofiler_pre(bcftools.out.contigs.join(align_reads_to_contigs.out), "pre")
        pydamage(align_reads_to_contigs.out)
        filter_contigs_damage(pydamage.out.csv.join(bcftools.out.contigs))
        quast_post(filter_contigs_damage.out.fasta, "post")
        damageprofiler_post(filter_contigs_damage.out.fasta.join(align_reads_to_contigs.out),"post")
        prokka(bcftools.out.contigs)

    emit:
        quast_pre = quast_pre.out
        quast_post = quast_post.out
        damageprofiler_pre = damageprofiler_pre.out
        damageprofiler_post = damageprofiler_post.out
        prokka = prokka.out

}
