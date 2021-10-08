#!/usr/bin/env python

import argparse
import re

import pysam


if __name__ == '__main__':
    Parser = argparse.ArgumentParser(description='Fill in missing contig info '
                                     'that was dropped when subsetting the BAM '
                                     'file.')
    Parser.add_argument('--bam', required=True,
                        help='BAM file of short-read sequencing data aligned against assembled contigs')
    Parser.add_argument('--fasta', required=True,
                        help='FastA file with assembled contigs')
    Parser.add_argument('--contiglist', required=True,
                        help='list of contigs of bin')
    Parser.add_argument('--output', required=True,
                        help='output BAM file with corrected header')
    Args = vars(Parser.parse_args())

    # Load BAM file
    bamfile = pysam.AlignmentFile(Args['bam'])
    header = bamfile.header.to_dict()

    # Generate subset of SQ tags from contig list
    extract_sn_ln = re.compile(r'(.+):1-([0-9]+)')
    with open(Args['contiglist'], "rt") as contigfile:
        sq_info = [extract_sn_ln.search(line.rstrip()).groups()
                   for line in contigfile]
    header_sq = [{'SN': contig[0], 'LN': int(contig[1])}
                 for contig in sq_info]

    # Add contigs of mates aligning to different contig than in bin to header
    contig_ids = [contig[0] for contig in sq_info]
    contig_lengths = {name: len(seq)
                      for name, seq in pyfastx.Fasta(Args['fasta'], build_index=False)}
    for read in bamfile:
        if read.next_reference_name not in contig_ids and read.next_reference_name is not None:
            contig_ids.append(read.next_reference_name)
            header_sq.append({'SN': read.next_reference_name, 'LN': contig_lengths[read.next_reference_name]})
    header['SQ'] = header_sq
    bamfile.close()

    # Write header and reads back to file and fix read referenceID
    with pysam.AlignmentFile(Args['bam']) as bamfile:
        with pysam.AlignmentFile(Args['output'], "wbu", header=header) as outfile:
            for read in bamfile:
                read.reference_id = contig_ids.index(read.reference_name)
                if read.is_paired and read.next_reference_name is not None:
                    read.next_reference_id = contig_ids.index(read.next_reference_name)
                outfile.write(read)
    pysam.index(Args['output'])
