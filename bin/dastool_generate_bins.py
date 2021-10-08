#!/usr/bin/env python

import argparse
from bgzip import BGZipWriter
import os
from pathlib import Path

import pandas as pd
import pyfastx


if __name__ == '__main__':
    Parser = argparse.ArgumentParser(description='Export refined bins of '
                                     'DASTool that have an estimated completeness '
                                     'of >= 50% unifying the naming scheme.')
    Parser.add_argument('--fasta', required=True,
                        help='FastA file with assembled contigs')
    Parser.add_argument('--summary', required=True,
                        help='DASTool summary file')
    Parser.add_argument('--contiglist', required=True,
                        help='DASTool contig list')
    Parser.add_argument('--samplename', required=True,
                        help='sample name')
    Parser.add_argument('--outputdir', required=True,
                        help='output directory for bgzip-compressed FastA files of refined bins')
    Parser.add_argument('--map', required=True,
                        help='map of new bin labels and DASTool bin labels')
    Args = vars(Parser.parse_args())

    os.makedirs(Args['outputdir'], exist_ok=True)

    if os.stat(Args['summary']).st_size > 0:  # DASTool produced valid bins
        # Generate labels for bins
        binsummary = pd.read_csv(Args['summary'], sep="\t")
        binsummary = binsummary.loc[binsummary['SCG_completeness'] >= 0.5]
        binsummary['label'] = [f"{Args['samplename']}_{i + 1:03d}" for i in binsummary.index]

        # Write bins to file
        for cbin in binsummary.itertuples():
            contigs = pd.read_csv(Args['contiglist'], sep="\t",
                                  header=None, names=['contig', 'bin'])
            contigs = set(contigs.loc[contigs['bin'] == cbin.bin, 'contig'].tolist())
            with open(f"{Args['outputdir']}/{cbin.label}.fasta.gz", "wb") as outfile:
                with BGZipWriter(outfile) as outstream:
                    for name, seq in pyfastx.Fasta(Args['fasta'], build_index=False):
                        if name in contigs:
                            outstream.write(str.encode(f">{name}\n{seq}\n"))

        binsummary[['label', 'bin']].to_csv(Args['map'], sep="\t", index=False)
    else:
        Path(Args['map']).touch()
