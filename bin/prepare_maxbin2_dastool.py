#!/usr/bin/env python

import argparse
from glob import glob
import os.path
import re

import pyfastx


if __name__ == '__main__':
    Parser = argparse.ArgumentParser(description='Prepare list of bin '
                                     'assignments for contigs by MaxBin2 as '
                                     'input into DASTool.')
    Parser.add_argument('-d', '--dir', required=True,
                        help='Output directory of MaxBin2')
    Parser.add_argument('-o', '--output', required=True,
                        help='list of bin assignments for each contig')
    Args = vars(Parser.parse_args())


    bins = glob(f"{Args['dir']}/*.fasta")
    with open(Args['output'], "wt") as outfile:
        for fafn in bins:
            binname = int(re.search(r'.+\.([0-9]+).fasta',
                                    os.path.basename(fafn)).group(1))
            for name, _ in pyfastx.Fasta(fafn, build_index=False):
                outfile.write(f"{name}\tbin_{binname:03d}\n")
