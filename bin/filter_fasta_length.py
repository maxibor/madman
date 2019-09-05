#!/usr/bin/env python

import argparse
import multiprocessing as mp
from functools import partial


def _get_args():
    '''This function parses and return arguments passed in'''
    parser = argparse.ArgumentParser(
        prog='fasta_filter_length',
        description='filterFastaByLength')
    parser.add_argument('fasta', help="path to input fasta file")
    parser.add_argument(
        '-min',
        default=1,
        help="Minimum sequence size. Default = 1")
    parser.add_argument(
        '-max',
        default=99999999999999999999999999999999,
        help="Maximim sequence size. Default = 99999999999999999999999999999999")
    parser.add_argument(
        '-p',
        default=1
        help='Number of process for multiprocessing. Default=1'
    )
    parser.add_argument(
        '-o',
        dest="output",
        default=None,
        help="Output file basename. Default = {basename}.filtered.fa")

    args = parser.parse_args()

    infile = args.fasta
    seqmin = int(args.min)
    seqmax = int(args.max)
    process = int(args.p)
    outfile = args.output

    return(infile, seqmin, seqmax, process, outfile)


def get_basename(file_name):
    if ("/") in file_name:
        basename = file_name.split("/")[-1].split(".")[0]
    else:
        basename = file_name.split(".")[0]
    return(basename)

def process_sequence(thekey, thedict, seqmin, seqmax, outdict):
    if len("".join(thedict[thekey])) >= seqmin and len("".join(thedict[thekey])) <= seqmax:
            outdict[thekey] = thedict[thekey]


if __name__ == "__main__":
    INFILE, SEQMIN, SEQMAX, PROCESS, OUTFILE = _get_args()

    basename = get_basename(INFILE)

    if not OUTFILE:
        OUTFILE = basename + ".filtered.fa"

    fastadict = {}
    fasta_filtered_dict = {}
    with open(INFILE, "r") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                seqname = line
                fastadict[seqname] = []
            else:
                fastadict[seqname].append(line)

    process_partial = partial(process_sequence, 
                              thedict = fasta_dict, 
                              seqmin=SEQMIN, 
                              seqmax=SEQMAX, 
                              outdict=fasta_filtered_dict)

    if (PROCESS > 1):
        with multiprocessing.Pool(PROCESS) as p:
            p.map(list(fasta_dict.keys(process_partial)))
    else :
        for seq in fastadict.keys():
            if len("".join(fastadict[seq])) >= SEQMIN and len("".join(fastadict[seq])) <= SEQMAX:
                fasta_filtered_dict[seq] = fastadict[seq]

    with open(OUTFILE, "w") as fw:
        for seq in fasta_filtered_dict.keys():
            fw.write(seq + "\n")
            fw.write("\n".join(fasta_filtered_dict[seq]) + "\n")