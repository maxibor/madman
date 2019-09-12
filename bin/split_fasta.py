#!/usr/bin/env python

import argparse
import multiprocessing as mp
from functools import partial


def _get_args():
    '''This function parses and return arguments passed in'''
    parser = argparse.ArgumentParser(
        prog='split_fasta',
        description='split fasta by sequence')
    parser.add_argument('fasta', help="path to input fasta file")
    parser.add_argument(
        '-p',
        default=1,
        help='Number of process for multiprocessing. Default=1'
    )

    args = parser.parse_args()

    infile = args.fasta
    process = int(args.p)

    return(infile, process)


def read_fasta(infile):
    """
    READ FASTA FILE INTO DICT
    Args:
        infile (str): Path to input fasta file
    
    Returns:
        seqdict (dict of lists): Dictionary of sequences, key as seqname,
            sequences in list. 
    """
    with open(infile, 'r') as f:
        seqdict ={}
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                seqname = line[1:].split()[0]
                seqdict[seqname] = []
                continue
            else:
                seqdict[seqname].append(line)
        return(seqdict)

def write_fasta(seqname, seqdict):
    """
    Writes fasta dict to file
    Args:
        seqname (str): name of sequence
        seqdict (dict of lists): Dictionary of sequences, key as seqname,
            sequences in list. 
    """

    with open(f"{seqname}.split.fa", "w") as fw:
        fw.write(f">{seqname}\n")
        for line in seqdict[seqname]:
            fw.write(f"{line}\n")


if __name__ == "__main__":
    INFILE, PROCESS = _get_args()

    fasta_dict = read_fasta(INFILE)

    write_partial = partial(write_fasta, seqdict = fasta_dict)

    if (PROCESS > 1):
        with mp.Pool(PROCESS) as p:
            p.map(write_partial, list(fasta_dict.keys()))
    else :
        for seqname in list(fasta_dict.keys()):
            write_fasta(seqname=seqname, seqdict=fasta_dict)