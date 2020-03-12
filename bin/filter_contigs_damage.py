#!/usr/bin/env python

import argparse
from textwrap import wrap
import pandas as pd


def _get_args():
    '''This function parses and return arguments passed in'''
    parser = argparse.ArgumentParser(
        prog='fasta_filter_length',
        description='filterFastaByLength')
    parser.add_argument('fasta', help="path to input fasta contig file")
    parser.add_argument('pydamage', help="path to input pydamage csv file")
    parser.add_argument(
        '-d',
        dest='damage',
        default=0.2,
        type=float,
        help="Mimimum amount of CtoT damage on the 5' end of the read. Default=0.2"
    )
    parser.add_argument(
        '-a',
        dest='alpha',
        default=0.05,
        type=float,
        help='alpha threshold. Default=0.05'
    )
    parser.add_argument(
        '-o',
        dest="output",
        default=None,
        type=str,
        help="Output file basename. Default = {basename}.filtered.fa")

    args = parser.parse_args()

    contigs = args.fasta
    pydamage = args.pydamage
    alpha = args.alpha
    mindamage = args.damage
    outfile = args.output

    return(contigs, pydamage, alpha, mindamage, outfile)


def get_basename(file_name):
    if ("/") in file_name:
        basename = file_name.split("/")[-1].split(".")[0]
    else:
        basename = file_name.split(".")[0]
    return(basename)

def parse_fasta(fasta_file):
    """Parse a fasta file

    Args:
        fasta_file (str): Path to fasta file
    Returns:
        (dict): Seqname as key, sequence as value
    """   

    fastadict = {}
    with open(fasta_file, "r") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                seqname = line
                fastadict[seqname] = []
            else:
                fastadict[seqname].append(line)
    return(fastadict)

def write_fasta(fasta_dict, outfile):
    """Write fasta to file

    Args:
        fasta_dict (dict): Seqname as key, sequence as value
        outfile (str): Path to output file
    """    
    with open(outfile, "w") as fw:
        for seq in list(fasta_dict.keys()):
            fw.write(seq + "\n")
            fw.write("\n".join(wrap(("".join(fasta_dict[seq])), width = 80)))
            fw.write("\n")
        

def get_ancient_contigs(pydamage_report, alpha=0.05, mindamage=0.2):
    """Get name of contigs passing Q-value for ancient damage

    Args:
        pydamage_report (str): path to pydamage report csv file
        alpha (float): alpha threshold
    Returns:
        (list): list of contigs name passing threshold
    """        

    d = pd.read_csv(pydamage_report, index_col='reference')
    return(list(d.query(f"qvalue <= {alpha} and geom_pmax >= {mindamage}").index))

def filter_contigs(all_contigs, ancient_contigs):
    """Filter contigs if in ancient contigs

    Args:
        all_contigs(dict): fasta dict of contigs, seqname as key, sequence as value
        ancient contigs(list): list of ancient contigs names
    Returns:
        (dict): ancient contigs, seqname as key, sequence as value
    """
    a_contigs = {}
    for c in all_contigs:
        cname = c.split()[0][1:]
        if cname in ancient_contigs:
            a_contigs[c] = all_contigs[c]
    return(a_contigs)

if __name__ == "__main__":
    CONTIGS, PYDAMAGE, ALPHA, MINDAMAGE OUTFILE = _get_args()

    if not OUTFILE:
        OUTFILE = basename + ".filtered.fa"

    basename = get_basename(PYDAMAGE)
    all_contigs = parse_fasta(CONTIGS)
    ancient_contigs_names = get_ancient_contigs(PYDAMAGE, ALPHA, MINDAMAGE)
    ancient_contigs = filter_contigs(all_contigs, ancient_contigs_names)
    write_fasta(ancient_contigs, OUTFILE)
    


