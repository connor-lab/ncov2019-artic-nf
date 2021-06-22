#!/usr/bin/env python3
from argparse import ArgumentParser, SUPPRESS
from Bio import SeqIO



def read_seqs(f,ref):
    '''return for sequence if it is not the reference genome, should only be two seqs'''
    for seq in SeqIO.parse(open(f,'rt'),'fasta'):
        if seq.id==ref: continue
        return seq

def rename_seq(seq,sample):
    seq.id='GPAS_SP3|{0}'.format(sample)
    return seq

def write_seqs(seq,out):
    SeqIO.write(seq,out,'fasta')

def run(opts):
    seq=read_seqs(opts.infile,opts.ref)
    seq=rename_seq(seq,opts.sample) 
    write_seqs(seq,opts.outfile)


if __name__ == "__main__":
    # args
    parser = ArgumentParser(description='create FN4 output fasta file from aligned MSA')
    parser.add_argument('-o', '--outfile', required=True,
                       help='output file')
    parser.add_argument('-i', '--infile', required=True,
                       help='input file')
    parser.add_argument('-r', '--ref', required=True,
                       help='reference id in MSA')
    parser.add_argument('-s', '--sample', required=True,
                       help='Sample name to add to header')

    opts, unknown_args = parser.parse_known_args()
    run(opts)
