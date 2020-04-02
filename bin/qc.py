#!/usr/bin/env python3

from Bio import SeqIO
import csv
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

"""
This script can incorporate as many QC checks as required
as long as it outputs a csv file containing a final column
headed with 'qc_pass' and rows for each sample indcating
'TRUE' if the overall QC check has passed or 'FALSE' if not.
"""

def make_depth_plot(depth_pos, samplename, window=200):
    df = pd.DataFrame({ 'position' : [pos[1] for pos in depth_pos], 'depth' : [dep[2] for dep in depth_pos]  })
    df['depth_moving_average'] = df.iloc[:,1].rolling(window=window).mean()
    plt.plot(df['depth_moving_average'],label='Depth')
    plt.legend(loc=2)
    plt.title(samplename)
    plt.savefig(samplename + '.depth.png')

def read_depth_file(bamfile):
    p = subprocess.Popen(['samtools', 'depth', '-a', '-d', '0', bamfile],
                       stdout=subprocess.PIPE)
    out, err = p.communicate()
    counter = 0

    pos_depth = []
    for ln in out.decode('utf-8').split("\n"):
       if ln:
          pos_depth.append(ln.split("\t"))
    
    return pos_depth


def collect_covered_pos(pos_depth, min_depth):
    counter = 0
    for contig, pos,depth in pos_depth:
        if int(depth) >= min_depth:
            counter = counter + 1
    
    return counter


def collect_largest_n_gap(fastafile):
    record = SeqIO.read(fastafile, "fasta")

    n_pos =  [i for i, letter in enumerate(record.seq.lower()) if letter == 'n']

    n_pos = [0] + n_pos

    n_pos = n_pos + [len(record.seq)] 

    n_gaps = [j-i for i, j in zip(n_pos[:-1], n_pos[1:])]

    return sorted(n_gaps)[-1]



def get_ref_length(ref):
    record = SeqIO.read(ref, "fasta")
    return len(record.seq)

def go(args):
    if args.illumina:
        depth = 10
    elif args.nanopore:
        depth = 20

    ref_length = get_ref_length(args.ref)
    depth_pos = read_depth_file(args.bam)

    depth_covered_bases = collect_covered_pos(depth_pos, depth)

    pct_covered_bases = depth_covered_bases / ref_length * 100

    largest_n_gap = collect_largest_n_gap(args.fasta)

    if largest_n_gap >= 10000 or pct_covered_bases > 50.0:
        qc_pass = "TRUE"
       
    else:
        qc_pass = "FALSE"


    qc_line = { 'sample_name' : args.sample, 
          'pct_covered_bases' : pct_covered_bases, 
           'longest_no_N_run' : largest_n_gap,
                       'fasta': args.fasta, 
                        'bam' : args.bam,
                    'qc_pass' : qc_pass}


    with open(args.outfile, 'w') as csvfile:
        header = qc_line.keys()
        writer = csv.DictWriter(csvfile, fieldnames=header)
        writer.writeheader()
        writer.writerow(qc_line)


    make_depth_plot(depth_pos, args.sample)

def main():
    import argparse

    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--nanopore', action='store_true')
    group.add_argument('--illumina', action='store_true')
    parser.add_argument('--outfile', required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--ref', required=True)
    parser.add_argument('--bam', required=True)
    parser.add_argument('--fasta', required=True)

    args = parser.parse_args()
    go(args)

if __name__ == "__main__":
    main()
