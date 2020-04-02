#!/usr/bin/env python3

from Bio import SeqIO
import csv
import subprocess

"""
This script can incorporate as many QC checks as required
as long as it outputs a csv file containing a final column
headed with 'qc_pass' and rows for each sample indcating
'TRUE' if the overall QC check has passed or 'FALSE' if not.
"""

def collect_covered_pos(bamfile, min_depth):
    p = subprocess.Popen(['samtools', 'depth', bamfile],
                       stdout=subprocess.PIPE)
    out, err = p.communicate()
    counter = 0
    for ln in out.decode('utf-8').split("\n"):
       if ln:
          contig, pos, depth = ln.split("\t")
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

    qc_lines = []

    depth_covered_bases = collect_covered_pos(args.bam, depth)

    pct_covered_bases = depth_covered_bases / ref_length * 100

    largest_n_gap = collect_largest_n_gap(args.fasta)

    if largest_n_gap >= 10000 or pct_covered_bases > 50.0:
        qc_pass = "TRUE"
       
    else:
        qc_pass = "FALSE"


    qc_lines.append({ 'sample_name' : args.sample, 
                'pct_covered_bases' : pct_covered_bases, 
                 'longest_no_N_run' : largest_n_gap,
                             'fasta': args.fasta, 
                              'bam' : args.bam,
                          'qc_pass' : qc_pass})


    with open(args.outfile, 'w') as csvfile:
        header = qc_lines[0].keys()
        writer = csv.DictWriter(csvfile, fieldnames=header)
        writer.writeheader()
        for line in qc_lines:
            writer.writerow(line)


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
