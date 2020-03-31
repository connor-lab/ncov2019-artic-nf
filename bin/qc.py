import sys

from Bio import SeqIO

record = SeqIO.read(sys.argv[1], "fasta")


n_pos =  [i for i, letter in enumerate(record.seq.lower()) if letter == 'n']

n_gaps = [j-i for i, j in zip(n_pos[:-1], n_pos[1:])]

print(sorted(n_gaps)[-1])
