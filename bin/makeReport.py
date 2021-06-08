#!/usr/bin/env python3
import pandas as pd
import sys

pango=pd.read_csv('pango.csv')
nextclade=pd.read_csv('nextclade.tsv', sep='\t')
aln2type=pd.read_csv('aln2type.csv')

pango['sample name']=sys.argv[1]
nextclade['sample name']=sys.argv[1]
aln2type['sample name']=sys.argv[1]

df=pango.merge(nextclade, on='sample name', how='left', suffixes=("_pango","_nextclade"))
df=df.merge(aln2type, left_on='sample name', right_on='sample_id', how='left', suffixes=(None,"_aln2type"))

df.to_csv('{0}_report.tsv'.format(sys.argv[1]), sep='\t', index=False)
