#!/usr/bin/env python3
import pandas as pd
import sys

pangolinTyping=pd.read_csv('pangolinTyping.csv')
nextclade=pd.read_csv('nextclade.tsv', sep='\t')

pangolinTyping['sample name']=sys.argv[1]
nextclade['sample name']=sys.argv[1]

df=pangolinTyping.merge(nextclade, on='sample name', how='left', suffixes=("_pango","_nextclade"))

df.to_csv('{0}_report.tsv'.format(sys.argv[1]), sep='\t', index=False)