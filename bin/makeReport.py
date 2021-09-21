#!/usr/bin/env python3
import pandas as pd
import sys

pangolinTyping=pd.read_csv('pangolinTyping.csv')
nextclade=pd.read_csv('nextclade.tsv', sep='\t')

pangolinTyping['sample name']=sys.argv[1]
nextclade['sample name']=sys.argv[1]

df=pangolinTyping.merge(nextclade, on='sample name', how='left', suffixes=("_pango","_nextclade"))


#removing head and trail characters
#df['taxon'] = [x[10:-26] for x in df['taxon']]

# removing specific characters
df['taxon'] = df['taxon'].replace(regex=True, to_replace="Consensus_", value="")
df['taxon'] = df['taxon'].replace(regex=True, to_replace="_subsample.primertrimmed.consensus_threshold_0.75_quality_20", value="")
#df.drop(df.index[df['taxon'] == 'Consensus_'], inplace = True)

df.to_csv('{0}_report.tsv'.format(sys.argv[1]), sep='\t', index=False)