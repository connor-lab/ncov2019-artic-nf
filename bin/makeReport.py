#!/usr/bin/env python3
import pandas as pd
import sys
import json

sample_name=sys.argv[1]

pango=pd.read_csv('pango.csv')
nextclade=pd.read_csv('nextclade.tsv', sep='\t')
aln2type=pd.read_csv('aln2type.csv')

pango['sample name']=sample_name
nextclade['sample name']=sample_name
aln2type['sample name']=sample_name

df=pango.merge(nextclade, on='sample name', how='left', suffixes=("_pango","_nextclade"))
df=df.merge(aln2type, on='sample name', how='left', suffixes=(None,"_aln2type"))

df.to_csv('{0}_report.tsv'.format(sys.argv[1]), sep='\t', index=False)

### convert to json
pango['program']='pango'
pango.set_index('program',inplace=True)
p=pango.to_dict(orient='index')

nextclade['program']='nextclade'
nextclade.set_index('program',inplace=True)
n=nextclade.to_dict(orient='index')

aln2type['program']='aln2type'
aln2type.set_index(['program','phe-label'],inplace=True)
a={level: aln2type.xs(level).to_dict('index') for level in aln2type.index.levels[0]}

d={sample_name:{}}
d[sample_name].update(p)
d[sample_name].update(n)
d[sample_name].update(a)

with open('{0}_report.json'.format(sample_name), 'w' ) as f:
    json.dump(d, f, indent=4, sort_keys=True)

