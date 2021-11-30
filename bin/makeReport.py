#!/usr/bin/env python3
import pandas as pd
import sys
import json
from Bio import SeqIO

sample_name=sys.argv[1]

pango=pd.read_csv('pango.csv')
nextclade=pd.read_csv('nextclade.tsv', sep='\t')
aln2type=pd.read_csv('aln2type.csv')

pango['sampleName']=sample_name
nextclade['sampleName']=sample_name
aln2type['sampleName']=sample_name

df=pango.merge(nextclade, on='sampleName', how='left', suffixes=("_pango","_nextclade"))
df=df.merge(aln2type, on='sampleName', how='left', suffixes=(None,"_aln2type"))

# versions
wf=open('workflow_commit.txt').read()
df['workflowCommit']=str(wf).strip()
df['manifestVerison']=sys.argv[2]
nextclade_version=open('nextclade_files/version.txt').read()
df['nextcladeVersion']=str(nextclade_version).strip()
aln2type_variant_commit=open('variant_definitions/aln2type_variant_git_commit.txt').read()
aln2type_variant_version=open('variant_definitions/aln2type_variant_version.txt').read()
aln2type_source_commit=open('variant_definitions/aln2type_commit.txt').read()
df['aln2typeVariantCommit']=str(aln2type_variant_commit).strip()
df['aln2typeVariantVersion']=str(aln2type_variant_version).strip()
df['aln2typeSourceVommit']=str(aln2type_source_commit).strip()

df.to_csv('{0}_report.tsv'.format(sys.argv[1]), sep='\t', index=False)

### convert to json
pango['program']='pango'
pango.set_index('program',inplace=True)
p=pango.to_dict(orient='index')

nextclade['program']='nextclade'
nextclade['nextcladeVersion']=str(nextclade_version).strip()
nextclade.set_index('program',inplace=True)
n=nextclade.to_dict(orient='index')
with open('nextclade.json','rt', encoding= 'utf-8') as inf:
    nj=json.load(inf)
n['nextcladeOutputJson']=nj

aln2type['program']='aln2type'
aln2type['label']=aln2type['phe-label']
aln2type['aln2typeVariantCommit']=str(aln2type_variant_commit).strip()
aln2type['aln2typeSourceCommit']=str(aln2type_source_commit).strip()
aln2type.set_index(['program','phe-label'],inplace=True)
a={level: aln2type.xs(level).to_dict('index') for level in aln2type.index.levels[0]}

w={'WorkflowInformation':{}}
w['WorkflowInformation']['workflowCommit']=str(wf).strip()
w['WorkflowInformation']['manifestVerison']=sys.argv[2]
w['WorkflowInformation']['sampleIdentifier']=sample_name

# add fasta to json
record = SeqIO.read('ref.fasta', "fasta")
w['WorkflowInformation']['referenceIdentifier']=record.id
#f={'FastaRecord':{'SeqId':record.id,
#    'SeqDescription': record.description,
#    'Sequence':str(record.seq),
#    'sampleName':sample_name}}


d={sample_name:{}}
d[sample_name].update(p)
d[sample_name].update(n)
d[sample_name].update(a)
d[sample_name].update(w)
#d[sample_name].update(f)

with open('{0}_report.json'.format(sample_name), 'w', encoding='utf-8') as f:
    json.dump(d, f, indent=4, sort_keys=True, ensure_ascii=False)

