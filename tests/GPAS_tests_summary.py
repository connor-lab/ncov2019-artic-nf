#!/usr/bin/env python3
import sys, os
import pandas as pd
from argparse import ArgumentParser, SUPPRESS
import unittest
from Bio import SeqIO

## helper functions
def getName(s):
    if len(s.split(' '))<2:
        return 'nonesampleprocess'
    n=s.split(' ')[-1][1:-1]
    if n in ['nanopore','reference','illumina','nCoV-2019.reference.fasta']:
        return 'nonesampleprocess'
    return n

def getProcess(s):
    n=s.split(' ')[0]
    return n


def readSummary(f):
    headers=['percentage','reads','reads_taxon','taxon','taxid','name']
    df = pd.read_csv( f,  sep='\t', names=headers )
    df['name']=df['name'].str.strip()
    reads=df['reads_taxon'].sum()
    human_reads=df[df['name']=='Homo sapiens']['reads'].sum()
    COVID_reads=df[df['name']=='Coronavirinae']['reads'].sum()
    return reads, human_reads, COVID_reads
## tests

class runTests(unittest.TestCase):
    def __init__(self, opts):
        self.opts=opts
        self.run()

    def checkNextflow(self, path):
        # read trace, check all processes worked
        df=pd.read_csv('{}/trace.txt'.format(path),sep='\t')
        df=df[['name','status']]
        df['sample']=df['name'].map(getName)
        df['process']=df['name'].map(getProcess)
        g=df.groupby(['sample','process','status']).count()
        g=g.reset_index()
        g=g.pivot(index='sample',columns='process', values='status')
        
        # create list of samples in run
        self.samples=df['sample'].unique().tolist()
        self.samples.remove('nonesampleprocess')
        return g
    
    def checkKraken(self, path):
        # read summary, check human reads are removed
        sums=os.listdir('{0}/kraken/'.format(path))
        sums=[s for s in sums if s.endswith('_summary.txt')]
        self.krakensamples=[s.replace('_summary.txt','') for s in sums]
        self.assertCountEqual(self.krakensamples, self.samples)
        results=[]
        for s in sums:
            sampleName=s.replace('_summary.txt','')
            if sampleName not in self.samples: continue
            r,h,c=readSummary('{0}/kraken/{1}'.format(path, s))
            results.append({'sample':sampleName,'KRAKEN:Human reads':h,
                'KRAKEN:reads':r,
                'KRAKEN:Coronavirinae reads':c})
        results=pd.DataFrame(results)
        return results
    
    def checkConsensus(self, path):
        # check file is provided, correct length, correct bases, correct INDELs
        seqs=os.listdir('{0}/consensus_seqs/'.format(path))
        seqs=[s for s in seqs if s.endswith('.fasta')]
        self.conensus_seqs=[s.replace('.fasta','') for s in seqs]
#        self.assertCountEqual(self.conensus_seqs, self.samples)

        results=[]
        for s in seqs:
            sample=s.replace('.fasta','')
            seq=SeqIO.read(open('{0}/consensus_seqs/{1}'.format(path,s),'rt'), 'fasta')

            l=len(seq.seq)
            p=seq.seq.count('-')
            N=seq.seq.count('N')
            A=seq.seq.count('A')
            C=seq.seq.count('C')
            G=seq.seq.count('G')
            T=seq.seq.count('T')
            
            results.append({'sample':sample,
                    'CONSENSUS:length':l,
                    'CONSENSUS:masked (N)':N,
                    'CONSENSUS:padded (-)':p,
                    'CONSENSUS:bases (acgt)':sum([A,C,G,T])})
            
        results=pd.DataFrame(results)
        return results
    
    def checkPango(self, path):
        # Check correct lineage is called
        tech=os.listdir('{0}/analysis/pango/'.format(path))[0]
        seqs=os.listdir('{0}/analysis/pango/{1}/'.format(path,tech))
        self.pango_reports=[s.replace('_lineage_report.csv','') for s in seqs]
        dfs=[]
        for s in seqs:
            sample=s.replace('_lineage_report.csv','')
            df=pd.read_csv('{0}/analysis/pango/{1}/{2}'.format(path,tech,s))
            df['sample']=sample
            dfs.append(df)
        df=pd.concat(dfs)
        results=df[['sample','lineage','scorpio_call']]
        return results
    
    def checkNextclade(self, path):
        # Check mutations and lineage
        tech=os.listdir('{0}/analysis/nextclade/'.format(path))[0]
        seqs=os.listdir('{0}/analysis/nextclade/{1}/'.format(path,tech))
        seqs=[s for s in seqs if s.endswith('.tsv')]
        self.pango_reports=[s.replace('.tsv','') for s in seqs]
        dfs=[]
        for s in seqs:
            sample=s.replace('.tsv','')
            df=pd.read_csv('{0}/analysis/nextclade/{1}/{2}'.format(path,tech,s),sep='\t')
            df=df.add_prefix('NEXTCLADE:')
            df['sample']=sample
            dfs.append(df)
        results=pd.concat(dfs)
        return results
    
    def checkAln2type(self, path):
        # check VOCs correct
        tech=os.listdir('{0}/analysis/aln2type/'.format(path))[0]
        seqs=os.listdir('{0}/analysis/aln2type/{1}/'.format(path,tech))
        self.pango_reports=[s.replace('.csv','') for s in seqs]
        dfs=[]
        for s in seqs:
            sample=s.replace('.csv','')
            df=pd.read_csv('{0}/analysis/aln2type/{1}/{2}'.format(path,tech,s))
            df=df.add_prefix('ALN2TYPE:')
            df['sample']=sample
            dfs.append(df)
        results=pd.concat(dfs)
        return results

    def readExpected(self, f):
        df=pd.read_csv(f,sep='\t')
        return df

    
    def run(self):
        nf=self.checkNextflow(self.opts.workpath)
        kr=self.checkKraken(self.opts.outpath)
        df=nf.merge(kr,on='sample',how='left')
        cr=self.checkConsensus(self.opts.outpath)
        df=df.merge(cr,on='sample',how='left')
        pa=self.checkPango(self.opts.outpath)
        df=df.merge(pa,on='sample',how='left')
        nc=self.checkNextclade(self.opts.outpath)
        df=df.merge(nc,on='sample',how='left')
        al=self.checkAln2type(self.opts.outpath)
        df=df.merge(al,on='sample',how='left')
        
        if self.opts.expectedtsv != None:
            ex=self.readExpected(self.opts.expectedtsv)
            dfc=ex.compare(df)
            dfc.to_csv(self.opts.expectcomparisonedtsv, sep='\t', index=False)
            if len(dfc) == 0:
                print('Workflow finished with no discrepencies')
            else:
                print('Workflow finished with some discrepencies, see {}'.format(self.opts.expectcomparisonedtsv))
        df.to_csv(self.opts.outputtsv, sep='\t', index=False)

if __name__ == "__main__":
        # args
        parser = ArgumentParser(description='check for output results')
        parser.add_argument('-w', '--workpath', required=True,
                help='nextflow work path')
        parser.add_argument('-i', '--outpath', required=True,
                help='nextflow output path')
        parser.add_argument('-t', '--outputtsv', required=True,
                help='output tsv')
        parser.add_argument('-e', '--expectedtsv', required=False,
                help='expected results tsv')
        parser.add_argument('-c', '--expectcomparisonedtsv', required=False,
                help='expected results tsv')
        opts, unknown_args = parser.parse_known_args()
        runTests(opts)
