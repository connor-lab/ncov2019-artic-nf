#!/usr/bin/env python3
import sys
import pandas as pd
from argparse import ArgumentParser, SUPPRESS


## helper functions
def getName(s):
    n=s.split(' ')[-1][1:-1]
    return n

def getProcess(s):
    n=s.split(' ')[0]
    return n

## tests
def checkNextflow(path):
    # read trace, check all processes worked
    df=pd.read_csv('{}/trace.txt'.format(path),sep='\t')
    df=df[['name','status']]
    df['sample']=df['name'].map(getName)
    df['process']=df['name'].map(getProcess)
    g=df.groupby(['sample','process','status']).count()
    g=g.reset_index()
    g=g.pivot(index='sample',columns='process', values='status')
    return g

def checkKraken(path):
    # read summary, check human reads are removed
    return results

def checkConsensus(path):
    # check file is provided, correct length, correct bases, correct INDELs
    return results

def checkPango(path):
    # Check correct lineage is called
    return results

def checkNextclade(path):
    # Check mutations and lineage
    return results

def checkAln2type(path):
    # check VOCs correct
    return results

def runTests(opts):
    df=checkNextflow(opts.path)
    print(df)

if __name__ == "__main__":
        # args
        parser = ArgumentParser(description='check for output results')
        parser.add_argument('-p', '--path', required=True,
                help='nextflow work path')
        opts, unknown_args = parser.parse_known_args()
        runTests(opts)
