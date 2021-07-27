#!/usr/bin/env python3
import os, sys, time,csv
import glob
import pandas as pd
import collections
import pathlib

data = sys.argv[1]
fastq_files = [f for f in os.listdir(data) if f.endswith('.fastq.gz')]

csv_files = glob.glob(os.path.join(data, "*.csv"))
for f in csv_files:
    #df = pd.read_csv(f, sep= ",", usecols=["files"], squeeze=False, encoding="UTF-8")
    df = pd.read_csv(f, sep= ",", usecols=["files"])
    ddf = df.files.to_list()
if collections.Counter(fastq_files) == collections.Counter(ddf):
    print("METADATA FILE AND FASTQ SAMPLE NAMES MATCHED...STARTING ANALYSIS")
else :
    print("PLEASE CHECK THE METADATA FILE AND FASTQ SAMPLE NAMES")