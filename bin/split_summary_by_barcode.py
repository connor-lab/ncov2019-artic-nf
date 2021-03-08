#!/usr/bin/env python3
import pandas as pd
import sys,os
from tqdm import tqdm
import multiprocessing as mp
import datetime as dt

dtypes = {
    'filename_fastq': 'object',
    'filename_fast5': 'object',
    'read_id': 'object',
    'run_id': 'category',
    'channel': 'int64',
    'mux': 'int64',
    'start_time': 'float64',
    'duration': 'float64',
    'num_events': 'int64',
    'passes_filtering': 'bool',
    'template_start': 'float64',
    'num_events_template': 'int64',
    'template_duration': 'float64',
    'sequence_length_template': 'int64',
    'mean_qscore_template': 'float64',
    'strand_score_template': 'float64',
    'median_template': 'float64',
    'mad_template': 'float64',
    'pore_type': 'category',
    'experiment_id': 'category',
    'sample_id': 'category',
    'end_reason': 'category',
    'alias': 'category',
    'type': 'category',
    'barcode_arrangement': 'category',
    'barcode_full_arrangement': 'category',
    'barcode_kit': 'category',
    'barcode_variant': 'category',
    'barcode_score': 'float64',
    'barcode_front_id': 'category',
    'barcode_front_score': 'float64',
    'barcode_front_refseq': 'category',
    'barcode_front_foundseq': 'category',
    'barcode_front_foundseq_length': 'int64',
    'barcode_front_begin_index': 'int64',
    'barcode_rear_id': 'category',
    'barcode_rear_score': 'float64',
    'barcode_rear_refseq': 'object',
    'barcode_rear_foundseq': 'object',
    'barcode_rear_foundseq_length': 'int64',
    'barcode_rear_end_index': 'int64',
    'bc_front': 'category',
    'bc_rear': 'category'}

def process_df(df, out_dir):
    gb = df.groupby('barcode_arrangement')
    for barcode in tqdm(gb.groups):
        #os.makedirs('{0}/{1}'.format(out_dir, barcode ), exist_ok=True)
        gb.get_group(barcode).to_csv('{0}/{1}.txt'.format(out_dir,barcode),
            sep='\t',
            mode="a",
            index=False)

def split_summary_by_barcode_chunk(summary_path, out_dir):
    '''Given a sequencing summary file path, write per barcode summaries to an output directory'''
    start_time = dt.datetime.now()
    os.makedirs(out_dir, exist_ok=True)
    pool = mp.Pool(4)
    df_iter = pd.read_csv(summary_path, sep='\t', dtype=dtypes, iterator=True, chunksize=100000)
    for df in df_iter:
        pool.apply_async(process_df(df, out_dir), [df])
    end_time = dt.datetime.now()

    elapsed_time = end_time - start_time
    print(f"Elapsed time {elapsed_time}")
seqsum=sys.argv[1]
outdir=sys.argv[2]
split_summary_by_barcode_chunk(seqsum, outdir)
