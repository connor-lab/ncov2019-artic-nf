#!/usr/bin/env python3

import argparse
import pandas as pd

def init_parser():
    '''
    Parser Arguments to pass to script from CL
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--base_dir',
        required=True,
        help='Nextflow out directory'
    )
    parser.add_argument(
        '--sample_info',
        required=True,
        help='Sampleinfo.csv file path'
    )
    parser.add_argument(
        '--prefix',
        required=True,
        help='Pass prefix string'
    )
    parser.add_argument(
        '--illumina',
        default=False,
        action='store_true',
        required=False,
        help='Specify if input is illumina data'
    )
    
    return parser

def parse_sample_csv(sample_csv, prefix):
    '''
    Take input sampleinfo.csv file and turn it into a usable format
    '''

    # Set up dataframe
    df_out = pd.DataFrame(columns=['Sample_Name', 'Project_ID', 'File_Forward', 'File_Reverse'])

    # Open file and populate df
    with open(sample_csv) as input_handle:
        for index, line in enumerate(input_handle):

            if index == 0:
                continue
            
            current_line_list = line.strip().split(',') # Order is [Name, run, barcode, project_id]


            # Set correct file name
            if len(current_line_list[2]) != 2: # Checking that barcode is 2 digits
                barcode = '0{}'.format(current_line_list[2])

                if len(barcode) != 2:
                    print('Error with barcodes in SampleInfo.csv')
                    quit()
            
            else:
                barcode = current_line_list[2]

            file_name = '{}_barcode{}.fastq'.format(prefix, barcode)
            

            # Set DataFrame for easy output csv

            df_out.at[index, 'Sample_Name'] = current_line_list[0] # Name from input sample info file
            df_out.at[index, 'Project_ID'] = current_line_list[3] # Project number from sample info file
            df_out.at[index, 'File_Forward'] = file_name

    return df_out


def main():
    '''
    Main script
    '''
    
    # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    base_dir = args.base_dir
    sample_csv = args.sample_info
    prefix = args.prefix

    with open('SampleList.csv', 'w') as handle:
        handle.write('[DATA]\n')

    # Start processing inputs
    if args.illumina: # Illumina needs to be paired I believe but I don't remember at the moment
        print('Illumina processing is not set up yet')
        quit()


    else: # Nanopore data

        df_out = parse_sample_csv(sample_csv, prefix)
        
        df_out.to_csv("SampleList.csv", mode='a', header=True, index=False)


if __name__ == "__main__":
    main()