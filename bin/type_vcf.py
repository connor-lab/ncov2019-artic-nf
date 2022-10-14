#!/usr/bin/env python3
import csv
import os
import sys
import re
import gzip
import io
import errno
import argparse
from subprocess import run, PIPE, DEVNULL
import yaml
import vcf

def parse_args(args=None):
    Description = 'Parse variant call files and call variation consequences.'

    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument('-i', '--id', dest='SAMPLE_ID', required=True, help="Sample ID")
    parser.add_argument('-y', '--yml', dest="YML_IN", required=False, help="YAML file with types and their variants")
    parser.add_argument('-ov', '--output_vcf', dest="VCF_OUT", required=True, help="Output VCF file with typed variants")
    parser.add_argument('-ot', '--output_typing_csv', dest="TYPING_CSV_OUT", required=True, help="Output CSV file with typed variants")
    parser.add_argument('-os', '--output_summary_csv', dest="SUMMARY_CSV_OUT", required=True, help="Output variant summary CSV file")
    parser.add_argument('-af', '--allele_freq_thresh', type=float, dest="ALLELE_FREQ_THRESH", default=0, help="Only output variants where allele frequency greater than this number (default: 0).")
    parser.add_argument('-dp', '--minimum_depth', type=int, dest="MIN_DEPTH", default=0, help="Only output variants where overall read depth is greater than this number (default: 0).")
    infile = parser.add_mutually_exclusive_group(required=True)
    infile.add_argument('-v', '--vcf', dest='VCF_IN', required=False, help="Input VCF file from either ARTIC pipeline (Nanopolish or Medaka).")
    infile.add_argument('-t', '--tab', dest='TSV_IN', required=False, help="Input iVar tsv file.")
    parser.add_argument('GFF_IN', help="Annotation GFF3 file.")
    parser.add_argument('REF_IN', help="Reference fasta file.")

    return parser.parse_args(args)

def ivar_variants_to_vcf_string(FileIn,RefIn):
    '''
     Credit to nf-core: https://github.com/nf-core/viralrecon

     Contribs from: https://github.com/drpatelh
                    https://github.com/saramonzon
                    https://github.com/svarona
    '''

    filename = os.path.splitext(os.path.basename(FileIn))[0]

    with open(RefIn,"r") as f: 
        contigName = []
        for line in f:
            if line.startswith(">"):
                contigName.append(line.split()[0][1:])

    if len(contigName) > 1:
        raise ValueError("Detected multiple contigs in reference file")

    header = ('##fileformat=VCFv4.2\n'
              '##source=iVar\n'
              '##contig=<ID='+contigName[0]+'>\n'
              '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'
              '##FILTER=<ID=PASS,Description="Result of p-value <= 0.05">\n'
              '##FILTER=<ID=FAIL,Description="Result of p-value > 0.05">\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
              '##FORMAT=<ID=REF_DP,Number=1,Type=Integer,Description="Depth of reference base">\n'
              '##FORMAT=<ID=REF_RV,Number=1,Type=Integer,Description="Depth of reference base on reverse reads">\n'
              '##FORMAT=<ID=REF_QUAL,Number=1,Type=Integer,Description="Mean quality of reference base">\n'
              '##FORMAT=<ID=ALT_DP,Number=1,Type=Integer,Description="Depth of alternate base">\n'
              '##FORMAT=<ID=ALT_RV,Number=1,Type=Integer,Description="Depth of alternate base on reverse reads">\n'
              '##FORMAT=<ID=ALT_QUAL,Number=1,Type=String,Description="Mean quality of alternate base">\n'
              '##FORMAT=<ID=ALT_FREQ,Number=1,Type=String,Description="Frequency of alternate base">\n')
    header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+filename+'\n'
    
    varList = []
    vcfString = []
    with open(FileIn) as f:
        for line in f:
            if not re.match("REGION",line):
                line = re.split("\t", line)
                CHROM=line[0]
                POS=line[1]
                ID='.'
                REF=line[2]
                ALT=line[3]
                var_type = 'SNP'
                if ALT[0] == '+':
                    ALT = REF + ALT[1:]
                    var_type = 'INS'
                elif ALT[0] == '-':
                    REF += ALT[1:]
                    ALT = line[2]
                    var_type = 'DEL'
                QUAL='.'
                pass_test=line[13]
                if pass_test == 'TRUE':
                    FILTER='PASS'
                else:
                    FILTER='FAIL'
                INFO='DP='+line[11]
                FORMAT='GT:REF_DP:REF_RV:REF_QUAL:ALT_DP:ALT_RV:ALT_QUAL:ALT_FREQ'
                SAMPLE='1:'+line[4]+':'+line[5]+':'+line[6]+':'+line[7]+':'+line[8]+':'+line[9]+':'+line[10]
                oline = CHROM+'\t'+POS+'\t'+ID+'\t'+REF+'\t'+ALT+'\t'+QUAL+'\t'+FILTER+'\t'+INFO+'\t'+FORMAT+'\t'+SAMPLE
                writeLine = True
                if (CHROM,POS,REF,ALT) in varList:
                    writeLine = False
                if re.match('^.*N+$', REF):
                    writeLine = False
                else:
                    varList.append((CHROM,POS,REF,ALT))
                if writeLine:
                    vcfString.append(oline)

    callsString = '\n'.join(vcfString)

    return header + callsString


def csq_annotate_vcf_string(vcfString,RefIn,GffIn):
    # bcftools csq -f fasta -g gff vcf
    p = run(['bcftools', 'csq', '-v', '0', '-f', RefIn, '-g', GffIn], stdout=PIPE,
            input=vcfString, encoding='ascii')

    return p.stdout

def extract_csq_info_from_vcf_string(csqVcf, minAF, minDP):

    v = io.StringIO(csqVcf)

    vcf_reader = vcf.Reader(v)
    bcsq_info = vcf_reader.infos['BCSQ']

    bcsq_keys_string = re.findall(r'Format: ?([\w+\|]+)', bcsq_info.desc)[0]
    bcsq_keys_string = re.sub(r'(\[|\]|\*)', '', bcsq_keys_string)
    bcsq_keys = [ x.lower() for x in bcsq_keys_string.split('|')]

    infos = []

    if vcf_reader.infos.get('SupportFraction'):
        vcf_type = "nanopolish"

    elif vcf_reader.formats.get('ALT_FREQ'):
        vcf_type = "ivar"

    elif vcf_reader.infos.get('AC'):
        vcf_type = "medaka"

    else:
        vcf_type = None
    

    for record in vcf_reader:

        if vcf_type == "nanopolish":
            if record.FILTER:
                continue

            if int(record.INFO['TotalReads']) < minDP:
                continue

            if float(record.INFO['SupportFraction']) < minAF:
                continue

        if vcf_type == "ivar":
            sampleid = record.samples[0].sample
            record_AF = float(record.genotype(sampleid)['ALT_FREQ'])

            if record.FILTER:
                continue

            if int(record.INFO['DP']) < minDP:
                continue

            if record_AF < minAF:
                continue

        if vcf_type == "medaka":

             if record.FILTER:
                 continue

             record_tot_DP = int(record.INFO['DP'])
             record_alt_DP = int(record.INFO['AC'][1])
             record_AF = record_alt_DP / record_tot_DP

             if record_tot_DP < minDP:
                 continue

             if record_AF < minAF:
                 continue

        variant_BCSQ = record.INFO.get('BCSQ')

        if variant_BCSQ:
            for variant in variant_BCSQ:
                variant_BCSQ_list = variant.split('|')

                if len(bcsq_keys) == len(variant_BCSQ_list):
                    variant_dict = dict(zip(bcsq_keys, variant_BCSQ_list))
                    if variant_dict not in infos:
                       infos.append(variant_dict)

    return infos


def get_variant_summary(info):

    sample_vars = []
 
    for variant in info:

        aa_r = re.compile("(?P<refpos>[0-9]+)(?P<refaa>[A-Z\*]+)*>*(?P<varpos>[0-9]+)*(?P<varaa>[A-Z\*]+)")
        aa_var = aa_r.match(variant['amino_acid_change']).groupdict()


        dna_r = re.compile("(?P<refpos>[0-9]+)(?P<refnucl>[A-Z\*]+)>(?P<varnucl>[A-Z\*]+)")
        dna_var = dna_r.match(variant['dna_change']).groupdict()

        if 'synonymous' not in variant['consequence'] and 'stop_retained' not in variant['consequence']:

            if 'frameshift' in variant['consequence']:

                complete_aa_variant_string = 'Frameshift.' + aa_var['refpos'] + aa_var['varaa']

            else:
                complete_aa_variant_string = aa_var['refaa'] + aa_var['refpos'] + aa_var['varaa']
  
        else:
            complete_aa_variant_string = 'Syn.' + aa_var['refpos'] + aa_var['varaa']
               

        complete_dna_variant_string = dna_var['refnucl'] + dna_var['refpos'] + dna_var['varnucl']

        variant_dict = { 'gene' : variant['gene'], 'aa_var' : complete_aa_variant_string, 'dna_var': complete_dna_variant_string }

        if variant_dict not in sample_vars:
            sample_vars.append(variant_dict)

    return sample_vars
    

def read_types_yaml(inFile):
    with open(inFile, 'r') as file:
        return yaml.full_load(file)

    
def type_vars_in_sample(types, sample_vars):

    types_assigned = []

    for typename, data in types.items():
        found_vars_from_type = []
        missing_vars_from_type = []
        additional_vars_from_type = sample_vars.copy()
        count_found = 0
        count_missing = 0

        # {gene, [var, var]}
        for gene in data['variants']:

            # { gene: gene, aa_var: var, dna_var: var }
            for sample_variant in sample_vars:

                if gene in sample_variant['gene']:

                    if sample_variant['aa_var'] in data['variants'][gene]:

                        data['variants'][gene].remove(sample_variant['aa_var'])
                        additional_vars_from_type.remove(sample_variant)

                        found_var_string = gene + "." + sample_variant['aa_var']
 
                        found_vars_from_type.append(found_var_string)

                        count_found += 1

   
            for missing_var in data['variants'][gene]:
                missing_var_string = gene + "." + missing_var
                if missing_var_string not in missing_vars_from_type:
                   missing_vars_from_type.append(missing_var_string)
                   count_missing += 1
        
        additional_vars_list = []

        for add_var in additional_vars_from_type:
            add_var_string = add_var['gene'] + '.' + add_var['aa_var']
            if not add_var_string in additional_vars_list:
                additional_vars_list.append(add_var_string)

        calc_coverage = count_found / ( count_missing + count_found )

        if calc_coverage >= data['coverage']:
            assigned_type = { 'type' : typename, 
                              'num_matching_vars' : count_found,
                              'num_missing_vars' : count_missing,
                              'type_coverage' : round(calc_coverage, 2),
                              'found_vars': ';'.join(found_vars_from_type),
                              'missing_vars': ';'.join(missing_vars_from_type),
                              'additional_vars' : ';'.join(additional_vars_list)}
            types_assigned.append(assigned_type)


    return types_assigned

def write_types_to_csv(types_assigned, sampleID, csvFileOut):

    fieldnames = list(types_assigned[0].keys())

    fieldnames.insert(0, 'sampleID')

    with open(csvFileOut, 'w', newline='') as csvfile:
        
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for row in types_assigned:
            row['sampleID'] = sampleID
            writer.writerow(row)

def read_vcf_to_vcf_string(FileIn):

    if FileIn.endswith('.gz'):
        with gzip.open(FileIn, 'rt') as f:
            vcfString = f.read()

    else:
        with open(FileIn, 'r') as f:
            vcfString = f.read()

    return vcfString

def write_csqAnnotatedVcfString_to_file(vcfOut, vcfString):
    with open(vcfOut, 'w') as f:
        f.write(vcfString)

def write_sample_vars_to_csv(summaryCsvOut, sampleID, sampleVars):
    fieldnames = list(sampleVars[0].keys())

    fieldnames.insert(0, 'sampleID')

    with open(summaryCsvOut, 'w', newline='') as csvfile:

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for row in sampleVars:
            row['sampleID'] = sampleID
            writer.writerow(row)

    

def main(args=None):
    args = parse_args(args)

    if args.TSV_IN:
        vcfString = ivar_variants_to_vcf_string(args.TSV_IN,args.REF_IN)

    elif args.VCF_IN:
        vcfString = read_vcf_to_vcf_string(args.VCF_IN)

    csqAnnotatedVcfString = csq_annotate_vcf_string(vcfString,args.REF_IN,args.GFF_IN)

    write_csqAnnotatedVcfString_to_file(args.VCF_OUT, csqAnnotatedVcfString)

    infos = extract_csq_info_from_vcf_string(csqAnnotatedVcfString,args.ALLELE_FREQ_THRESH, args.MIN_DEPTH)

    sample_vars = get_variant_summary(infos)

    if sample_vars:
        write_sample_vars_to_csv(args.SUMMARY_CSV_OUT, args.SAMPLE_ID, sample_vars)

    if args.YML_IN:
        types = read_types_yaml(args.YML_IN)
        types_assigned = type_vars_in_sample(types, sample_vars)

        if types_assigned:
            write_types_to_csv(types_assigned, args.SAMPLE_ID, args.TYPING_CSV_OUT)

if __name__ == '__main__':
    sys.exit(main())

