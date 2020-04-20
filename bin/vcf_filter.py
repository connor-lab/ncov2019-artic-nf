#!/usr/bin/env python

import vcf
import sys
import subprocess
from operator import attrgetter
from collections import defaultdict
#from .vcftagprimersites import read_bed_file

class NanoporeFilter:
    def __init__(self, *args, **kwargs):
        pass

    def check_filter(self, v):
        total_reads = float(v.INFO['TotalReads'])
        qual = v.QUAL
        strandbias = float(v.INFO['StrandFisherTest'])

        if qual / total_reads <= 4:
            return False

        strand_fraction_by_strand = v.INFO['SupportFractionByStrand']
        if float(strand_fraction_by_strand[0]) < 0.5: 
            return False

        if float(strand_fraction_by_strand[1]) < 0.5:
            return False

        #if strandbias >= 100:
        #    return False

        if total_reads < 20:
            return False

        return True

class MedakaFilter:
    def __init__(self, *args, **kwargs):
        pass

    def check_filter(self, v):
        if v.num_het:
            return False
        return True

class LongshotFilter:
    def __init__(self, *args, **kwargs):
        pass

    def check_filter(self, v):
        depth = v.INFO['DP']
        if depth < 20:
            return False

        if v.num_het:
            return False
        return True

class IlluminaFilter:
        def __init__(self, *args, **kwargs):
            # Option to have different filters for variants overlapping primer sites and those not
            self.overlapping_primer_site = kwargs.get('run_primer_site_filtering')

        def check_filter(self, v):
            if v.is_indel:
                return self.check_filter_indels(v)
            else:
                return self.check_filter_snps(v)
      
        def check_filter_snps(self, v):
            # Filters for snps called with lofreq
            if self.overlapping_primer_site:
                return NotImplementedError("Primer site specific filtering not yet implemented for indels")
            else:
                # TODO Other filters here
                # Consider annotating as fail and returning as true, will be removed from vcf later
                # Return as fail needed to populate fail.vcf- consider impact of
                pass
            return True
        
        def check_filter_indels(self, v):
            # Filters for indels called with lofreq
            if self.overlapping_primer_site:
                return NotImplementedError("Primer site specific filtering not yet implemented for indels")
            else:
                # TODO Other filters here
                # Consider annotating as fail and returning as true, will be removed from vcf later
                # Return as fail needed to populate fail.vcf- consider impact of
                pass
            return True           

def go(args):
    vcf_reader = vcf.Reader(filename=args.inputvcf)
    vcf_writer = vcf.Writer(open(args.outputvcf, 'w'), vcf_reader)
    vcf_writer_filtered = vcf.Writer(open(args.output_fail_vcf, 'w'), vcf_reader)
    if args.nanopolish:
        filter = NanoporeFilter()
    elif args.medaka:
        filter = MedakaFilter()
    elif args.longshot:
        filter = LongshotFilter()
    elif args.illumina:
        filter = IlluminaFilter(run_primer_site_filtering=args.primer_site)
    else:
        print("Please specify a VCF type, i.e. --nanopolish, --medaka or --illumina. For illumina \
                the optional argument --primer_site is available to run different options for primer site \
                    variants. Requires a vcf pre-filtered to include only primer site variants.\n")
        raise SystemExit

    for v in vcf_reader:
        if filter.check_filter(v):
            vcf_writer.write_record(v)
        else:
            vcf_writer_filtered.write_record(v)

def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--nanopolish', action='store_true')
    parser.add_argument('--medaka', action='store_true')
    parser.add_argument('--longshot', action='store_true')
    parser.add_argument('--illumina', action='store_true')
    parser.add_argument('--primer_site', action='store_true', default=False)
    parser.add_argument('inputvcf')
    parser.add_argument('outputvcf')
    parser.add_argument('output_fail_vcf')

    args = parser.parse_args()
    
    if getattr(args, 'primer_site', None) and not getattr(args, 'illumina', None):
        parser.print_help()
        parser.exit(1, "Primer site specific vcf option is only valid for illumina filtering. \
                        Only use --strand_bias with --illumina option.")

    go(args)

if __name__ == "__main__":
    main()