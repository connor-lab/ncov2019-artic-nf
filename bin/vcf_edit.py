#!/usr/bin/env python

import vcf
import argparse
import textwrap
import os
from collections import namedtuple

# Dictionary of IUPAC codes
iupac_dict = {('A', 'G'): 'R', ('C', 'T'): 'Y', ('C', 'G'): 'S', ('A', 'T'): 'W', ('G', 'T'): 'K', ('A', 'C'): 'M',
              ('C', 'G', 'T'): 'B', ('A', 'G', 'T'): 'D', ('A', 'C', 'T'): 'H', ('A', 'C', 'G'): 'V',
              ('A', 'C', 'G', 'T'): 'N', 'uncertainty': 'N'}

# This will use a lot of memory for large vcfs

# Default values
default_minimum_depth = 20
minimum_allele_support_proportion = 0.25
confident_allele_support_proportion = 0.75


def get_args():
    argument_parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=textwrap.dedent(
            '''
            summary:
            Describe software!!!
            '''))

    # Add arguments
    # REQUIRED: Path to config file containing authorisation information
    argument_parser.add_argument(
        '-v', '--path_to_vcf', action='store',
        help=textwrap.dedent(
            '''
            File path to vcf file to be edited. REQUIRED.
            '''))
    argument_parser.add_argument(
        '-c', '--caller', action='store',
        help=textwrap.dedent(
            '''
             Specify variant caller that has been used to generate the VCF. Currently supported
             options are bcftools and lofreq. REQUIRED.
            '''))

   # OPTIONAL:
    optional_options = argument_parser.add_mutually_exclusive_group()
    optional_options.add_argument(
        '-o', '--output_file', action='store', default=False,
        help=textwrap.dedent(
            '''
            Add output filename. Default is {input filename}_edited.vcf. OPTIONAL.
            '''))
    optional_options.add_argument(
        '-d', '--min_depth', action='store', default=False,
        help=textwrap.dedent(
            '''
            Add user configurable minimum depth. Bases below this depth will be set to N Default value 20. OPTIONAL.
            '''))

    return argument_parser.parse_args()


def write_vcf(records_to_write, read_vcf, output_vcf="output.vcf"):
    vcf_writer = vcf.Writer(open(output_vcf, 'w'), read_vcf)
    for record in records_to_write:
        vcf_writer.write_record(record)


def apply_iupac(bases_dict):
    # Order in alphabetical order for lookup
    base_list = []
    [base_list.append(b) for b in bases_dict.keys()]
    base_list.sort()
    bases = tuple(base_list)
    find_iupac = bases
    called = iupac_dict.get(find_iupac)
    if not called:
        # Combination of nucleotides not in dictionary- call uncertain- set to N
        called = iupac_dict.get('uncertainty')
    return called


def format_data_format_creator(format_dict, base_pos, iupac_uncertainty=False):
    format_titles = []
    format_values = []
    for format_title, format_value in format_dict.items():
        # Only ever one entry for GT, not one entry per allele as for others
        format_titles.append(format_title)
        try:
            format_values.append(format_value[base_pos])
        except (IndexError, TypeError):
            # Single value only
            format_values.append(format_value)
    ca = namedtuple('CallData', format_titles)
    caller_data = ca(*format_values)
    return caller_data


def reads_to_proportions(total_depth, depths_per_base):
    proportions = []
    try:
        for base_depth in depths_per_base:
            proportions.append(base_depth/total_depth)
    except TypeError:
        proportions.append(depths_per_base/total_depth)
    return proportions


def get_proportions(info_field):
    try:
        proportions = info_field.get('AF')
    except KeyError:
        raise KeyError("No allele frequency annotation AF found")
    try:
        # For multiple alleles
        proportions = [x for x in proportions]
    except TypeError:
        # For single alleles
        proportions = [proportions]
    return proportions


def parse_out_sample_format(format_field, call_data):
    sample_format_data_dict = {}
    format_field_annotations = format_field.split(':')
    for i, annotation in enumerate(format_field_annotations):
        sample_format_data_dict[annotation] = call_data[i]
    return sample_format_data_dict


def get_total_high_qual_bases_at_site(info_field):
    try:
        high_quality_bases = sum(info_field['DP4'])
    except TypeError:
        try:
            high_quality_bases = int(info_field['DP4'])
        except (TypeError, ValueError):
            raise ValueError(f"Total number of high quality bases could not be determined from DP4 value. Proportion"
                             f" of bases supporting each read cannot be determined.")
    return high_quality_bases


def get_total_supporting_reads(reads):
    try:
        high_quality_bases = sum(reads)
    except TypeError:
        try:
            high_quality_bases = int(reads)
        except (TypeError, ValueError):
            raise ValueError(f"Total number of high quality bases could not be determined from DP4 value. Proportion"
                             f" of bases supporting each read cannot be determined.")
    return high_quality_bases


def load_vcf(vcf_file_path):
    try:
        read_vcf = vcf.Reader(open(vcf_file_path, 'r'))
    except TypeError:
        raise FileNotFoundError("VCF file could not be loaded. Check path to VCF has been correctly entered on"
                                " the command line (flag -v) and that it is a VCF file.")
    # Update vcf metadata to add processing
    read_vcf.metadata.update({"vcf_edit": ['VCF edited to apply IUPAC uncertainty codes to ALT snp calls.']})
    return read_vcf


def parse_vcf():
    return None


def parse_vcf_lofreq(vcf, minimum_depth):
    vcf_record = []
    # Update ALT base for SNPs
    # Do not process indel calls
    for record in vcf:
        # Don't process indels
        # Inbuild is-indel function finds reference calls (no alt) to be indels- bespoke identify indels
        is_indel = record.INFO.get('INDEL', False)
        # Do not process indels further
        if is_indel:
            # Do not propagate indels called through this workflow
            #continue
            # Output indels with allele frequency > 0.75
            indel_frequency = get_proportions(record.INFO)
            if len(indel_frequency) == 1:
                indel_frequency = indel_frequency[0]
            else:
                raise NotImplementedError("Currently no support for multiallelic indels")
            if indel_frequency > 0.75: #TODO Hardcoded temporily for working
                vcf_record.append(record)
            else:
                # Add filter
                record.FILTER.append("lowaf_indel")
                vcf_record.append(record)
            # Move on to next record (important due to processing for ref sites below)
            continue

        # Process all sites that are single nucleotides
        # Use DP for depth filtering- consistency with calculation for REF sites (if not emitted by caller)
        '''
        try:
            total_bases_at_site = record.INFO['DP']
        except KeyError:
            raise KeyError(f"Could not determine read coverage from this VCF at site {record.POS}. No DP annotation.")
        '''

        # Handle cases where there is only one DP4 value- calculate total reads supporting high quality bases using DP4
        #total_high_qual_bases_at_site = get_total_high_qual_bases_at_site(record.INFO)

        # Sites with minimum depth of coverage over all alleles at below threshold value are uncertain and set to N
        '''
        if total_bases_at_site < minimum_depth:
            record.ALT = iupac_dict.get('uncertainty')
            vcf_record.append(record)
            continue
        '''

        # Calculate high quality depth at site
        # Do not process reference where no alt was called sites further (note, these sites have no PL annotation)
        if record.is_indel:
            vcf_record.append(record)
            continue

        # Multisample vcfs not supported- TODO if needed
        if len(record.samples) > 1:
            raise NotImplementedError(f"Multisample VCFs not currently supported. "
                                      f" Only one sample per record supported.")

        # For variant sites passing minimum depth threshold, obtain proportion of reads supporting each base
        # Do not add an alt base if there is not one (and it is therefore None)
        all_bases = [str(base) for base in record.ALT if base is not None]
        all_bases.insert(0, record.REF)
        alt_bases = [str(base) for base in record.ALT if base is not None]

        # No FORMAT field in lofreq VCF
        # Filtering on proportion of supporting bases
        bases_to_process = {}
        iupac_uncertainty = False
        # Obtain proportion of reads supporting each base from INFO field
        supporting_proportions = get_proportions(record.INFO)
        running_total_low_alts = 0
        for i, b in enumerate(alt_bases):
            # If allele frequency is <0.25 pass through reference base- unknown AF to support REF base- set to .
            if supporting_proportions[i] < minimum_allele_support_proportion:
                record.INFO['AF'] = '.'
                record.ALT = record.REF

            # Pass through only bases representing 0.75 or more of the total reads
            elif supporting_proportions[i] > confident_allele_support_proportion:
                # Only one base can meet the above criteria
                # Edit INFO fields to remove proportions where bases have been removed (no FORMAT field)
                # Update INFO entries
                record.INFO['AF'] = supporting_proportions[i]

            # Sites with bases with support from >0.25 but <0.75 of the total reads are uncertain
            elif minimum_allele_support_proportion < supporting_proportions[i] < confident_allele_support_proportion:
                bases_to_process[b] = supporting_proportions[i]
                iupac_uncertainty = True

            # Count support for any other alleles with AF <0.25 (running total) at this site
            else:
                running_total_low_alts += supporting_proportions[i]

        # Process uncertainty sites
        if iupac_uncertainty:
            # Calculate allele frequency supporting any of the bases included in the iupac code- update
            # Case where there are multiple ALT alleles both with AF >0.25 and <0.75
            # If the support for all ALTs is <0.75 (including low support fraction), include REF (as it is >0.25)
            if (sum(bases_to_process.values()) + running_total_low_alts) < confident_allele_support_proportion:
                # Update bases to process dict to ensure REF is included in IUPAC (unknown support for call)
                bases_to_process[record.REF] = '.'
                record.INFO['AF'] = '.'
            else:
                # Edit INFO field to supply evidence for ambiguity base (summation of frequencies supporting)
                record.INFO['AF'] = sum(bases_to_process.values())
            record.ALT = apply_iupac(bases_to_process)
        # Update vcf records
        vcf_record.append(record)
    return vcf_record


def parse_vcf_bcftools(vcf, minimum_depth):
    vcf_record = []
    # Update ALT base for SNPs
    # Do not process indel calls
    for record in vcf:
        # Don't process indels
        # Inbuild is-indel function finds reference calls (no alt) to be indels- bespoke identify indels
        is_indel = record.INFO.get('INDEL', False)
        # Do not process indels further
        if is_indel:
            # Do not propagate indels called through this workflow
            continue

        # Process all sites that are single nucleotides
        # Use DP for depth filtering- consistency with calculation for REF sites (if not emitted by caller)
        '''
        try:
            total_bases_at_site = record.INFO['DP']
        except KeyError:
            raise KeyError(f"Could not determine read coverage from this VCF at site {record.POS}. No DP annotation.")
        '''

        # Handle cases where there is only one DP4 value- calculate total reads supporting high quality bases using DP4
        total_high_qual_bases_at_site = get_total_high_qual_bases_at_site(record.INFO)

        # Sites with minimum depth of coverage below threshold value are uncertain and set to N
        '''
        if total_bases_at_site < minimum_depth:
            record.ALT = iupac_dict.get('uncertainty')
            vcf_record.append(record)
            continue
        '''

        # Calculate high quality depth at site
        # Do not process reference where no alt was called sites further (note, these sites have no PL annotation)
        if record.is_indel:
            vcf_record.append(record)
            continue

        # Multisample vcfs not supported- TODO if needed
        if len(record.samples) > 1:
            raise NotImplementedError(f"Multisample VCFs not currently supported. "
                                      f" Only one sample per record supported.")

        # For variant sites passing minimum depth threshold, obtain proportion of reads supporting each base
        # Do not add an alt base if there is not one (and it is therefore None)
        all_bases = [str(base) for base in record.ALT if base is not None]
        all_bases.insert(0, record.REF)

        # Obtain correct FORMAT data for site (correct present fields associated with correct entry)
        format_field_dictionary = parse_out_sample_format(record.FORMAT, record.samples[0].data)
        try:
            supporting_reads = format_field_dictionary['AD']
        except KeyError:
            raise KeyError(f"Allele depth annotation not found for VCF. Unable to determine supporting reads per"
                           f" allele. Cannot filter VCFs without the AD annotation.")
        # Sanity Check
        total_supporting_reads = get_total_supporting_reads(supporting_reads)
        if total_high_qual_bases_at_site != total_supporting_reads:
            raise ValueError(f"Total number of high quality bases could not be unambiguously determined. Proportion"
                             f" of bases supporting each read cannot be determined.")

        # Obtain proportion of reads supporting each base from numbers
        supporting_proportions = reads_to_proportions(total_high_qual_bases_at_site, supporting_reads)

        # Filtering on proportion of supporting bases
        bases_to_process = {}
        iupac_uncertainty = False
        # Get AD from INFO field, once per site
        site_info_allele_depth = record.INFO.get('AD')
        for i, b in enumerate(all_bases):
            # Obtain allele depth (supporting reads) for this base
            info_allele_depth = site_info_allele_depth[i]
            format_allele_depth = supporting_reads[i]
            # Sanity check- both allele depths same
            if info_allele_depth != format_allele_depth:
                raise ValueError(
                    f"Total number of high quality bases could not be unambiguously determined. Proportion"
                    f" of bases supporting each read cannot be determined.")

            # Pass through only bases representing 0.75 or more of the total reads
            if supporting_proportions[i] > confident_allele_support_proportion:
                # Only one base can meet the above criteria
                # Edit INFO and FORMAT fields to remove reads where bases have been removed
                # Update INFO entries
                record.INFO['AD'] = format_allele_depth

                # Update FORMAT entries
                caller_data = format_data_format_creator(format_field_dictionary, i)
                record.samples[0].data = caller_data

                # If final call is the same as the reference base alt should be set to .
                if b == record.REF:
                    record.ALT = '.'
                else:
                    record.ALT = b

            # Sites with bases with support from >0.25 but <0.75 of the total reads are uncertain
            elif minimum_allele_support_proportion < supporting_proportions[i] < confident_allele_support_proportion:
                bases_to_process[b] = supporting_reads[i]
                iupac_uncertainty = True

        # Process uncertainty sites
        if iupac_uncertainty:
            record.ALT = apply_iupac(bases_to_process)
            # Calculate depth of reads supporting any of the bases included in the iupac code- update
            updated_format_allele_depth = (sum(bases_to_process.values()))

            # Edit INFO and FORMAT fields to supply evidence for ambiguity base (summation of reads supporting)
            record.INFO['AD'] = updated_format_allele_depth

            # Update FORMAT entries- retain unknown genotype and allele depth information only- single allele only (0)
            record.FORMAT = 'GT:AD'
            new_format_field_dictionary = {'GT': '.', 'AD': updated_format_allele_depth}
            caller_data = format_data_format_creator(new_format_field_dictionary, 0, iupac_uncertainty=True)
            record.samples[0].data = caller_data

        # Update vcf records
        vcf_record.append(record)
    return vcf_record


def main():
    __version__ = '0.0.1'
    __updated__ = '30/03/2020'
    args = get_args()
    # Load VCF file
    vcf_data = load_vcf(args.path_to_vcf)
    # Set output filename
    if args.output_file:
        outfile = args.output_file
    else:
        outfile = f"{os.path.splitext(args.path_to_vcf)[0]}_edited.vcf"
    # Set minimum depth to call a base (below this it is set to N)
    if args.min_depth:
        min_depth = args.min_depth
    else:
        min_depth = default_minimum_depth

    # Update/Filter vcf data where required
    if args.caller == 'bcftools':
        updated_vcf_data = parse_vcf_bcftools(vcf_data, min_depth)
    elif args.caller == 'lofreq':
        updated_vcf_data = parse_vcf_lofreq(vcf_data, min_depth)
    else:
        raise NotImplementedError(f"No current support for variant callers other than bcftools and lofreq")

    # Write updated vcf file
    write_vcf(updated_vcf_data, vcf_data, outfile)


if __name__ == '__main__':
    main()
