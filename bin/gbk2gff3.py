#!/usr/bin/env python3

# based on work by Damien Farrell https://dmnfarrell.github.io/bioinformatics/bcftools-csq-gff-format
import sys

def GFF_bcftools_format(in_handle, out_handle):
    """Convert a bacterial genbank file from NCBI to a GFF3 format that can be used in bcftools csq.
    see https://github.com/samtools/bcftools/blob/develop/doc/bcftools.txt#L1066-L1098.
    Args:
        in_file: genbank file
        out_file: name of GFF file
    """

    from BCBio import GFF
    #in_handle = open(in_file)
    #out_handle = open(out_file, "w")
    from Bio.SeqFeature import SeqFeature
    from Bio.SeqFeature import FeatureLocation
    from copy import copy, deepcopy
    from Bio import SeqIO

    for record in SeqIO.parse(in_handle, "genbank"):
        #make a copy of the record as we will be changing it during the loop
        new = copy(record)
        new.features = []
        #loop over all features
        for feat in record.features:
            q = feat.qualifiers
            #remove some unecessary qualifiers
            for label in ['note','translation','product','experiment']:
                if label in q:
                    del q[label]


            if(feat.type == "CDS"):
                #use the CDS feature to create the new lines
                tag = q['gene'][0] #q['locus_tag'][0]
                protein_id = q['protein_id'][0]
                q['ID'] = 'CDS:%s' %protein_id
                q['biotype'] = 'protein_coding'

                for i, new_loc in enumerate(feat.location.parts if(hasattr(feat.location, "parts")) else (feat.location,)):
                        new_feat = deepcopy(feat)
                        tr_id = 'transcript:%s' %(protein_id+"_"+str(i))
                        new_feat.qualifiers['Parent'] = tr_id
                        new_feat.location = new_loc
                        new.features.append(new_feat)
                        #create mRNA feature
                        m = SeqFeature(feat.location, type='mRNA',strand=feat.strand)
                        q2 = m.qualifiers
                        q2['ID'] = tr_id
                        q2['Parent'] = 'gene:%s' %tag
                        q2['biotype'] = 'protein_coding'
                        new.features.append(m)

            elif(feat.type == "gene"):
                tag = q['gene'][0]
                #edit the gene feature
                q=feat.qualifiers
                q['ID'] = 'gene:%s' %tag
                q['biotype'] = 'protein_coding'
                q['Name'] = q['gene']
                new.features.append(feat)
        #write the new features to a GFF
        GFF.write([new], out_handle)
        return

if __name__ == "__main__":
   GFF_bcftools_format(sys.stdin, sys.stdout)
