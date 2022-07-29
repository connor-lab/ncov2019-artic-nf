#!/usr/bin/env python
import os
import sys
import re

version_regex = {
    "artic": r"artic (\S+)",
    "bcftools": r"bcftools (\S+)",
    "bwa": r"Version: (\S+)",
    "constellations": r"constellations v(\S+)",
    "fastqc": r"FastQC v(\S+)",
    "freebayes": r"version:  v(\d\.\d\.\d+)",
    "gofasta": r"gofasta version (\S+)",
    "ivar": r"iVar version (\S+)",
    "mafft": r"v(\S+)",
    "medaka": r"medaka (\S+)",
    "minimap2": r"(\S+)",
    "multiqc": r"multiqc, version (\S+)",
    "muscle": r"MUSCLE (\S+)",
    "nanopolish": r"nanopolish version (\S+)",
    "nextclade": r"(\S+)",
    "pangolin-data": r"pangolin-data (\S+)",
    "pangolin": r"pangolin (\S+)",
    "picard": r"Version:(\S+)",
    "porechop": r"(\S+)",
    "python": r"Python (\S+)",
    "samtools": r"samtools (\S+)",
    "scorpio": r"scorpio (\S+)",
    "snakemake": r"(\S+)",
    "trim_galore": r"version (\S+)",
    "usher": r"UShER \((\S+)\)",
}
# Create csv file for all versions
output = sys.argv[1]
with open(output, "w") as out:
    for tool, regex in version_regex.items():
        # Search for version number
        version_file = "version_{}.txt".format(tool)
        if os.path.isfile(version_file):
            with open(version_file) as f:
                contents = f.read()
                match = re.search(regex, contents)
                if match:
                    # Add version number to output
                    out.write("{},{}\n".format(tool, match.group(1)))
