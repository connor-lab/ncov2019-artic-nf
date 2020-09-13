# # nextflow run main.nf -profile docker --illumina --prefix "ctr26" --directory data --schemeRepoURL "--recurse-submodules https://github.com/artic-network/artic-ncov2019"

# nextflow run main.nf -profile docker --illumina --prefix "ena" --directory illumina --schemeRepoURL "https://github.com/ctr26/primer-schemes"

FILENAME = "$1"

INPUT = "/data/${FILENAME}_input/"
OUTDIR="results"

nextflow run ncov2019-artic-nf/main.nf -profile docker --nanopolish --prefix "ena" --basecalled_fastq $INPUT --fast5-directory $INPUT --sequencing_summary sequencing_summary.txt --schemeRepoURL "https://github.com/ctr26/primer-schemes"