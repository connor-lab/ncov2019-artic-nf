# # nextflow run main.nf -profile docker --illumina --prefix "ctr26" --directory data --schemeRepoURL "--recurse-submodules https://github.com/artic-network/artic-ncov2019"

# nextflow run main.nf -profile docker --illumina --prefix "ena" --directory illumina --schemeRepoURL "https://github.com/ctr26/primer-schemes"

INPUT = "$1"

nextflow run main.nf -profile k8s --medaka --prefix "ena" --basecalled_fastq $INPUT

# nextflow run main.nf -profile docker --medaka --prefix "ena" --basecalled_fastq $INPUT