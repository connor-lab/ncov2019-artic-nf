# nextflow run connor-lab/ncov2019-artic-nf -profile docker --nanopolish --prefix "ctr26" --basecalled_fastq data --sequencing_summary sequencing_summary.txt --schemeRepoURL "--recurse-submodules https://github.com/artic-network/artic-ncov2019"


# nextflow run main.nf -profile docker --illumina --prefix "ctr26" --directory data --schemeRepoURL "--recurse-submodules https://github.com/artic-network/artic-ncov2019"


INPUT = "$1"

nextflow run main.nf -profile k8s --illumina --prefix "ena" --directory $INPUT
