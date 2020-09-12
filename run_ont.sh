# nextflow run main.nf -profile docker --illumina --prefix "ctr26" --directory data --schemeRepoURL "--recurse-submodules https://github.com/artic-network/artic-ncov2019"

nextflow run main.nf -profile docker --illumina --prefix "ena" --directory illumina --schemeRepoURL "https://github.com/ctr26/primer-schemes"

nextflow run main.nf -profile docker --nanopolish --prefix "ctr26" --basecalled_fastq data --sequencing_summary sequencing_summary.txt --schemeRepoURL "https://github.com/ctr26/primer-schemes"