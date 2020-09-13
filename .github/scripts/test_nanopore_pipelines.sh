#!/bin/bash
set -eo pipefail
export PATH=/opt/conda/bin:$PATH

# run current pull request code
singularity --version
# write test log as github Action artifact
echo "Nextflow run current PR in --nanopolish mode (no barcodes).." >> artifacts/test_artifact.log
NXF_VER=20.03.0-edge nextflow run ./main.nf \
       -profile singularity \
       --nanopolish \
       --sequencing_summary $PWD/.github/data/nanopore/20200311_1427_X4_FAK72834_a3787181/sequencing_summary_FAK72834_298b7829.txt \
       --basecalled_fastq $PWD/.github/data/nanopore/20200311_1427_X4_FAK72834_a3787181/fastq_pass/ \
       --fast5_pass $PWD/.github/data/nanopore/20200311_1427_X4_FAK72834_a3787181/fast5_pass/ \
       --prefix 20200311_1427_X4_FAK72834_a3787181 
cp .nextflow.log artifacts/nanopolish.nextflow.log
rm -rf results && rm -rf work && rm -rf .nextflow*

echo "Nextflow run current PR in --nanopolish mode (barcodes).." >> artifacts/test_artifact.log
NXF_VER=20.03.0-edge nextflow run ./main.nf \
       -profile singularity \
       --nanopolish \
       --sequencing_summary $PWD/.github/data/nanopore/20200311_1427_X1_FAK72834_a3787181/sequencing_summary_FAK72834_298b7829.txt \
       --basecalled_fastq $PWD/.github/data/nanopore/20200311_1427_X1_FAK72834_a3787181/fastq_pass/ \
       --fast5_pass $PWD/.github/data/nanopore/20200311_1427_X1_FAK72834_a3787181/fast5_pass/ \
       --prefix 20200311_1427_X1_FAK72834_a3787181 
cp .nextflow.log artifacts/nanopolish_barcodes.nextflow.log
rm -rf results && rm -rf work && rm -rf .nextflow*

echo "Nextflow run current PR in --medaka mode (no barcodes).." >> artifacts/test_artifact.log
NXF_VER=20.03.0-edge nextflow run ./main.nf \
       -profile singularity \
       --medaka \
       --basecalled_fastq $PWD/.github/data/nanopore/20200311_1427_X4_FAK72834_a3787181/fastq_pass/ \
       --prefix 20200311_1427_X4_FAK72834_a3787181
cp .nextflow.log artifacts/medaka.nextflow.log
rm -rf results && rm -rf work && rm -rf .nextflow*

echo "Nextflow run current PR in --medaka mode (barcodes).." >> artifacts/test_artifact.log
NXF_VER=20.03.0-edge nextflow run ./main.nf \
       -profile singularity \
       --medaka \
       --basecalled_fastq $PWD/.github/data/nanopore/20200311_1427_X1_FAK72834_a3787181/fastq_pass/ \
       --prefix 20200311_1427_X1_FAK72834_a3787181

cp .nextflow.log artifacts/medaka_barcodes.nextflow.log
rm -rf results && rm -rf work && rm -rf .nextflow*
