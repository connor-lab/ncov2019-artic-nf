#!/bin/bash
set -eo pipefail
export PATH=/opt/conda/bin:$PATH

# run current pull request code
singularity --version

# Clone variant_definitions repo
git clone https://github.com/phe-genomics/variant_definitions.git

# write test log as github Action artifact
echo "Nextflow run current PR in --nanopolish mode with typing.." >> artifacts/test_artifact.log
nextflow run ./main.nf \
       -profile singularity \
       --gb $PWD/typing/NC_045512.2.gb \
       --variant_definitions  $PWD/variant_definitions/variant_yaml \
       --nanopolish \
       --sequencing_summary $PWD/.github/data/nanopore/20200311_1427_X4_FAK72834_a3787181/sequencing_summary_FAK72834_298b7829.txt \
       --basecalled_fastq $PWD/.github/data/nanopore/20200311_1427_X4_FAK72834_a3787181/fastq_pass/ \
       --fast5_pass $PWD/.github/data/nanopore/20200311_1427_X4_FAK72834_a3787181/fast5_pass/ \
       --prefix 20200311_1427_X4_FAK72834_a3787181_typing
cp .nextflow.log artifacts/nanopolish_typing.nextflow.log
rm -rf results && rm -rf work && rm -rf .nextflow*

echo "Nextflow run current PR in --medaka mode with typing .." >> artifacts/test_artifact.log
nextflow run ./main.nf \
       -profile singularity \
       --medaka \
       --gb $PWD/typing/NC_045512.2.gb \
       --variant_definitions  $PWD/variant_definitions/variant_yaml \
       --basecalled_fastq $PWD/.github/data/nanopore/20200311_1427_X4_FAK72834_a3787181/fastq_pass/ \
       --prefix 20200311_1427_X4_FAK72834_a3787181
cp .nextflow.log artifacts/medaka_typing.nextflow.log
rm -rf results && rm -rf work && rm -rf .nextflow*

echo Nextflow run current PR in --illumina mode with typing.. >> artifacts/test_artifact.log
nextflow run ./main.nf \
       -profile singularity \
       --gb $PWD/typing/NC_045512.2.gb \
       --variant_definitions  $PWD/variant_definitions/variant_yaml \
       --directory $PWD/.github/data/fastqs/ \
       --illumina \
       --prefix test
cp .nextflow.log artifacts/illumina_typing.nextflow.log
rm -rf results && rm -rf work && rm -rf .nextflow*
