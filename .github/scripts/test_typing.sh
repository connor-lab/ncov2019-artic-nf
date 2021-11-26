#!/bin/bash
set -eo pipefail
export PATH=/opt/conda/bin:$PATH

# run current pull request code
singularity --version
# write test log as github Action artifact
echo "Nextflow run current PR in --nanopolish mode with typing.." >> artifacts/test_artifact.log
NXF_VER=21.04.0 nextflow run main.nf \
       -profile singularity \
       --gff typing/MN908947.3.gff \
       --yaml typing/SARS-CoV-2.types.yaml \
       --nanopolish \
       --sequencing_summary .github/data/nanopore/20200311_1427_X4_FAK72834_a3787181/sequencing_summary_FAK72834_298b7829.txt \
       --basecalled_fastq .github/data/nanopore/20200311_1427_X4_FAK72834_a3787181/fastq_pass/ \
       --fast5_pass .github/data/nanopore/20200311_1427_X4_FAK72834_a3787181/fast5_pass/ \
       --prefix "20200311_1427_X4_FAK72834_a3787181_typing"
cp .nextflow.log artifacts/nanopolish_typing.nextflow.log
rm -rf results && rm -rf work && rm -rf .nextflow*

echo "Nextflow run current PR in --medaka mode with typing .." >> artifacts/test_artifact.log
NXF_VER=21.04.0 nextflow run main.nf \
       -profile singularity \
       --medaka \
       --gff typing/MN908947.3.gff \
       --yaml typing/SARS-CoV-2.types.yaml \
       --basecalled_fastq .github/data/nanopore/20200311_1427_X4_FAK72834_a3787181/fastq_pass/ \
       --prefix "20200311_1427_X4_FAK72834_a3787181"
cp .nextflow.log artifacts/medaka_typing.nextflow.log
rm -rf results && rm -rf work && rm -rf .nextflow*

echo Nextflow run current PR in --illumina mode with typing.. >> artifacts/test_artifact.log
NXF_VER=21.04.0 nextflow run main.nf \
       -profile singularity \
       --gff typing/MN908947.3.gff \
       --yaml typing/SARS-CoV-2.types.yaml \
       --directory .github/data/fastqs/ \
       --illumina \
       --prefix "test"
cp .nextflow.log artifacts/illumina_typing.nextflow.log
rm -rf results && rm -rf work && rm -rf .nextflow*
