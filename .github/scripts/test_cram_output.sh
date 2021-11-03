#!/bin/bash
set -eo pipefail
export PATH=/opt/conda/bin:$PATH

echo test --outCram and --ref inputs >> artifacts/test_artifact.log

export REF_FILE=$(find work_singularity_profile | grep -v scheme | grep 'ref.fa$' | xargs readlink | grep ncov2019-artic-nf | head -n1 | sed s'/ncov2019-artic-nf\/work\//ncov2019-artic-nf\/work_singularity_profile\//'g)
export BED_FILE=$(find work_singularity_profile | grep -v scheme | grep 'bed$' | xargs readlink | grep ncov2019-artic-nf | head -n1 | sed s'/ncov2019-artic-nf\/work\//ncov2019-artic-nf\/work_singularity_profile\//'g)
echo ref file: $REF_FILE
echo bed file: $BED_FILE

# run current pull request code
# --sanger profile: there are only 2 available cpus in the github runner execution
sed -i s'/cpus = 4/cpus = 2/'g conf/coguk/sanger.config
singularity --version
echo Nextflow run --illumina mode with --ref, --bed --cram and outCram... >> artifacts/test_artifact.log
NXF_VER=21.04.0 nextflow run main.nf \
       -profile sanger,singularity \
       --ref $REF_FILE \
       --bed $BED_FILE \
       --cram \
       --outCram \
       --directory .github/data/crams/ \
       --illumina \
       --prefix "test"
cp .nextflow.log artifacts/Outcram_ref_bed.nextflow.log

# get first cram file that macthes in given folder 
cram_exists=($(find results/*_bamToCram -name '*.cram'))

if [ -f "$cram_exists" ]; then
      echo Cram file generated successfully >> artifacts/test_artifact.log
      # clean tests
      rm -rf results && rm -rf work && rm -rf .nextflow*
      exit 0
else
      echo Cram file generation failed >> artifacts/test_artifact.log
      # clean-up for following tests
      echo cleaning tests
      rm -rf results && rm -rf work && rm -rf .nextflow*
      exit 1
fi
