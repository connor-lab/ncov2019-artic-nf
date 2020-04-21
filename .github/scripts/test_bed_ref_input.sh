#!/bin/bash
set -eo pipefail
export PATH=/opt/conda/bin:$PATH

echo test --bed and --ref inputs.. >> artifacts/test_artifact.log

export REF_FILE=$(find work_singularity_profile | grep -v scheme | grep 'ref.fa$' | xargs readlink | grep ncov2019-artic-nf | head -n1 | sed s'/ncov2019-artic-nf\/work\//ncov2019-artic-nf\/work_singularity_profile\//'g)
export BED_FILE=$(find work_singularity_profile | grep -v scheme | grep 'bed$' | xargs readlink | grep ncov2019-artic-nf | head -n1 | sed s'/ncov2019-artic-nf\/work\//ncov2019-artic-nf\/work_singularity_profile\//'g)
echo ref file: $REF_FILE
echo bed file: $BED_FILE
echo ref file: $REF_FILE >> artifacts/test_artifact.log
echo bed file: $BED_FILE >> artifacts/test_artifact.log

# run current pull request code
singularity --version
echo Nextflow run --illumina mode with --ref, --bed .. >> artifacts/test_artifact.log
NXF_VER=20.03.0-edge nextflow run ./main.nf \
       -profile singularity \
       --ref $REF_FILE \
       --bed $BED_FILE \
       --directory $PWD/.github/data/fastqs/ \
       --illumina \
       --prefix test
cp .nextflow.log artifacts/ref_bed.nextflow.log

# check that git clone did not ran:
find work -name .command.sh \
    | xargs cat | grep 'git clone' \
    && bash -c "echo test failed\: git clone ran and shouldn\'t have in --ref, --bed mode >> artifacts/test_artifact.log && exit 1"
echo "ran with --bed and --ref: did NOT use git clone" >> artifacts/test_artifact.log

echo Compare ouputs of current PR with or without using --bed/--ref.. >> artifacts/test_artifact.log
find results results_singularity_profile \
     -name "test.qc.csv" \
     -o -name "*.fq.gz" \
     -o -name "*.bam" \
     -o -name "scheme" | xargs rm -rf
git diff --stat --no-index results results_singularity_profile > diffs.txt
if [ -s diffs.txt ]
then
  echo "test failed: differences found for pull request with or without --ref and --bed" >> artifacts/test_artifact.log
  echo "see diffs_bed_ref.txt" >> artifacts/test_artifact.log
  cp diffs.txt artifacts/diffs_bed_ref.txt  
  exit 1
else
  echo no differences found for pull request with or without --ref and --bed  >> artifacts/test_artifact.log 
fi

# clean-up for following tests
rm -rf results && rm -rf work && rm -rf .nextflow*
