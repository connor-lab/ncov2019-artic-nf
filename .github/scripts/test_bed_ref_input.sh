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
NXF_VER=21.04.0 nextflow run main.nf \
       -profile singularity \
       --ref $REF_FILE \
       --bed $BED_FILE \
       --directory .github/data/fastqs/ \
       --illumina \
       --prefix "test"
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

rm -rf results && rm -rf work && rm -rf .nextflow*

# check it works if the ref is a symlink to a file with no index files next to it
export REAL_REF=$REF_FILE
export REF_FILE=$(find work_singularity_profile -type f -name "ref.fa.bwt" | xargs dirname | xargs realpath )/ref.fa
echo ref file: $REF_FILE
echo ref file: $REF_FILE >> artifacts/test_artifact.log
# ($REF_FILE is a symlink, but since we changed the directory name, it is the wrong symlink; recreate it)
rm $REF_FILE
ln -s $REAL_REF $REF_FILE

echo Nextflow run --illumina mode with symlinked --ref, --bed .. >> artifacts/test_artifact.log
NXF_VER=21.04.0 nextflow run main.nf \
       -profile singularity \
       --ref $REF_FILE \
       --bed $BED_FILE \
       --directory .github/data/fastqs/ \
       --illumina \
       --prefix "test"
cp .nextflow.log artifacts/symref_bed.nextflow.log

# check that bwa index did not ran:
find work -name .command.sh \
    | xargs cat | grep 'bwa index' \
    && bash -c "echo test failed\: bwa index ran and shouldn\'t have in --ref, --bed mode >> artifacts/test_artifact.log && exit 1"
echo "ran with --bed and --ref: did NOT run bwa index" >> artifacts/test_artifact.log

# clean-up for following tests
rm -rf results && rm -rf work && rm -rf .nextflow*
