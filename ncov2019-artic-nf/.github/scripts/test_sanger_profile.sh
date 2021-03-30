#!/bin/bash
set -eo pipefail
export PATH=/opt/conda/bin:$PATH

# test --sanger profile
# there are only 2 available cpus in the github runner execution
sed -i s'/cpus = 4/cpus = 2/'g conf/coguk/sanger.config
echo run pipeline in --illumina mode with --sanger profile.. >> artifacts/test_artifact.log
NXF_VER=20.03.0-edge nextflow run ./main.nf \
       -profile sanger,singularity \
       --directory $PWD/.github/data/fastqs/ \
       --illumina \
       --prefix test
cp .nextflow.log ./artifacts/sanger.profile.nextflow.log
cp -r results results_sanger_profile

# check that sanger profile activated 4 cpus on bwa mem:
find work -name .command.err \
    | xargs cat | grep '\[main\] CMD: bwa mem -t 2' \
    && echo "sanger profile: bwa started with 4 cpus" >> artifacts/test_artifact.log \
	|| bash -c "echo test failed\: bwa mem should have used 4 cpus with --sanger profile >> artifacts/test_artifact.log && exit 1"

# check that sanger profile activated params.allowNoprimer = false:
find work -name .command.sh \
    | xargs cat | grep 'ivar trim -e' \
	    && bash -c "echo test failed\: should not have run \'ivar trim -e\' with --sanger profile >> artifacts/test_artifact.log && exit 1"
echo "sanger profile: did NOT ran ivar trim -e" >> artifacts/test_artifact.log 

# clean-up for following tests
rm -rf results && rm -rf work && rm -rf .nextflow*
