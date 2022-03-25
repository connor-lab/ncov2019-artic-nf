#!/bin/bash

sudo mkdir -p /work/tmp
sudo chown ubuntu:ubuntu /work/tmp
echo \
"SARS-CoV-2_reference_ox,mmm-artic-ill-s11511-1
SARS-CoV-2_reference_ox,mmm-artic-ill-s12220-1
SARS-CoV-2_reference_ox,mmm-artic-ill-s12368-1
SARS-CoV-2_reference_ox,mmm-artic-ill-s16621-1" \
    > /work/tmp/illumina_data_short.csv

echo \
"SARS-CoV-2_reference_ox,mmm-artic-ont-s11511-1
SARS-CoV-2_reference_ox,mmm-artic-ont-s12220-4
SARS-CoV-2_reference_ox,mmm-artic-ont-s12368-1
SARS-CoV-2_reference_ox,mmm-artic-ont-s16621-3" \
    > /work/tmp/ONT_data_short.csv

pushd /data/pipelines/ncov2019-artic-nf
git_version=$(git describe --tags)
popd

# ont_viridian_test
test_name=ont_viridian_short
echo Running ${test_name} test workflow
mkdir -p /work/runs/${test_name}_test
cd /work/runs/${test_name}_test

nextflow kuberun \
        oxfordmmm/ncov2019-artic-nf \
        -with-trace trace.txt -with-report report.html -with-dag dag.png \
        -r ${git_version} \
        --prefix nanopore \
        -profile oke \
        --objstore /work/tmp/ONT_data_short.csv \
        --varCaller viridian \
        --refmap '"{}"' \
        --run_uuid ${test_name}_test \
        --TESToutputMODE true \
        --outdir /work/output/${test_name}_test \
        > nextflow.txt

sudo chown ubuntu:ubuntu /work/output/${test_name}_test

python3 /data/pipelines/ncov2019-artic-nf/tests/GPAS_tests_summary.py \
        -w /work/runs/${test_name}_test \
        -i /work/output/${test_name}_test/ \
        -t /work/output/${test_name}_test/${test_name}_summary.tsv  \
        -e /data/pipelines/ncov2019-artic-nf/tests/${test_name}_expected.tsv \
        -c /work/output/${test_name}_test/${test_name}_comparison.tsv

# illumina_Viridian_test
test_name=illumina_viridian_short
echo Running ${test_name} test workflow
mkdir -p /work/runs/${test_name}_test
cd /work/runs/${test_name}_test


nextflow kuberun \
        oxfordmmm/ncov2019-artic-nf \
        -with-trace trace.txt -with-report report.html -with-dag dag.png \
        -r ${git_version} \
        --readpat '*{1,2}.fastq.gz' \
        --illumina --prefix illumina \
        -profile oke \
        --objstore /work/tmp/illumina_data_short.csv \
        --varCaller viridian \
        --refmap '"{}"' \
        --run_uuid ${test_name}_test \
        --TESToutputMODE true \
        --outdir /work/output/${test_name}_test \
        > nextflow.txt

sudo chown ubuntu:ubuntu /work/output/${test_name}_test

python3 /data/pipelines/ncov2019-artic-nf/tests/GPAS_tests_summary.py \
       -w /work/runs/${test_name}_test \
       -i /work/output/${test_name}_test/ \
       -t /work/output/${test_name}_test/${test_name}_summary.tsv  \
       -e /data/pipelines/ncov2019-artic-nf/tests/${test_name}_expected.tsv \
       -c /work/output/${test_name}_test/${test_name}_comparison.tsv


