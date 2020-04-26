#!/bin/bash
set -eo pipefail
export PATH=/opt/conda/bin:$PATH

mkdir conda_cache_dir

# first NF run will create the conda env in the cache dir
echo run pipeline with conda --cache to create cache.. >> artifacts/test_artifact.log
NXF_VER=20.03.0-edge nextflow run main.nf \
       -profile conda \
       --cache $PWD/conda_cache_dir \
       --directory $PWD/.github/data/fastqs/ \
       --illumina \
       --prefix test
cp .nextflow.log ./artifacts/cache_creation.conda.profile.nextflow.log

cat .nextflow.log | grep 'Conda create complete env=/home/runner/work/ncov2019-artic-nf/ncov2019-artic-nf/environment-illumina.yml path=/home/runner/work/ncov2019-artic-nf/ncov2019-artic-nf/conda_cache_dir/artic-ncov2019-illumina-' \
    && echo "Conda env created in cache dir" >> artifacts/test_artifact.log \
	|| bash -c "echo test failed\: Conda environment not created as expected >> artifacts/test_artifact.log && exit 1"
    
rm -rf results && rm -rf work && rm -rf .nextflow*
# second NF run will use the conda env created in the previous run
echo re-run pipeline with conda --cache.. >> artifacts/test_artifact.log
NXF_VER=20.03.0-edge nextflow run main.nf \
       -profile conda \
       --cache $PWD/conda_cache_dir \
       --directory $PWD/.github/data/fastqs/ \
       --illumina \
       --prefix test
cp .nextflow.log ./artifacts/cache_use.conda.profile.nextflow.log

cat .nextflow.log | grep 'Conda found local env for environment=/home/runner/work/ncov2019-artic-nf/ncov2019-artic-nf/environment-illumina.yml; path=/home/runner/work/ncov2019-artic-nf/ncov2019-artic-nf/conda_cache_dir/artic-ncov2019-illumina-' \
    && echo "Conda env found in cache dir" >> artifacts/test_artifact.log \
	|| bash -c "echo test failed\: Conda environment not not found in cache dir >> artifacts/test_artifact.log && exit 1"

# clean-up for following tests
rm -rf results && rm -rf work && rm -rf .nextflow*
