#!/bin/bash

# ont_artic_test_a03b0d59

mkdir -p /work/runs/ont_artic_test_a03b0d59
cd /work/runs/ont_artic_test_a03b0d59
nextflow run \
        /data/pipelines/ncov2019-artic-nf/main.nf \
        -with-trace -with-report -with-timeline -with-dag dag.png \
        --prefix nanopore \
        -profile singularity \
        -process.executor slurm \
        --objstore false \
        --catsup /data/inputs/s3/oracle-test/a03b0d59-0908-4597-8600-25f9e787c762 \
        --bucket catsup-test \
        --varCaller medaka \
        --refmap '"{}"' \
        --pipeline_name oxforduni-ncov2019-artic-nf-nanopore \
        --run_uuid ffdd1e7f-2aaa-43a7-a230-f6b991bf4631 \
        --head_node_ip 10.0.1.2 \
        --outdir /work/output/ont_artic_test_a03b0d59 \
        > nextflow.txt

# ont_viridian_test_a03b0d59

mkdir -p /work/runs/ont_viridian_test_a03b0d59
cd /work/runs/ont_viridian_test_a03b0d59
nextflow  run \
        /data/pipelines/ncov2019-artic-nf/main.nf \
        -with-trace -with-report -with-timeline -with-dag dag.png \
        --prefix nanopore \
        -profile singularity \
        -process.executor slurm \
        --objstore false \
        --catsup /data/inputs/s3/oracle-test/a03b0d59-0908-4597-8600-25f9e787c762 \
        --bucket catsup-test \
        --varCaller viridian \
        --refmap '"{}"' \
        --pipeline_name oxforduni-ncov2019-artic-nf-nanopore \
        --run_uuid b6a04e93-e031-4a80-9ece-0a279f9b1fe4 \
        --head_node_ip 10.0.1.2 \
        --outdir /work/output/ont_viridian_test_a03b0d59 \
        > nextflow.txt

# illumina_artic_test_bc19173

mkdir -p /work/runs/illumina_artic_test_bc19173
cd /work/runs/illumina_artic_test_bc19173

nextflow run /data/pipelines/ncov2019-artic-nf/main.nf \
        -with-trace -with-report -with-timeline -with-dag dag.png \
        --readpat '*{1,2}.fastq.gz' \
        --illumina --prefix illumina \
        -profile singularity \
        -process.executor slurm \
        --objstore false \
        --catsup /data/inputs/s3/oracle-test/bc19173b-2efe-401e-8c75-2b393af77b8e \
        --bucket catsup-test \
        --varCaller iVar \
        --refmap '"{}"' \
        --pipeline_name oxforduni-ncov2019-artic-nf-illumina \
        --run_uuid 19f03473-156a-4cec-a947-f7cfd1a03947 \
        --head_node_ip 10.0.1.2 \
        --outdir /work/output/illumina_artic_test_bc19173 \
        > nextflow.txt

# illumina_Viridian_test_bc19173

mkdir -p /work/runs/illumina_viridian_test_bc19173
cd /work/runs/illumina_viridian_test_bc19173

nextflow run /data/pipelines/ncov2019-artic-nf/main.nf \
        -with-trace -with-report -with-timeline -with-dag dag.png \
        --readpat '*{1,2}.fastq.gz' \
        --illumina --prefix illumina \
        -profile singularity \
        -process.executor slurm \
        --objstore false \
        --catsup /data/inputs/s3/oracle-test/bc19173b-2efe-401e-8c75-2b393af77b8e \
        --bucket catsup-test \
        --varCaller viridian \
        --refmap '"{}"' \
        --pipeline_name oxforduni-ncov2019-artic-nf-illumina \
        --run_uuid 387691ae-1f78-444d-a317-23443472b188 \
        --head_node_ip 10.0.1.2 \
        --outdir /work/output/illumina_viridian_test_bc19173 \
        > nextflow.txt

