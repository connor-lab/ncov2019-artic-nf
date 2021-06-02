[toc]
# gms-artic in ngp-gms

## 1. login to the system
ssh xxxxxx@xx.x.gu.se \
password \
$

==NOTE: sge process is not supported on the /gms-storage-aws/, so need to run the pipeline locally from /home==

## 2. Run Illumina pipeline
```
$ [your home]
$ cp /gms-storage-aws/gms-artic/ .
$ nextflow run main.nf -profile singularity,sge \
    --illumina --prefix "test_illumina"     \
    --directory .github/data/fastqs/    \
    --outdir illumina_test
```

## 3. Run Nanopore Pipeline
```
$ [your home]
$ cp /gms-storage-aws/gms-artic/ .
$ nextflow run main.nf -profile singularity,sge \
    --nanopolish --prefix "test_nanopore" \
    --basecalled_fastq .github/data/nanopore/20200311_1427_X1_FAK72834_a3787181/fastq_pass/ \
    --fast5_pass .github/data/nanopore/20200311_1427_X1_FAK72834_a3787181/fast5_pass/ \
    --sequencing_summary .github/data/nanopore/20200311_1427_X1_FAK72834_a3787181/sequencing_summary_FAK72834_298b7829.txt \
     --outdir nanopore_test
```
