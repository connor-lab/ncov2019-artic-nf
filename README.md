# GMS-artic (ncov2019-artic-nf)

A nextflow pipeline with a GMS touch for running the ARTIC network's fieldbioinformatics tools (https://github.com/artic-network/fieldbioinformatics).

------------
#### Major changes

* The illumina and nanopore tracks automatically run pangolin

##### gms-artic in ngp-gms

###### 2. Run Illumina pipeline
```
$ nextflow run main.nf -profile singularity,sge \
    --illumina --prefix "test_illumina"     \
    --directory .github/data/fastqs/    \
    --outdir illumina_test
```

###### 3. Run Nanopore Pipeline
####### deafult is "midnight" protocol
```
$ nextflow run main.nf -profile singularity \
    --nanopolish --prefix "midnight" \
    --basecalled_fastq /home/test/fastq_pass/ \
    --fast5_pass /home/test/fast5_pass/ \
    --sequencing_summary /home/test/sequencing_summary_FAP82331_657703c9.txt \
    --scheme-directory primer_schemes/midnight/nCoV-2019/V1/ \
    --outdir /home/test/midnight_test -with-report midnight
```
####### To run the artic pipeline

```
$ 

