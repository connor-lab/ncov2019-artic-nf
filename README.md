# GMS-artic (ncov2019-artic-nf)

#Test
A nextflow pipeline with a GMS touch for running the ARTIC network's fieldbioinformatics tools (https://github.com/artic-network/fieldbioinformatics).

------------
#### Major changes

* The illumina and nanopore tracks automatically run pangolin and nextclade.
* Generates report for base changes.

###### 1. gms-artic in ngp-gms

*for nanopore analysis (default is "midnight")*
```
    sample_name
         |___ fast5_pass/
         |___ fastq_pass/
         |___ sequencing_summary.txt
```
*for illumina analysis*
```
    sample_name     
         |___ fastq/
```
#### Manual running of analysis pipeline
###### 2. Run Illumina pipeline
```
$ nextflow run main.nf -profile singularity,sge \
    --illumina --prefix "test_illumina"     \
    --directory .github/data/fastqs/    \
    --outdir illumina_test
```

###### 3. Run Nanopore Pipeline
###### **Deafult is "midnight" protocol**
```
$ nextflow run main.nf -profile singularity \
    --nanopolish --prefix "midnight" \
    --basecalled_fastq /home/test/fastq_pass/ \
    --fast5_pass /home/test/fast5_pass/ \
    --sequencing_summary /home/test/sequencing_summary_FAP82331_657703c9.txt \
    --scheme-directory primer_schemes/midnight/nCoV-2019/V1/ \
    --outdir /home/test/midnight_test -with-report midnight
```
###### **To run the artic pipeline, please change the [nanopore.config](https://github.com/JD2112/gms-artic/blob/master/conf/nanopore.config) 'min_length' (default = 400) and 'max_length' (default = 700)**

```
$ nextflow run main.nf -profile singularity,sge \
    --nanopolish --prefix "test_nanopore" \
    --basecalled_fastq .github/data/nanopore/20200311_1427_X1_FAK72834_a3787181/fastq_pass/ \
    --fast5_pass .github/data/nanopore/20200311_1427_X1_FAK72834_a3787181/fast5_pass/ \
    --sequencing_summary .github/data/nanopore/20200311_1427_X1_FAK72834_a3787181/sequencing_summary_FAK72834_298b7829.txt \
    --outdir nanopore_test
```
