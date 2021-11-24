# GMS-artic (ncov2019-artic-nf)

A nextflow pipeline with a GMS touch for running the ARTIC network's fieldbioinformatics tools (https://github.com/artic-network/fieldbioinformatics).

______________
#### Release tags
v1.5.0
Date: 2021-11-24
* check and update pangolin typing everyday automatically
* build docker containers with updated pangolin version
* Nanopore Medaka workflow contains pangolin version 

v1.4.0
Date: 2021-09-13
* Nanopore Midnight protocol set as the default for nanopore analysis
* pangolinTyping updated to latest release

v1.3.0
Date: 2020-12-23
* Add script to type variants and consequences
* Update deps for typing script
* Rename type_variants.py so its not confused with AR group tool
* Add NF processes for typing and making summary csvs
* Add example yaml and gff files for typing
* Add a typing workflow that can be attached to nanopore or illumina data
* Update config options for typing
* Update conda deps for typing script
* Add typing workflow to illumina pipeline
* Add named outputs to typing process and update workflow
* Add typing workflow to nanopore pipeline
* Update help with new options
* Tidy up commented line in qc.nf
* Make output of typing csvs optional for when no variants pass filters
* Output boxplots from nanopore pipeline
* Add rudimentary tests of typing code
* Add N501Y.V2 variants to types.yaml because why not
* Add comments to example typing.yaml file

v1.2.0
Date: 2020-11-18
* Update ivar version to 1.3
* Pin artic version to 1.1.3
* Test against v1.1.1
* Remove depth plot png files from tests

v1.1.1
Date: 2020-09-26
* Fix samtools sort reversion
* Move BAM header cleaning to trimPrimerSequences process

v1.1.0
Date: 2020-09-24
* Move fastq globbing into a module and join paths properly
* Use makeFastqSearchPath to make fileglobpath in main.nf
* Add PHW specific config for sample names
* Restore test against previous release
* Sometimes paths have a trailing slash or two, remove

v1.0.0
Date: 2020-09-15
* Add new option to remove samples with low reads after guppyplex 
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
