
![logo](workflow-image/logo.png)

A nextflow pipeline with a GMS touch for running the ARTIC network's fieldbioinformatics tools (https://github.com/artic-network/fieldbioinformatics).

### Table of contents - 
- [Version updates](#Version-updates)
- [How to run in NGP server](#How-to-run-in-NGP-server)
  - [Datafile structure](#Datafile-structure)
  - [Pipeline run command](#Manual-running-of-analysis-pipeline)
    - [Illumina pipeline](#Run-Illumina-pipeline)
    - [Nanopore pipeline](#Run-Nanopore-Pipeline)
  - [Primer scheme parameters setup](#Primer-scheme-parameters-setup)
- [Useful information](#Useful-information) 
------------
## Version updates
### v2.0.0
#### Major updates
- Docker container separated for Pangolin typing -
    - Illumina container: [gms-artic-illumina](https://hub.docker.com/repository/docker/genomicmedicinesweden/gms-artic-illumina)
    - Nanopore container: [gms-artic-nanopore](https://hub.docker.com/repository/docker/genomicmedicinesweden/gms-artic-nanopore)
    - Pangolin container: [gms-artic-pangolin](https://hub.docker.com/repository/docker/genomicmedicinesweden/gms-artic-pangolin)
    - pycoQC container  : pycoqc
- Added separate package version files for each workflow -
    - versions      : for Illumina and Nanopore
    - pangoversion  : for pangolin typing
- Illumina analysis additional features -
    - flagstat
    - depth
    - VEP annotation
- Illumina results works for sc2reporter visualization
- Nanopore analysis additional features (artic & medaka)-
    - [fastqc](https://github.com/s-andrews/FastQC)
    - [multiqc](https://multiqc.info)
    - [pycoQC](https://github.com/a-slide/pycoQC)

### v1.8.0
#### Minor updates

- Pangolin v4 support
- Updated Picard arguments
- FastQC commands can be added from config
- Added version of pangolin to build_dockerfile

#### Bug fixes
- Fixed build_dockerfile
- Fixed R issue
- Fixed mamba issue

#### Major changes

* The illumina and nanopore tracks automatically run pangolin and nextclade.
* Generates report for base changes.

## How to run in NGP server

### Datafile structure
#### gms-artic in ngp-gms

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
### Manual running of analysis pipeline
#### Run Illumina pipeline
```
$ nextflow run main.nf -profile singularity,sge \
    --illumina --prefix "test_illumina"     \
    --directory .github/data/fastqs/    \
    --outdir illumina_test
```

#### Run Nanopore Pipeline
###### **Deafult is "midnight" protocol**
```
$ nextflow run main.nf -profile singularity \
    --nanopolish --prefix "midnight" \
    --basecalled_fastq /home/test/fastq_pass/ \
    --fast5_pass /home/test/fast5_pass/ \
    --sequencing_summary /home/test/sequencing_summary_FAP82331_657703c9.txt \
    --scheme-directory midnight-primer/V1/ \
    --outdir /home/test/midnight_test -with-report midnight
```

### Primer scheme parameters setup
###### --scheme: To use the primer list, add --scheme to the CLI, eg., use 'nCoV-2019/V3' for artic primers or 'midnight-primer/V1'

```
--scheme nCoV-2019/V3/
--scheme midnight-primers/V1/
--scheme eden-primers/V1/

```

## Useful information
###### **To run the artic pipeline, please change the [nanopore.config](https://github.com/JD2112/gms-artic/blob/master/conf/nanopore.config) 'min_length' (default = 400) and 'max_length' (default = 700)**

```
$ nextflow run main.nf -profile singularity,sge \
    --nanopolish --prefix "test_nanopore" \
    --basecalled_fastq .github/data/nanopore/20200311_1427_X1_FAK72834_a3787181/fastq_pass/ \
    --fast5_pass .github/data/nanopore/20200311_1427_X1_FAK72834_a3787181/fast5_pass/ \
    --sequencing_summary .github/data/nanopore/20200311_1427_X1_FAK72834_a3787181/sequencing_summary_FAK72834_298b7829.txt \
    --outdir nanopore_test
```
#### To update your container image to the latest version from [dockerhub](https://hub.docker.com/orgs/genomicmedicinesweden/repositories), please delete your local image first before running the analysis pipeline. 
