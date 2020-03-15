# ncov2019-artic-nf
A Nextflow pipeline for running the ARTIC network's fieldbioinformatics tools (https://github.com/artic-network/fieldbioinformatics), with a focus on ncov2019 

#### Introduction

------------

This Nextflow pipeline automates the ARTIC network [nCoV-2019 novel coronavirus bioinformatics protocol](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html "nCoV-2019 novel coronavirus bioinformatics protocol"). It should handle both barcoded and non-barcoded nanopore run data, although the barcoded workflow isn't yet tested.

##### Quick-start
`nextflow run connor-lab/ncov2019-artic-nf [-c /path/to/additional_config.nf]  [-profile conda,singularity,docker,slurm] [--minimap] [--barcode] --directory /path/to/reads`

##### Installation
An up-to-date version of Nextflow is required because the pipeline is written in DSL2. Following the instructions at https://www.nextflow.io/ to download and install Nextflow should get you a recent-enough version. 

##### Containers
This repo contains both [Singularity]("https://sylabs.io/guides/3.0/user-guide/index.html") and Dockerfiles. You can build the Singularity containers locally by running `scripts/build_singularity_containers.sh` and use them with `-profile singularity` The containers will be available from Docker/Singularityhub shortly.

##### Conda
The repo contains an environment.yml file which automatically builds the correct conda env if `-profile conda` is specifed in the command. Although you'll need `conda` installed, this is probably the easiest way to run this pipeline.

##### Executors
By default, the pipeline just runs on the local machine. You can specify `-profile slurm` to use a SLURM cluster. Additional [config]("https://www.nextflow.io/docs/latest/executor.html#slurm") can be added with `-c /path/to/config.nf`. 

##### Profiles
You can use multiple profiles at once, separating them with a comma. This is described in the Nextflow [documentation](https://www.nextflow.io/docs/latest/config.html#config-profiles) 

##### Config
Configuration options are set in `conf/base.config`. They are described and set to sensible defaults (as suggested in the [nCoV-2019 novel coronavirus bioinformatics protocol](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html "nCoV-2019 novel coronavirus bioinformatics protocol"))

##### Options
The only required option is `--directory`, which should point to a nanopore output directory containing:
- drift_correction_ACG216_1e627889.csv
- duty_time_ACG216_1e627889.csv
- fast5_fail
- fast5_pass
- fastq_fail
- fastq_pass
- final_summary_ACG216_1e627889.txt
- report_ACG216_20200306_1325_cfb936bc.md
- report_ACG216_20200306_1325_cfb936bc.pdf
- sequencing_summary_ACG216_1e627889.txt
- throughput_ACG216_1e627889.csv

Use `--minimap` to swap to minimap for mapping, `--barcode` if you ran with barcodes, and `--medaka` to run the experimental medaka version of the pipeline.
