#### Quick-start

##### Illumina
`nextflow run connor-lab/ncov2019-artic-nf [-profile conda,singularity,docker,slurm,lsf] --illumina --prefix "output_file_prefix" --directory /path/to/reads`

You can also use cram file input by passing the --cram flag.
You can also specify cram file output by passing the --outCram flag.

For production use at large scale, where you will run the workflow many times, you can avoid cloning the scheme repository, creating an ivar bed file and indexing the reference every time by supplying both --ivarBed /path/to/ivar-compatible.bed and --alignerRefPrefix /path/to/bwa-indexed/ref.fa.

Alternatively you can avoid just the cloning of the scheme repository to remain on a fixed revision of it over time by passing --schemeRepoURL /path/to/own/clone/of/github.com/artic-network/artic-ncov2019. This removes any internet access from the workflow except for the optional upload steps.

##### Nanopore
###### Nanopolish
`nextflow run connor-lab/ncov2019-artic-nf [-profile conda,singularity,docker,slurm,lsf] --nanopolish --prefix "output_file_prefix" --basecalled_fastq /path/to/directory --fast5_pass /path/to/directory --sequencing_summary /path/to/sequencing_summary.txt`

###### Medaka
 `nextflow run connor-lab/ncov2019-artic-nf [-profile conda,singularity,docker,slurm,lsf] --medaka --prefix "output_file_prefix" --basecalled_fastq /path/to/directory --fast5_pass /path/to/directory --sequencing_summary /path/to/sequencing_summary.txt`

#### Installation
An up-to-date version of Nextflow is required because the pipeline is written in DSL2. Following the instructions at https://www.nextflow.io/ to download and install Nextflow should get you a recent-enough version. 

#### Containers
This repo contains both [Singularity]("https://sylabs.io/guides/3.0/user-guide/index.html") and Dockerfiles. You can build the Singularity containers locally by running `scripts/build_singularity_containers.sh` and use them with `-profile singularity` The containers will be available from Docker/Singularityhub shortly.

#### Conda
The repo contains a environment.yml files which automatically build the correct conda env if `-profile conda` is specifed in the command. Although you'll need `conda` installed, this is probably the easiest way to run this pipeline.

--cache /some/dir can be specified to have a fixed, shared location to store the conda build for use by multiple runs of the workflow.

#### Executors
By default, the pipeline just runs on the local machine. You can specify `-profile slurm` to use a SLURM cluster, or `-profile lsf` to use an LSF cluster. In either case you may need to also use one of the COG-UK institutional config profiles (phw or sanger), or provide queue names to use in your own config file.

#### Profiles
You can use multiple profiles at once, separating them with a comma. This is described in the Nextflow [documentation](https://www.nextflow.io/docs/latest/config.html#config-profiles) 

#### Config
Common configuration options are set in `conf/base.config`. Workflow specific configuration options are set in `conf/nanopore.config` and `conf/illumina.config` They are described and set to sensible defaults (as suggested in the [nCoV-2019 novel coronavirus bioinformatics protocol](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html "nCoV-2019 novel coronavirus bioinformatics protocol"))

##### Options
- `--outdir` sets the output directory.
- `--bwa` to swap to bwa for mapping (nanopore only).

##### Workflows

###### Nanopore
Use `--nanopolish` or `--medaka` to run these workflows. `--basecalled_fastq` should point to a directory created by `guppy_basecaller` (if you ran with no barcodes), or `guppy_barcoder` (if you ran with barcodes). It is imperative that the following `guppy_barcoder` command be used for demultiplexing:

```
guppy_barcoder --require_barcodes_both_ends -i run_name -s output_directory --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"
```

###### Illumina
The Illumina workflow leans heavily on the excellent [ivar](https://github.com/andersen-lab/ivar) for primer trimming and consensus making. This workflow will be updated to follow ivar, as its also in very active development! Use `--illumina` to run the Illumina workflow. Use `--directory` to point to an Illumina output directory usually coded something like: `<date>_<machine_id>_<run_no>_<some_zeros>_<flowcell>`. The workflow will recursively grab all fastq files under this directory, so be sure that what you want is in there, and what you don't, isn't! 

Important config options are:

| Option | Description |
|:-------|------------:|
|allowNoprimer | Allow reads that don't have primer sequence? Ligation prep = false, nextera = true|
|illuminaKeepLen | Length of illumina reads to keep after primer trimming|
|illuminaQualThreshold | Sliding window quality threshold for keeping reads after primer trimming (illumina)|
|mpileupDepth | Mpileup depth for ivar|
|ivarFreqThreshold | ivar frequency threshold for variant|
|ivarMinDepth | Minimum coverage depth to call variant|

#### QC
A script to do some basic COG-UK QC is provided in `bin/qc.py`. This currently tests if >50% of reference bases are covered by >10 reads (Illumina) or >20 reads (Nanopore), OR if there is a stretch of more than 10 Kb of sequence without N - setting qc_pass in `<outdir>/<prefix>.qc.csv` to TRUE. `bin/qc.py` can be extended to incorporate any QC test, as long as the script outputs a csv file a "qc_pass" last column, with samples TRUE or FALSE.

#### Output
A subdirectory for each process in the workflow is created in `--outdir`. A `qc_pass_climb_upload` subdirectory containing files important for [COG-UK](https://github.com/COG-UK) is created. 
