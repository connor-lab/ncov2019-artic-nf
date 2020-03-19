#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {articDownloadScheme } from '../modules/artic.nf' params(params)
include {readsTrimming} from '../modules/illumina.nf' params(params)

workflow ncovIllumina{
    take:
      ch_filePairs

    main:
      articDownloadScheme()

      readsTrimming(ch_filePairs)

}

