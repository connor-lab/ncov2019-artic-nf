#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {articDownloadScheme } from '../modules/artic.nf' params(params)
include {makeIvarBedfile} from '../modules/illumina.nf' params(params)
include {readTrimming} from '../modules/illumina.nf' params(params)
include {readMapping} from '../modules/illumina.nf' params(params)
include {trimPrimerSequences} from '../modules/illumina.nf' params(params)
include {makeConsensus} from '../modules/illumina.nf' params(params)

workflow ncovIllumina{
    take:
      ch_filePairs

    main:
      articDownloadScheme()

      makeIvarBedfile(articDownloadScheme.out)

      readTrimming(ch_filePairs)

      readMapping(articDownloadScheme.out.combine(readTrimming.out))

      trimPrimerSequences(makeIvarBedfile.out.combine(readMapping.out))

      makeConsensus(trimPrimerSequences.out)
}
