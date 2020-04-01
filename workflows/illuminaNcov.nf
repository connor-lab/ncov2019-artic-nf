#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {articDownloadScheme } from '../modules/artic.nf' 
include {makeIvarBedfile} from '../modules/illumina.nf' 
include {readTrimming} from '../modules/illumina.nf' 
include {readMapping} from '../modules/illumina.nf' 
include {trimPrimerSequences} from '../modules/illumina.nf' 
include {makeConsensus} from '../modules/illumina.nf' 

// import subworkflows
include {CLIMBrsync} from './upload.nf'


workflow sequenceAnalysis {
    take:
      ch_filePairs

    main:
      articDownloadScheme()

      makeIvarBedfile(articDownloadScheme.out)

      readTrimming(ch_filePairs)

      readMapping(articDownloadScheme.out.combine(readTrimming.out))

      trimPrimerSequences(makeIvarBedfile.out.combine(readMapping.out))

      makeConsensus(trimPrimerSequences.out.ptrim)

    emit:
      bams = trimPrimerSequences.out.mapped
      fastas = makeConsensus.out
      
}

workflow ncovIllumina {
    take:
      ch_filePairs

    main:
      sequenceAnalysis(ch_filePairs)
      
      if ( params.upload ) {
        
        Channel.fromPath("${params.CLIMBkey}")
               .set{ ch_CLIMBkey }
      
        CLIMBrsync(sequenceAnalysis.out.bams, sequenceAnalysis.out.fastas, ch_CLIMBkey )
      }
}

