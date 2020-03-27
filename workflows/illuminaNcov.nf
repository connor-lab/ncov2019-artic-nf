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

include {collateSamples} from '../modules/upload.nf' params(params)
include {uploadToCLIMB} from '../modules/upload.nf' params(params)

include {cramToFastq} from '../modules/illumina.nf' params(params)

workflow sequenceAnalysis {
    take:
      ch_filePairs

    main:
      articDownloadScheme()

      makeIvarBedfile(articDownloadScheme.out)

      readTrimming(ch_filePairs)

      readMapping(articDownloadScheme.out.combine(readTrimming.out))

      trimPrimerSequences(makeIvarBedfile.out.combine(readMapping.out))

      makeConsensus(trimPrimerSequences.out)

    emit:
      bams = trimPrimerSequences.out
      fastas = makeConsensus.out
      
}

workflow CLIMBrsync {
    take:
      ch_sequenceAnalysisBAMs
      ch_sequenceAnalysisFastas
      ch_CLIMBkey

    main:
      collateSamples(ch_sequenceAnalysisBAMs.join(ch_sequenceAnalysisFastas, by: 0))
      uploadToCLIMB(ch_CLIMBkey.combine(collateSamples.out.collect().toList()))
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

workflow ncovIlluminaCram {
    take:
      ch_cramDirectory
    main:
      cramToFastq(ch_cramDirectory)
      ncovIllumina(cramToFastq.out)
}

