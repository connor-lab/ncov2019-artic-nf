// ARTIC ncov workflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {articDownloadScheme} from '../modules/artic.nf' params(params)
include {articGather} from '../modules/artic.nf' params(params)
include {articDemultiplex} from  '../modules/artic.nf' params(params)
include {nanopolishIndex} from  '../modules/artic.nf' params(params)
include {articMinIONMedaka} from  '../modules/artic.nf' params(params)
include {articRemoveUnmappedReads} from '../modules/artic.nf' params(params)

include {collateSamples} from '../modules/upload.nf' params(params)
include {uploadToCLIMB} from '../modules/upload.nf' params(params)


// workflow component for artic pipeline
workflow sequenceAnalysis {
    take:
      ch_runDirectory
    
    main:
      articDownloadScheme()

      articGather(ch_runDirectory)
      
      if(params.barcode) {
          articDemultiplex(articGather.out.gathered)
          
          articMinIONMedaka(articDemultiplex.out.flatten()
                                          .combine(articDownloadScheme.out)
                                          .combine(ch_runDirectory))
          
          articRemoveUnmappedReads(articMinIONMedaka.out.sorted_bam)

      } else {
          articMinIONMedaka(articGather.out.fastq
                                     .combine(articDownloadScheme.out)
                                     .combine(ch_runDirectory))
          
          articRemoveUnmappedReads(articMinIONMedaka.out.sorted_bam)
      }

    emit:
      bams = articRemoveUnmappedReads.out
      fastas = articMinIONMedaka.out.consensus_fasta
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

workflow articNcovMedaka {
    take:
      ch_runDirectory

    main:
      sequenceAnalysis(ch_runDirectory)

      if ( params.upload ) {

        Channel.fromPath("${params.CLIMBkey}")
               .set{ ch_CLIMBkey }

        CLIMBrsync(sequenceAnalysis.out.bams, sequenceAnalysis.out.fastas, ch_CLIMBkey )
      }
}

