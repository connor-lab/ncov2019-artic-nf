// ARTIC ncov workflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {articDownloadScheme} from '../modules/artic.nf' params(params)
include {articGuppyPlex} from '../modules/artic.nf' params(params)
include {articMinION} from  '../modules/artic.nf' params(params)
include {articRemoveUnmappedReads} from '../modules/artic.nf' params(params)

include {collateSamples} from '../modules/upload.nf' params(params)
include {uploadToCLIMB} from '../modules/upload.nf' params(params)


// workflow component for artic pipeline
workflow sequenceAnalysis {
    take:
      ch_runDirectory
      ch_runFastqDirs
    
    main:
      articDownloadScheme()

      articGuppyPlex(ch_runDirectory.combine(ch_runFastqDirs.flatten()))

      articMinION(articGuppyPlex.out.fastq
                             .combine(articDownloadScheme.out)
                             .combine(ch_runDirectory))

      articRemoveUnmappedReads(articMinION.out.sorted_bam)

    emit:
      bams = articRemoveUnmappedReads.out
      fastas = articMinION.out.consensus_fasta
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

workflow articNcovNanopore {
    take:
      ch_runDirectory
      ch_fastqDirs

    main:
      sequenceAnalysis(ch_runDirectory, ch_fastqDirs)
/*
      if ( params.upload ) {

        Channel.fromPath("${params.CLIMBkey}")
               .set{ ch_CLIMBkey }

        CLIMBrsync(sequenceAnalysis.out.bams, sequenceAnalysis.out.fastas, ch_CLIMBkey )
      }
*/
}

