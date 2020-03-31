// ARTIC ncov workflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {articDownloadScheme} from '../modules/artic.nf' 
include {articGuppyPlex} from '../modules/artic.nf' 
include {articMinION} from  '../modules/artic.nf' 
include {articRemoveUnmappedReads} from '../modules/artic.nf' 

include {collateSamples} from '../modules/upload.nf' 
include {uploadToCLIMB} from '../modules/upload.nf' 


// workflow component for artic pipeline
workflow sequenceAnalysis {
    take:
      ch_runFastqDirs
      ch_fast5Pass
      ch_seqSummary
    
    main:
      articDownloadScheme()
      
      articGuppyPlex(ch_runFastqDirs.flatten())

      articMinION(articGuppyPlex.out.fastq
                             .combine(articDownloadScheme.out)
                             .combine(ch_fast5Pass)
                             .combine(ch_seqSummary))

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
      ch_seqSummary

    main:
      sequenceAnalysis(ch_runDirectory, ch_fastqDirs, ch_seqSummary)
/*
      if ( params.upload ) {

        Channel.fromPath("${params.CLIMBkey}")
               .set{ ch_CLIMBkey }

        CLIMBrsync(sequenceAnalysis.out.bams, sequenceAnalysis.out.fastas, ch_CLIMBkey )
      }
*/
}

