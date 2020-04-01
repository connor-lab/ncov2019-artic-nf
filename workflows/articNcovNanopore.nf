// ARTIC ncov workflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {articDownloadScheme} from '../modules/artic.nf' 
include {articGuppyPlex} from '../modules/artic.nf' 
include {articMinION} from  '../modules/artic.nf' 
include {articRemoveUnmappedReads} from '../modules/artic.nf' 

include {makeQCCSV} from '../modules/qc.nf'
include {writeQCSummaryCSV} from '../modules/qc.nf'

include {collateSamples} from '../modules/upload.nf'


// import subworkflows
include {CLIMBrsync} from './upload.nf'


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
                             .combine(articDownloadScheme.out.scheme)
                             .combine(ch_fast5Pass)
                             .combine(ch_seqSummary))

      articRemoveUnmappedReads(articMinION.out.mapped)

      makeQCCSV(articMinION.out.ptrim.join(articMinION.out.consensus_fasta, by: 0)
                           .combine(articDownloadScheme.out.reffasta))

      makeQCCSV.out.splitCsv()
                   .unique()
                   .branch {
                       header: it[-1] == 'qc_pass'
                       fail: it[-1] == 'FALSE'
                       pass: it[-1] == 'TRUE'
                   }
                   .set { qc }

     writeQCSummaryCSV(qc.header.concat(qc.pass).concat(qc.fail).toList())

     collateSamples(qc.pass.map{ it[0] }
                           .join(articMinION.out.consensus_fasta, by: 0)
                           .join(articRemoveUnmappedReads.out))



    emit:
      qc_pass = collateSamples.out

}
     

workflow articNcovNanopore {
    take:
      ch_runDirectory
      ch_fastqDirs
      ch_seqSummary

    main:
      sequenceAnalysis(ch_runDirectory, ch_fastqDirs, ch_seqSummary)
      if ( params.upload ) {

        Channel.fromPath("${params.CLIMBkey}")
               .set{ ch_CLIMBkey }

        CLIMBrsync(sequenceAnalysis.out.qc_pass, ch_CLIMBkey )
      }
}

