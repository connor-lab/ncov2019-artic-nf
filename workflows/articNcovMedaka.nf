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

// workflow component for artic pipeline
workflow articNcovMedaka {
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
}
