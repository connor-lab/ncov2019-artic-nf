// ARTIC ncov workflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {articDownloadScheme} from '../modules/artic.nf' params(params)
include {articGather} from '../modules/artic.nf' params(params)
include {articDemultiplex} from  '../modules/artic.nf' params(params)
include {nanopolishIndex} from  '../modules/artic.nf' params(params)
include {articMinION} from  '../modules/artic.nf' params(params)


// workflow component for artic pipeline
workflow articNcov {
    take:
      ch_runDirectory
    
    main:
      articDownloadScheme()

      articGather(ch_runDirectory)
      
      nanopolishIndex(articGather.out.gathered
                                     .combine(ch_runDirectory))

      if(params.barcode) {
          articDemultiplex(articGather.out.gathered)
          
          articMinION(articGather.out.fastq
                                     .combine(articDemultiplex.out.flatten())
                                     .combine(nanopolishIndex.out.toList())
                                     .combine(articDownloadScheme.out)
                                     .combine(ch_runDirectory))

      } else {
          articMinION(articGather.out.fastq
                                     .combine(articGather.out.fastq)
                                     .combine(nanopolishIndex.out.toList())
                                     .combine(articDownloadScheme.out)
                                     .combine(ch_runDirectory))
      }
}
