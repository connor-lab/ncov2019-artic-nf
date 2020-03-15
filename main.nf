#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import subworkflows
include {articNcovNanopolish} from './workflows/articNcovNanopolish.nf' params(params)
include {articNcovMedaka} from './workflows/articNcovMedaka.nf' params(params)


// main workflow
workflow {
   runDirectory = "${params.directory}"

   Channel.fromPath( "${runDirectory}" )
          .set{ ch_runDirectory }


   main:
     if( params.medaka ) {
         articNcovMedaka(ch_runDirectory)
     } else {
         articNcovNanopolish(ch_runDirectory)
     }
}

