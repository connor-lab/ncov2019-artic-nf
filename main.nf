#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import subworkflows
include {articNcov} from './workflows/articNcov.nf' params(params)


// main workflow
workflow {
   runDirectory = "${params.directory}"

   Channel.fromPath( "${runDirectory}" )
          .set{ ch_runDirectory }


   main:
     articNcov(ch_runDirectory)
}

