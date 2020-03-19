#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import subworkflows
include {articNcovNanopolish} from './workflows/articNcovNanopolish.nf' params(params)
include {articNcovMedaka} from './workflows/articNcovMedaka.nf' params(params)
include {ncovIllumina} from './workflows/illuminaNcov.nf' params(params)

// main workflow
workflow {
  
   runDirectory = "${params.directory}"
    
   if ( params.illumina ) {
  
       illuminaSuffixes = ['*_R{1,2}_001', '*_R{1,2}', '*_{1,2}' ]
       fastq_exts = ['.fastq.gz', '.fq.gz']
  
       Channel.fromFilePairs( params.fastqSearchPath, flat: true)
              .set{ ch_filePairs }
   }
   else {
       Channel.fromPath( "${runDirectory}" )
              .set{ ch_runDirectory }
   }

   main:
     if( params.medaka ) {
         articNcovMedaka(ch_runDirectory)
     } else if ( params.nanopolish ) {
         articNcovNanopolish(ch_runDirectory)
     } else if ( params.illumina ) {
         ncovIllumina(ch_filePairs)
     } else {
         println("Please select a workflow with --nanopolish, --illumina or --medaka")
     }
}

