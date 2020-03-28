#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import subworkflows
include {articNcovNanopore} from './workflows/articNcovNanopore.nf' params(params)
include {ncovIllumina} from './workflows/illuminaNcov.nf' params(params)

// main workflow
workflow {
  
   runDirectory = "${params.directory}"
   if ( params.illumina ) {
  
       Channel.fromFilePairs( params.fastqSearchPath, flat: true)
              .set{ ch_filePairs }
   }
   else {
       // Need the base run directory
       Channel.fromPath( "${runDirectory}" )
              .set{ ch_runDirectory }

       
       // Check to see if we have barcodes
       def nanoporeBarcodeDirs = new FileNameByRegexFinder().getFileNames(params.directory, /.*\/fastq_pass\/.*barcode\d\d$/)
       
       
       if( nanoporeBarcodeDirs ) {
            // Yes, barcodes!
            Channel.fromPath( "${runDirectory}/fastq_pass/barcode*", type: 'dir', maxDepth: 1 )
                   .filter{ it.listFiles().size() > 10 }
                   .map{ it.getBaseName().toString() }
                   .set{ ch_fastqDirs }
       } else {
            // No, no barcodes
            Channel.fromPath( "${runDirectory}/fastq_pass/", type: 'dir', maxDepth: 1 )
                    .set{ ch_fastqDirs }
      }
   }

   main:
     if ( params.nanopolish || params.medaka ) {
         articNcovNanopore(ch_runDirectory, ch_fastqDirs)
     } else if ( params.illumina ) {
         ncovIllumina(ch_filePairs)
     } else {
         println("Please select a workflow with --nanopolish, --illumina or --medaka")
     }
     
}

