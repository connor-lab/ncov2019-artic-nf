#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import subworkflows
include {articNcovNanopore} from './workflows/articNcovNanopore.nf' 
include {ncovIllumina} from './workflows/illuminaNcov.nf'
include {ncovIlluminaCram} from './workflows/illuminaNcov.nf'


if ( params.illumina ) {
   if (! params.directory ) {
       println("Please supply a directory containing fastqs with --directory")
       System.exit(1)
   }
   if ( (params.ivarBed && ! params.alignerRefPrefix) || (!params.ivarBed && params.alignerRefPrefix) ) {
       println("ivarBed and alignerRefPrefix must be supplied together")
       System.exit(1)
   }
} else if ( params.medaka || params.nanopolish ) {
   if (! params.basecalled_fastq ) {
       println("Please supply a directory containing basecalled fastqs with --basecalled_fastq (this is the output directory from guppy_barcoder or guppy_basecaller)")
   }
   if (! params.fast5_pass ) {
       println("Please supply a directory containing fast5 files with --fast5_pass (this is the fast5_pass directory)")
   }
   if (! params.sequencing_summary ) {
       println("Please supply the path to the sequencing_summary.txt file from your run with --sequencing_summary")
       System.exit(1)
   }
   if ( params.ivarBed || params.alignerRefPrefix ) {
       println("ivarBed and alignerRefPrefix only work in illumina mode")
       System.exit(1)
   }
} else {
       println("Please select a workflow with --nanopolish, --illumina or --medaka")
       System.exit(1)
}

if ( ! params.prefix ) {
     println("Please supply a prefix for your output files with --prefix")
     System.exit(1)
}


// main workflow
workflow {
   if ( params.illumina ) {
       if (params.cram) {
           Channel.fromPath( "${params.directory}/**.cram" )
              .set{ ch_cramDirectory }
       }
       else {
	   Channel.fromFilePairs( params.fastqSearchPath, flat: true)
	       .set{ ch_filePairs }
       }
   }
   else {
       Channel.fromPath( "${params.basecalled_fastq}" )
              .set{ ch_runDirectory }

       Channel.fromPath( "${params.fast5_pass}" )
              .set{ ch_fast5Pass }

       Channel.fromPath( "${params.sequencing_summary}" )
              .set{ ch_seqSummary }

 
       // Check to see if we have barcodes
       def nanoporeBarcodeDirs = new FileNameByRegexFinder().getFileNames(params.basecalled_fastq, /.*barcode\d\d$/)
       
       
       if( nanoporeBarcodeDirs ) {
            // Yes, barcodes!
            Channel.fromPath( "${params.basecalled_fastq}/barcode*", type: 'dir', maxDepth: 1 )
                   .filter{ it.listFiles().size() > 5 }
                   .set{ ch_fastqDirs }
       } else {
            // No, no barcodes
            Channel.fromPath( "${params.basecalled_fastq}", type: 'dir', maxDepth: 1 )
                    .set{ ch_fastqDirs }
      }
   }

   main:
     if ( params.nanopolish || params.medaka ) {
         articNcovNanopore(ch_fastqDirs, ch_fast5Pass, ch_seqSummary)
     } else if ( params.illumina ) {
         if ( params.cram ) {
            ncovIlluminaCram(ch_cramDirectory)
         }
         else {
            ncovIllumina(ch_filePairs)
         }
     } else {
         println("Please select a workflow with --nanopolish, --illumina or --medaka")
     }
     
}

