#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// include modules
include {printHelp} from './modules/help.nf'
include {makeFastqSearchPath} from './modules/util.nf'

// import subworkflows
include {articNcovNanopore} from './workflows/articNcovNanopore.nf' 
include {ncovIllumina} from './workflows/illuminaNcov.nf'
include {ncovIlluminaCram} from './workflows/illuminaNcov.nf'
include {ncovIlluminaObj} from './workflows/illuminaNcov.nf'
include {ncovAnalysis} from './workflows/analysis'

if (params.varCaller == 'medaka'){
    params.medaka=true
}
else if (params.varCaller == 'viridian' && !params.illumina){
    params.viridian=true
} else {
    params.viridian = false
}

if (params.help){
    printHelp()
    exit 0
}

if (params.profile){
    println("Profile should have a single dash: -profile")
    System.exit(1)
}

if ( params.illumina ) {
   if ( !params.directory && !params.objstore && !params.catsup) {
       println("Please supply a directory containing fastqs or CRAMs with --directory. Specify --cram if supplying a CRAMs directory")
       println("Use --help to print help")
       System.exit(1)
   }
   if ( (params.bed && ! params.ref) || (!params.bed && params.ref) ) {
       println("--bed and --ref must be supplied together")
       System.exit(1)
   }
} else if ( params.nanopolish ) {
   if (! params.basecalled_fastq ) {
       println("Please supply a directory containing basecalled fastqs with --basecalled_fastq. This is the output directory from guppy_barcoder or guppy_basecaller - usually fastq_pass. This can optionally contain barcodeXX directories, which are auto-detected.")
   }
   if (! params.fast5_pass ) {
       println("Please supply a directory containing fast5 files with --fast5_pass (this is the fast5_pass directory)")
   }
   if (! params.sequencing_summary ) {
       println("Please supply the path to the sequencing_summary.txt file from your run with --sequencing_summary")
       System.exit(1)
   }
   if ( params.bed || params.ref ) {
       println("ivarBed and alignerRefPrefix only work in illumina mode")
       System.exit(1)
   }
} else if ( params.medaka ) {
   if (! params.basecalled_fastq && !params.objstore && !params.catsup ) {
       println("Please supply a directory containing basecalled fastqs with --basecalled_fastq. This is the output directory from guppy_barcoder or guppy_basecaller - usually fastq_pass. This can optionally contain barcodeXX directories, which are auto-detected.")
   }
} else if ( params.viridian ) {
   if (! params.basecalled_fastq && !params.objstore && !params.catsup ) {
       println("Please supply a directory containing basecalled fastqs with --basecalled_fastq. This is the output directory from guppy_barcoder or guppy_basecaller - usually fastq_pass. This can optionally contain barcodeXX directories, which are auto-detected.")
   }
} else if ( params.analysis ) {
   if ( params.consensus_seqs ) {
       println("Please supply a directory containing consensus with --consensus_seqs")
}
} else {
       println("Please select a workflow with --nanopolish, --illumina or --medaka, or use --help to print help")
       System.exit(1)
}


if ( ! params.prefix ) {
     println("Please supply a prefix for your output files with --prefix")
     println("Use --help to print help")
     System.exit(1)
} else {
     if ( params.prefix =~ /\// ){
         println("The --prefix that you supplied contains a \"/\", please replace it with another character")
         System.exit(1)
     }
} 


Channel.from("$workflow.manifest.version")
	.view()
	.collectFile(name: "${params.outdir}/workflow_manifest_version.txt")

println "Manifest's pipeline version: $workflow.manifest.version"
println "Project : $workflow.projectDir"

// main workflow
workflow {
   if (params.analysis) {
           Channel.fromPath( "${params.consensus_seqs}/consensus_seqs/*.fasta")
                  .map { file -> tuple(file.simpleName, file) }
                  .view()
                  .set{ ch_consensusFiles }
       }
   else if ( params.illumina ) {
       if (params.cram) {
           Channel.fromPath( "${params.directory}/**.cram" )
                  .map { file -> tuple(file.baseName, file) }
                  .set{ ch_cramFiles }
       }
       else if (params.objstore != false) {
           Channel.fromPath( "${params.objstore}" )
                  .splitCsv()
                  .map { row -> tuple(row[0], row[1], row[1]) }
                  .set{ ch_objFiles }
       }
       else if (params.catsup != false) {
           Channel.fromPath( "${params.catsup}/sp3data.csv" )
                  .splitCsv(header: true)
                  .map { row -> tuple("${params.bucket}", "${row.submission_uuid4}/${row.sample_uuid4}", "${row.sample_uuid4}") }
                  .unique()
                  .set{ ch_objFiles }
       }
       else {
           fastqSearchPath = makeFastqSearchPath( params.illuminaPrefixes, params.illuminaSuffixes, params.fastq_exts )

	   Channel.fromFilePairs( fastqSearchPath, flat: true)
	          .filter{ !( it[0] =~ /Undetermined/ ) }
	          .set{ ch_filePairs }
       }
   }
   else {
       // Check to see if we have barcodes
       nanoporeBarcodeDirs = file("${params.basecalled_fastq}/barcode*", type: 'dir', maxdepth: 1 )
       nanoporeNoBarcode = file("${params.basecalled_fastq}/*.fastq", type: 'file', maxdepth: 1)

       if (params.objstore) {
           Channel.fromPath( "${params.objstore}" )
                  .splitCsv()
                  .map { row -> tuple(row[0], row[1], row[1]) }
                  .set{ ch_objFiles }
       }
       else if (params.catsup != false) {
           Channel.fromPath( "${params.catsup}/sp3data.csv" )
                  .splitCsv(header: true)
                  .map { row -> tuple("${params.bucket}", "${row.submission_uuid4}/${row.sample_uuid4}", "${row.sample_uuid4}") }
                  .view()
                  .unique()
                  .set{ ch_objFiles }
       }
       else{
       if( nanoporeBarcodeDirs ) {
            // Yes, barcodes!
            Channel.fromPath( nanoporeBarcodeDirs )
                   .filter( ~/.*barcode[0-9]{1,4}$/ )
                   .filter{ d ->
                            def count = 0
                            for (x in d.listFiles()) {
                                if (x.isFile()) {
                                    count += x.countFastq()
                                }
                            }
                            count > params.minReadsPerBarcode
                   }.set{ ch_fastqDirs }
       } else if ( nanoporeNoBarcode ){
            // No, no barcodes
            Channel.fromPath( "${params.basecalled_fastq}", type: 'dir', maxDepth: 1 )
                    .set{ ch_fastqDirs }
      } else {
            println("Couldn't detect whether your Nanopore run was barcoded or not. Use --basecalled_fastq to point to the unmodified guppy output directory.")
            System.exit(1)
      }
   }
   }

   main:
     if ( params.nanopolish || params.medaka || (params.viridian && !params.illumina )) {
	if ( params.objstore || params.catsup ) {
	     articNcovNanopore(ch_objFiles)
	}
	else {
         articNcovNanopore(ch_fastqDirs)
	}
     } else if ( params.illumina ) {
         if ( params.cram ) {
            ncovIlluminaCram(ch_cramFiles)
         }
         else if ( params.objstore || params.catsup ) {
            ncovIlluminaObj(ch_objFiles)
         }
         else {
            ncovIllumina(ch_filePairs)
         }
      } else if (params.analysis) {
            ncovAnalysis(ch_consensusFiles)
      } else {
         println("Please select a workflow with --nanopolish, --illumina or --medaka")
     }
     
}

