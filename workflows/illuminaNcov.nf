#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {articDownloadScheme } from '../modules/artic.nf' 
include {readTrimming} from '../modules/illumina.nf' 
include {indexReference} from '../modules/illumina.nf'
include {readMapping} from '../modules/illumina.nf' 
include {trimPrimerSequences} from '../modules/illumina.nf' 
include {callVariants} from '../modules/illumina.nf'
include {makeConsensus} from '../modules/illumina.nf' 
include {callConsensusFreebayes} from '../modules/illumina.nf'

include {pangolinTyping} from '../modules/typing.nf' 
include {nextclade} from '../modules/typing.nf'
include {getVariantDefinitions} from '../modules/analysis.nf'
include {makeReport} from '../modules/analysis.nf'
include {cramToFastq} from '../modules/illumina.nf'

include {makeQCCSV} from '../modules/qc.nf'
include {writeQCSummaryCSV} from '../modules/qc.nf'
include {fastqc} from '../modules/qc.nf'
include {statsCoverage} from '../modules/qc.nf'
include {statsInsert} from '../modules/qc.nf'
include {statsAlignment} from '../modules/qc.nf'
include {multiqc} from '../modules/qc.nf'


include {bamToCram} from '../modules/out.nf'

include {collateSamples} from '../modules/upload.nf'

// import subworkflows
include {Genotyping} from './typing.nf'

workflow prepareReferenceFiles {
    // Get reference fasta
    if (params.ref) {
      Channel.fromPath(params.ref)
              .set{ ch_refFasta }
    } else {
      articDownloadScheme()
      articDownloadScheme.out.reffasta
                          .set{ ch_refFasta }
    }


    /* Either get BWA aux files from reference 
       location or make them fresh */
    
    if (params.ref) {
      // Check if all BWA aux files exist, if not, make them
      bwaAuxFiles = []
      refPath = new File(params.ref).getAbsolutePath()
      new File(refPath).getParentFile().eachFileMatch( ~/.*.bwt|.*.pac|.*.ann|.*.amb|.*.sa/) { bwaAuxFiles << it }
     
      if ( bwaAuxFiles.size() == 5 ) {
        Channel.fromPath( bwaAuxFiles )
               .set{ ch_bwaAuxFiles }

        ch_refFasta.combine(ch_bwaAuxFiles.collect().toList())
                   .set{ ch_preparedRef }
      } else {
        indexReference(ch_refFasta)
        indexReference.out
                      .set{ ch_preparedRef }
      }
    } else {
      indexReference(ch_refFasta)
      indexReference.out
                    .set{ ch_preparedRef }
    }
  
    /* If bedfile is supplied, use that,
       if not, get it from ARTIC github repo */ 
 
    if (params.bed ) {
      Channel.fromPath(params.bed)
             .set{ ch_bedFile }

    } else {
      articDownloadScheme.out.bed
                         .set{ ch_bedFile }
    }

    emit:
      bwaindex = ch_preparedRef
      bedfile = ch_bedFile
      reffasta = ch_refFasta
}


workflow sequenceAnalysis {
    take:
      ch_filePairs
      ch_preparedRef
      ch_bedFile

    main:
      fastqc(ch_filePairs)

      readTrimming(ch_filePairs)

      readMapping(readTrimming.out.trim.combine(ch_preparedRef))

      trimPrimerSequences(readMapping.out.combine(ch_bedFile))

      freebayes_out = callConsensusFreebayes(trimPrimerSequences.out.ptrim.combine(ch_preparedRef.map{ it[0] }))     
      freebayes_consensus_out = freebayes_out[0]

      callVariants(trimPrimerSequences.out.ptrim.combine(ch_preparedRef.map{ it[0] })) 

      makeConsensus(trimPrimerSequences.out.ptrim)

      makeQCCSV(trimPrimerSequences.out.ptrim.join(makeConsensus.out.consensus_fasta, by: 0)
                                   .combine(ch_preparedRef.map{ it[0] }))

      makeQCCSV.out.csv.splitCsv()
                       .unique()
                       .branch {
                           header: it[-1] == 'qc_pass'
                           fail: it[-1] == 'FALSE'
                           pass: it[-1] == 'TRUE'
    		       }
                       .set { qc }

      writeQCSummaryCSV(qc.header.concat(qc.pass).concat(qc.fail).toList())

      statsCoverage(trimPrimerSequences.out.ptrim.combine(ch_preparedRef.map{ it[0] }))

      statsInsert(trimPrimerSequences.out.ptrim.combine(ch_preparedRef.map{ it[0] }))

      statsAlignment(readMapping.out.combine(ch_preparedRef.map{ it[0] }))

      multiqc(fastqc.out.collect(), readTrimming.out.log.collect(), statsCoverage.out.collect(),
              statsInsert.out.stats.collect(), statsAlignment.out.collect())
      
      collateSamples(qc.pass.map{ it[0] }
                           .join(makeConsensus.out, by: 0)
                           .join(trimPrimerSequences.out.mapped))     
      
      pangolinTyping(makeConsensus.out.consensus_fasta)
      
      nextclade(makeConsensus.out.consensus_fasta)

      getVariantDefinitions()

      makeReport(pangolinTyping.out.combine(nextclade.out,by:0))

      makeReport.out.tsv.collectFile(name:'analysisReport.tsv',
		  storeDir:"${params.outdir}/AnalysisReport/${params.prefix}" , 
		  keepHeader:true,
		  skip:1)


      if (params.outCram) {
        bamToCram(trimPrimerSequences.out.mapped.map{it[0] } 
                        .join (trimPrimerSequences.out.ptrim.combine(ch_preparedRef.map{ it[0] })) )

      }

    emit:
      qc_pass = collateSamples.out
      variants = callVariants.out.variants

}

workflow ncovIllumina {
    take:
      ch_filePairs

    main:
      // Build or download fasta, index and bedfile as required
      prepareReferenceFiles()
      
      // Actually do analysis
      sequenceAnalysis(ch_filePairs, prepareReferenceFiles.out.bwaindex, prepareReferenceFiles.out.bedfile)

      // Do some typing if we have the correct files
      if ( params.gff ) {
          Channel.fromPath("${params.gff}")
                 .set{ ch_refGff }

          Channel.fromPath("${params.yaml}")
                 .set{ ch_typingYaml }

          Genotyping(sequenceAnalysis.out.variants, ch_refGff, prepareReferenceFiles.out.reffasta, ch_typingYaml) 

      }
}

workflow ncovIlluminaCram {
    take:
      ch_cramFiles
    main:
      // Convert CRAM to fastq
      cramToFastq(ch_cramFiles)

      // Run standard pipeline
      ncovIllumina(cramToFastq.out)
}

