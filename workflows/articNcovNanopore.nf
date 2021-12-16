// ARTIC ncov workflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {articDownloadScheme} from '../modules/artic.nf' 
include {articGuppyPlex} from '../modules/artic.nf' 
include {articMinIONNanopolish} from  '../modules/artic.nf' 
include {articMinIONMedaka} from  '../modules/artic.nf'
include {articRemoveUnmappedReads} from '../modules/artic.nf' 
include {splitSeqSum} from '../modules/artic.nf' 
include {getObjFilesONT} from '../modules/artic'
include {download_primers} from '../modules/analysis.nf'
include {articMinIONViridian} from '../modules/artic'
include {viridianONTPrimers} from '../modules/viridian.nf'
include {viridianONTAuto} from '../modules/viridian.nf'


include {makeQCCSV} from '../modules/qc.nf'
include {writeQCSummaryCSV} from '../modules/qc.nf'

include {bamToCram} from '../modules/out.nf'

include {collateSamples} from '../modules/upload.nf'

//analysis
include {pango} from '../modules/analysis.nf'
include {nextclade} from '../modules/analysis.nf'
include {getVariantDefinitions} from '../modules/analysis.nf'
include {aln2type} from '../modules/analysis.nf'
include {makeReport} from '../modules/analysis.nf'
include {uploadToBucket} from '../modules/upload.nf'

// import subworkflows
include {CLIMBrsync} from './upload.nf'
include {Genotyping} from './typing.nf'
include {downstreamAnalysis} from './analysis.nf'

// workflow component for artic pipeline
workflow sequenceAnalysisNanopolish {
    take:
      ch_runFastqDirs
      ch_fast5Pass
      ch_seqSummary
    
    main:
      articDownloadScheme()
      
      articGuppyPlex(ch_runFastqDirs.flatten())

      splitSeqSum(ch_seqSummary)

      splitSeqSum.out.flatten()
                             .map {file -> tuple(file.baseName, file) }
			     .set{ barcodeSeqSums }

      articMinIONNanopolish(articGuppyPlex.out
                                          .combine(barcodeSeqSums, by:0)
                                          .combine(articDownloadScheme.out.scheme)
                                          .combine(ch_fast5Pass, by:0) )

      articRemoveUnmappedReads(articMinIONNanopolish.out.mapped)

      makeQCCSV(articMinIONNanopolish.out.ptrim
                                     .join(articMinIONNanopolish.out.consensus_fasta, by: 0)
                                     .combine(articDownloadScheme.out.reffasta))

      makeQCCSV.out.csv.splitCsv()
                       .unique()
                       .branch {
                           header: it[-1] == 'qc_pass'
                           fail: it[-1] == 'FALSE'
                           pass: it[-1] == 'TRUE'
                       }
                       .set { qc }

     writeQCSummaryCSV(qc.header.concat(qc.pass).concat(qc.fail).toList())

     collateSamples(qc.pass.map{ it[0] }
                           .join(articMinIONNanopolish.out.consensus_fasta, by: 0)
                           .join(articRemoveUnmappedReads.out))

     if (params.outCram) {
        bamToCram(articMinIONNanopolish.out.ptrim.map{ it[0] } 
                        .join (articDownloadScheme.out.reffasta.combine(ch_preparedRef.map{ it[0] })) )

      }


    emit:
      qc_pass = collateSamples.out
      reffasta = articDownloadScheme.out.reffasta
      vcf = articMinIONNanopolish.out.vcf

}

workflow sequenceAnalysisMedaka {
    take:
      ch_runFastqDirs

    main:
      articDownloadScheme()

      if (params.objstore || params.catsup) {
          getObjFilesONT(ch_runFastqDirs)
          articGuppyPlex(getObjFilesONT.out.fqs)
      }
      else {
          articGuppyPlex(ch_runFastqDirs.flatten())
      }

      articMinIONMedaka(articGuppyPlex.out
                                      .combine(articDownloadScheme.out.scheme))

      articRemoveUnmappedReads(articMinIONMedaka.out.mapped)

      makeQCCSV(articMinIONMedaka.out.ptrim.join(articMinIONMedaka.out.consensus_fasta, by: 0)
                           .combine(articDownloadScheme.out.reffasta))

      makeQCCSV.out.csv.splitCsv()
                       .unique()
                       .branch {
                           header: it[-1] == 'qc_pass'
                           fail: it[-1] == 'FALSE'
                           pass: it[-1] == 'TRUE'
                       }
                       .set { qc }

     writeQCSummaryCSV(qc.header.concat(qc.pass).concat(qc.fail).toList())

     collateSamples(qc.pass.map{ it[0] }
                           .join(articMinIONMedaka.out.consensus_fasta, by: 0)
                           .join(articRemoveUnmappedReads.out))

     downstreamAnalysis(articMinIONMedaka.out.consensus_fasta, articMinIONMedaka.out.vcf, articDownloadScheme.out.reffasta, articDownloadScheme.out.bed)     

     if (params.outCram) {
        bamToCram(articMinIONMedaka.out.ptrim.map{ it[0] } 
                        .join (articDownloadScheme.out.reffasta.combine(ch_preparedRef.map{ it[0] })) )

      }
    emit:
      qc_pass = collateSamples.out
      reffasta = articDownloadScheme.out.reffasta
      vcf = articMinIONMedaka.out.vcf

}

workflow sequenceAnalysisViridian {
    take:
      ch_runFastqDirs

    main:
      articDownloadScheme()
      

      getObjFilesONT(ch_runFastqDirs)


      if (params.primers == 'auto') {
      viridianONTAuto(getObjFilesONT.out.fqs
                                      .combine(articDownloadScheme.out.scheme))
      viridian=viridianONTAuto
      }
      else if (params.primers != 'auto') {
      download_primers(params.primers)
      viridianONTPrimers(getObjFilesONT.out.fqs
                                      .combine(articDownloadScheme.out.scheme)
                                      .combine(download_primers.out))
      viridian=viridianONTPrimers
      }

      // analysis
      downstreamAnalysis(viridian.out.consensus, viridian.out.vcfs,articDownloadScheme.out.reffasta,articDownloadScheme.out.bed)     
      
      if (params.uploadBucket != false) {
         uploadToBucket(viridian.out.consensus.combine(viridian.out.bam, by:0)
				.combine(viridian.out.vcfs, by:0))
      } 
}

workflow articNcovNanopore {
    take:
      ch_fastqDirs
    
    main:
      if ( params.nanopolish ) {
          Channel.fromPath( "${params.fast5_pass}/*/*.fast5" )
                 .map {file -> tuple(file.parent.name, file) }
                 .groupTuple(by: 0)
                 .set{ ch_fast5Pass }

          Channel.fromPath( "${params.sequencing_summary}" )
                 .set{ ch_seqSummary }

          sequenceAnalysisNanopolish(ch_fastqDirs, ch_fast5Pass, ch_seqSummary)

          sequenceAnalysisNanopolish.out.vcf.set{ ch_nanopore_vcf }

          sequenceAnalysisNanopolish.out.reffasta.set{ ch_nanopore_reffasta }

      } else if ( params.medaka ) {
          sequenceAnalysisMedaka(ch_fastqDirs)

          sequenceAnalysisMedaka.out.vcf.set{ ch_nanopore_vcf }

          sequenceAnalysisMedaka.out.reffasta.set{ ch_nanopore_reffasta }
           
      } else if ( params.viridian ) {
          sequenceAnalysisViridian(ch_fastqDirs)
      }

      if ( params.gff ) {
          Channel.fromPath("${params.gff}")
                 .set{ ch_refGff }

          Channel.fromPath("${params.yaml}")
                 .set{ ch_typingYaml }

          Genotyping(ch_nanopore_vcf, ch_refGff, ch_nanopore_reffasta, ch_typingYaml)

      }
      if ( params.upload ) {

        Channel.fromPath("${params.CLIMBkey}")
               .set{ ch_CLIMBkey }

        CLIMBrsync(sequenceAnalysis.out.qc_pass, ch_CLIMBkey )
      }
}

