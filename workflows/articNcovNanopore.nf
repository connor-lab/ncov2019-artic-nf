// ARTIC ncov workflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {articGuppyPlex} from '../modules/artic.nf'
include {articMinIONMedaka} from  '../modules/artic.nf'
include {articRemoveUnmappedReads} from '../modules/artic.nf'

include {makeQCCSV} from '../modules/qc.nf'
include {writeQCSummaryCSV} from '../modules/qc.nf'

include {bamToCram} from '../modules/out.nf'

include {collateSamples} from '../modules/upload.nf'


// import subworkflows
include {Genotyping} from './typing.nf'


workflow sequenceAnalysisMedaka {
    take:
      ch_runFastqDirs

    main:

      // prepare reference fasta
      ch_reffasta = Channel.fromPath( "${params.schemeDir}/${params.scheme}/SARS-CoV-2/${params.schemeVersion}/*.reference.fasta" )
      ch_scheme = Channel.fromPath( "${params.schemeDir}/${params.scheme}" )

      articGuppyPlex(ch_runFastqDirs.flatten())

      articMinIONMedaka(articGuppyPlex.out.fastq
                                      .filter{ it.countFastq() > params.minReadsArticGuppyPlex }
                                      .combine(ch_scheme))

      articRemoveUnmappedReads(articMinIONMedaka.out.mapped)

      makeQCCSV(articMinIONMedaka.out.ptrim.join(articMinIONMedaka.out.consensus_fasta, by: 0)
                           .combine(ch_reffasta))

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

     if (params.outCram) {
        bamToCram(articMinIONMedaka.out.ptrim.map{ it[0] } 
                        .join (ch_reffasta.combine(ch_preparedRef.map{ it[0] })) )

      }
    emit:
      qc_pass = collateSamples.out
      reffasta = ch_reffasta
      vcf = articMinIONMedaka.out.vcf

}


workflow articNcovNanopore {
    take:
      ch_fastqDirs
    
    main:
      if ( params.medaka ) {
          sequenceAnalysisMedaka(ch_fastqDirs)

          sequenceAnalysisMedaka.out.vcf.set{ ch_nanopore_vcf }

          sequenceAnalysisMedaka.out.reffasta.set{ ch_nanopore_reffasta }
      }

      if ( params.gff ) {
          Channel.fromPath("${params.gff}")
                 .set{ ch_refGff }

          Channel.fromPath("${params.yaml}")
                 .set{ ch_typingYaml }

          Genotyping(ch_nanopore_vcf, ch_refGff, ch_nanopore_reffasta, ch_typingYaml)

      }
}

