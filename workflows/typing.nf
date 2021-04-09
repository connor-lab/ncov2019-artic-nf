
include {alignSeqs} from '../modules/typing.nf'
include {typeVariants} from '../modules/typing.nf'
include {mergeTypingCSVs} from '../modules/typing.nf'

workflow Genotyping {
    take:
      ch_consensus
      ch_refFasta
      ch_variantDefinitions

    main:
      alignSeqs(ch_consensus.combine(ch_refFasta))
      
      if ( params.gb ) {
        Channel.fromPath("${params.gb}")
               .set{ ch_gb }   
      } else {
        Channel.fromPath("${baseDir}/conf/dummyfile")
               .set{ ch_gb }
      }
      
      alignSeqs.out.map{ 
                        [ 
                          it[0], 
                          it[1].splitFasta(record: [id: true])[0]['id'], 
                          it[1] 
                          ]  
                          }.combine(ch_variantDefinitions)
                           .combine(ch_gb)
                           .set{ch_alignedSeqs}
      
      typeVariants(ch_alignedSeqs)
      
      mergeTypingCSVs(typeVariants.out.typing_csv.toList().map{ [ it ] }.combine(typeVariants.out.variants_csv.toList().map{ [ it ] }))
}
