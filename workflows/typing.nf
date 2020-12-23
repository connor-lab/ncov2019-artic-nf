
include {typeVariants} from '../modules/typing.nf'
include {mergeTypingCSVs} from '../modules/typing.nf'

workflow Genotyping {
    take:
      ch_variantCalls
      ch_refGff
      ch_refFasta
      ch_typingYaml

    main:
      typeVariants(ch_variantCalls.combine(ch_refGff).combine(ch_refFasta).combine(ch_typingYaml))
      mergeTypingCSVs(typeVariants.out.typing_csv.toList().map{ [ it ] }.combine(typeVariants.out.variants_csv.toList().map{ [ it ] }))
}
