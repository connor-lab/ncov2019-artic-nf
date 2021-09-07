#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {pango} from '../modules/analysis.nf'
include {nextclade} from '../modules/analysis.nf'
include {download_nextclade_files} from '../modules/analysis.nf'
include {getVariantDefinitions} from '../modules/analysis.nf'
include {aln2type} from '../modules/analysis.nf'
include {makeReport} from '../modules/analysis.nf'
include {articDownloadScheme} from '../modules/artic.nf'
include {FN4_upload} from '../modules/analysis.nf'

workflow ncovAnalysis {
    take:
       ch_consensusFiles
    
    main: 
    articDownloadScheme()

    downstreamAnalysis(ch_consensusFiles, articDownloadScheme.out.reffasta,articDownloadScheme.out.bed)
}

workflow downstreamAnalysis {
    take:
      consensus
      ch_preparedRef
      ch_bedFile

    main:
    println "Git info: $workflow.repository - $workflow.revision [$workflow.commitId]"

    pango(consensus)

    download_nextclade_files()

    nextclade(consensus.combine(ch_preparedRef).combine(download_nextclade_files.out))

    getVariantDefinitions()

    aln2type(consensus.combine(getVariantDefinitions.out.defs).combine(ch_preparedRef).combine(ch_bedFile))

    makeReport(pango.out.combine(aln2type.out, by:0).combine(nextclade.out.tsv,by:0))

    makeReport.out.tsv.collectFile(name:'analysisReport.tsv',
              storeDir:"${params.outdir}/analysis/report/${params.prefix}" ,
              keepHeader:true,
              skip:1)

    FN4_upload(consensus.combine(getVariantDefinitions.out).combine(ch_preparedRef).combine(ch_bedFile))

}
    
