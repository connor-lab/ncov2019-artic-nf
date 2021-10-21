process getVariantDefinitions {
    
    output:
    path('variant_definitions') 

    script:
    """
    git clone https://github.com/phe-genomics/variant_definitions
    """
}

process makeReport {
    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: 'copy', pattern: "${sampleName}_report.tsv"

    input:
    tuple(sampleName, path('pangolinTyping.csv'), path('nextclade_tree.json'), path('nextclade.tsv'),
		path('nextclade.json'))

    output:
    path "${sampleName}_report.tsv", emit: tsv

    script:
    """
    makeReport.py ${sampleName}
    """
}
