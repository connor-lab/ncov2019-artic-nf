process pango {
    tag { sampleName }

    publishDir "${params.outdir}/analysis/pango/${params.prefix}", mode: 'copy'

    input:
    tuple(sampleName,  path(fasta))

    output:
    tuple(sampleName, path("${sampleName}_lineage_report.csv"))

    script:
    """
    pangolin ${fasta}
    mv lineage_report.csv ${sampleName}_lineage_report.csv
    """
}


process nextclade {
    tag { sampleName }

    publishDir "${params.outdir}/analysis/nextclade/${params.prefix}", mode: 'copy'

    input:
    tuple(sampleName,  path(fasta))

    output:
    tuple(sampleName, path("${sampleName}_tree.json"),
	path("${sampleName}.tsv"),path("${sampleName}.json"))

    script:
    """
    nextclade --input-fasta ${fasta} \
        --output-tree ${sampleName}_tree.json \
        --output-tsv ${sampleName}.tsv \
        --output-json ${sampleName}.json
    """
}
