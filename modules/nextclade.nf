process nextclade {
    label 'nextclade'
    publishDir "${params.outdir}/${params.lineagedir}/${sampleName}/", mode: 'copy', pattern: "${sampleName}_clade.tsv"
    input:
        tuple val(sampleName), path(consensus_fasta)
    output:
        tuple val(sampleName), path("${sampleName}_clade.tsv")
    script:
    """
    nextclade --input-fasta ${consensus_fasta} --output-tsv tmp.tsv
    cat tmp.tsv | tr -d "\r" > ${sampleName}_clade.tsv
    """
}
