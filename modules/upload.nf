process collateSamples {
    tag { sampleName }

    publishDir "${params.outdir}/qc_pass_climb_upload/${params.prefix}", pattern: "${sampleName}", mode: 'copy'

    input:
    tuple(sampleName, path(bam), path(fasta))

    output:
    path("${sampleName}")

    script:
    """
    mkdir ${sampleName}
    mv ${bam} ${fasta} ${sampleName}
    """
}
