
process copyIridaSamples {

    tag { fastq }

    publishDir "${params.outdir}", pattern: "irida_upload", mode: "symlink"

    input:
    file(fastq)

    output:
    path("irida_upload")

    script:
    """
    mkdir -p irida_upload
    mv ${fastq} irida_upload
    """
}

process generateIridaReport {

    publishDir "${params.outdir}/irida_upload", pattern: "SampleList.csv", mode: "copy"

    //conda 'environments/extras.txt'

    input:
    path(sample_upload)
    file(samplecsv)

    output:
    file("SampleList.csv")

    script:
    """
    irida_samples.py --sample_info ${samplecsv} --prefix ${params.prefix} --sample_dir ${sample_upload}
    """
}