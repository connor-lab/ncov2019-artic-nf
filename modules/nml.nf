process generateIridaReport {

    publishDir "${params.outdir}", pattern: "irida_upload", mode: "copy"

    //conda 'environments/extras.txt'

    input:
    file(fastq)
    file(samplecsv)

    output:
    path("irida_upload")

    script:
    """
    mkdir irida_upload
    mv ${fastq} irida_upload
    irida_samples.py --sample_info ${samplecsv} --prefix ${params.prefix} --sample_dir irida_upload
    """
}
