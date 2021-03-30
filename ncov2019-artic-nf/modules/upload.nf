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

process prepareUploadDirectory {
    tag { params.prefix }

    input:
    path("${params.prefix}/*")

    output:
    path("${params.prefix}")

    script:
    """
    echo "dummy" > dummyfile
    """
}


process uploadToCLIMB {
    tag { params.prefix }

    input:
    tuple(path(sshkey), path(uploadDir))

    output:

    script:
    """
    rsync -Lav -e "ssh -i ${sshkey} -l ${params.CLIMBUser}" ${uploadDir} ${params.CLIMBHostname}:upload/
    """
}

