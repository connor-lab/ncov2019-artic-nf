process collateSamples {
    tag { sampleName }

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

process uploadToCLIMB {
    tag { params.prefix }

    input:
    tuple(path(sshkey), path("${params.prefix}/*"))

    output:

    script:
    """
    rsync -Lav -e "ssh -i ${sshkey} -l ${params.CLIMBUser}" ${params.prefix} ${params.CLIMBHostname}:upload/
    """
}

