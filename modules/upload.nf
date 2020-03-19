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
    tag { params.runPrefix }

    input:
    tuple(path(sshkey), path("${params.runPrefix}/*"))

    output:

    script:
    """
    echo "rsync -Lav -e "ssh -i ${sshkey} -l ${params.CLIMBUser}" ${params.runPrefix} ${params.CLIMBHostname}:upload/"
    """
}

