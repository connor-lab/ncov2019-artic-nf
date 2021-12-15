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
    mv ${fasta} ${sampleName}
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

process uploadToBucket {
    tag {prefix}

    input:
    tuple( val(prefix), path("${prefix}.fasta"), path("${prefix}.bam"),"${prefix}.vcf")

    script:
    bucketName=params.uploadBucket
    """
    mkdir ${prefix}
    cp ${prefix}.fasta ${prefix}/
    cp "${prefix}.bam" ${prefix}/
    cp "${prefix}.vcf" ${prefix}/
    gzip ${prefix}/${prefix}.fasta
 
    oci os object put \
	-bn $bucketName \
	--force \
        --auth instance_principal \
	--file ${prefix}/${prefix}.fasta.gz

    oci os object put \
	-bn $bucketName \
	--force \
        --auth instance_principal \
	--file ${prefix}/${prefix}.bam

    oci os object put \
	-bn $bucketName \
	--force \
        --auth instance_principal \
	--file ${prefix}/${prefix}.bam
    """ 
}
