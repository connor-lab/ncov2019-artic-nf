process makeQCCSV {
    tag { sampleName }

    input:
    tuple sampleName, path(bam), path(fasta), path(ref)

    output:
    path("${params.prefix}.${sampleName}.qc.csv")

    script:
    if ( params.illumina ) {
       qcSetting = "--illumina"
    } else {
       qcSetting = "--nanopore"
    }

    """
    qc.py ${qcSetting} --outfile ${params.prefix}.${sampleName}.qc.csv --sample ${sampleName} --ref ${ref} --bam ${bam} --fasta ${fasta}
    """
}


process writeQCSummaryCSV {
    tag { params.prefix }

    input:
    val lines

    exec:
    new File("${params.outdir}/${params.prefix}.qc.csv").withWriter { writer ->
        for ( line in lines ) {
            writer.writeLine(line.join(','))
         }   
    }
}

