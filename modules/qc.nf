process makeQCCSV {
    tag { sampleName }

    publishDir "${params.outdir}/qc_plots", pattern: "${sampleName}.depth.png", mode: 'copy'

    input:
    tuple sampleName, path(bam), path(fasta), path(ref)

    output:
    path "${params.prefix}.${sampleName}.qc.csv", emit: csv
    path "${sampleName}.depth.png"

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
    file("${params.outdir}/${params.prefix}.qc.csv").withWriter { writer ->
        for ( line in lines ) {
            writer.writeLine(line.join(','))
         }   
    }
}

process fastqc {
    tag { sampleName }

    publishDir "${params.outdir}/fastqc", mode: 'copy', overwrite: true

    input:
    tuple(sampleName, path(forward), path(reverse))

    output:
    file "*fastqc*"

    """
    fastqc ${forward} ${reverse} --format fastq --threads ${task.cpus}
    """

}

process mappingStatistics {
    tag { sampleName }

    label 'largemem'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "*.txt", mode: 'copy'

    input:
        tuple(sampleName, path(bam), path(ref))

    output:
        path "*.txt"

    script:
    """
    picard ${params.picardJavaSettings} CollectWgsMetrics -I ${sampleName}.mapped.primertrimmed.sorted.bam \
    -O ${sampleName}.mapped.primertrimmed.sorted.metrics.txt -R ${ref} ${params.wgsMetricsOptions}      
    """
}

process multiqc {   
    tag { params.prefix }
    
    label 'largemem'
    
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    
    input:
    path trimLogList
    path mapLogList
    path qcLogList

    output:
    file '*multiqc.html'
    file '*multiqc_data/multiqc_data.json'

    script:
    """
    multiqc . --filename ${params.prefix}_multiqc.html --data-format json \
    ${params.multiqcOptions}
    """
}
