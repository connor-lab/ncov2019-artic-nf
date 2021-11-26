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

    publishDir "${params.outdir}/QCStats/${task.process.replaceAll(":","_")}", mode: 'copy', overwrite: true

    input:
    tuple sampleName, path(forward), path(reverse)

    output:
    file "*fastqc*"

    """
    fastqc ${forward} ${reverse} --format fastq --threads ${task.cpus}
    """
}

process statsCoverage {
    tag { sampleName }

    label 'largemem'

    publishDir "${params.outdir}/QCStats/${task.process.replaceAll(":","_")}", pattern: "*.txt", mode: 'copy'

    input:
    tuple sampleName, path(bam), path(ref)

    output:
    path "*.txt"

    """
    picard ${params.picardJavaSettings} CollectWgsMetrics -I ${sampleName}.mapped.primertrimmed.sorted.bam \
    -O ${sampleName}.mapped.primertrimmed.sorted.metrics.txt -R ${ref} ${params.wgsMetricsOptions}      
    """
}

process statsInsert {
    tag { sampleName }

    label 'largemem'

    publishDir "${params.outdir}/QCStats/${task.process.replaceAll(":","_")}", pattern: "*.txt", mode: 'copy'
    publishDir "${params.outdir}/QCStats/${task.process.replaceAll(":","_")}", pattern: "*.pdf", mode: 'copy'

    input:
    tuple sampleName, path(bam), path(ref)

    output:
    path "${sampleName}_insert_size.metrics.txt", emit: stats
    path "${sampleName}_insert_size.distribution.pdf", optional: true

    """
    mkdir ./tmp
    if [[ \$( samtools view -F 12 ${bam} | cut -f1 | sort -T ./tmp | uniq -c | awk '\$1 > 1 { count++ } END { print count }' ) > 0 ]]
    then
       picard CollectInsertSizeMetrics I=${bam} O=${sampleName}_insert_size.metrics.txt \
       H=${sampleName}_insert_size.distribution.pdf
    else
       echo "Skipping sample ${sampleName} - no usable paired reads"
       touch ${sampleName}_insert_size.metrics.txt
    fi
    """
}

process statsAlignment {
    tag { sampleName }

    publishDir "${params.outdir}/QCStats/${task.process.replaceAll(":","_")}", pattern: "*.txt", mode: 'copy'

    input:
    tuple sampleName, path(bam), path(ref)

    output:
    path "${sampleName}_alignment.metrics.txt"

    """
    picard CollectAlignmentSummaryMetrics R=${ref} I=${bam} O=${sampleName}_alignment.metrics.txt
    """
}

process multiqc {
    tag { params.prefix }
    
    label 'largemem'
    
    publishDir "${params.outdir}/QCStats/${task.process.replaceAll(":","_")}", mode: 'copy'
    
    input:
    path qcLogList
    path trimLogList
    path mapLogList
    path insertLogList
    path alignLogList

    output:
    file '*multiqc.html'
    file '*multiqc_data/multiqc_data.json'

    """
    multiqc . --filename ${params.prefix}_multiqc.html --data-format json \
    ${params.multiqcOptions}
    """
}
