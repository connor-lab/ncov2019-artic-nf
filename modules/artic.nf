// ARTIC processes

process articDownloadScheme{
    tag params.schemeRepoURL

    label 'internet'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "scheme", mode: "copy"

    output:
    file("scheme")

    script:
    """
    git clone ${params.schemeRepoURL} scheme
    """
}

process articGather {
    tag params.runPrefix

    label 'largemem'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${params.runPrefix}_fastq_pass.fastq", mode: "copy"
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${params.runPrefix}_sequencing_summary.txt", mode: "copy"

    input:
    file(runDirectory)

    output:
    tuple file("${params.runPrefix}_fastq_pass.fastq"), file("${params.runPrefix}_sequencing_summary.txt"), emit: gathered
    path "${params.runPrefix}_fastq_pass.fastq", emit: fastq

    script:
    if ( params.barcode ) 
        """
        artic gather \
        --min-length ${params.min_length} \
        --max-length ${params.max_length} \
        --prefix ${params.runPrefix} \
        --directory ${runDirectory}

        cat ${params.runPrefix}_barcode*.fastq ${params.runPrefix}_unclassified.fastq > ${params.runPrefix}_fastq_pass.fastq
        """
    else
        """
        artic gather \
        --min-length ${params.min_length} \
        --max-length ${params.max_length} \
        --prefix ${params.runPrefix} \
        --directory ${runDirectory}
        """

 
}

process articDemultiplex {
    tag params.runPrefix

    label 'largecpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${params.runPrefix}_fastq_pass-NB*.fastq", mode: "copy"

    input:
    tuple file(fastqPass), file(sequencingSummary)

    output:
    file "${params.runPrefix}_fastq_pass-NB*.fastq"

    script:
    """
    artic demultiplex --threads ${task.cpus} ${fastqPass}
    """
}

process nanopolishIndex {
   tag params.runPrefix

   label 'largemem'

   cpus 1

   publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${fastqPass}.index*", mode: "copy"

   input:
   tuple file(fastqPass), file(sequencingSummary), file(runDirectory)

   output:
   file "${fastqPass}.index*"

   script:
   """
   ln -s ${runDirectory}/fast5_pass .
   nanopolish index -s ${sequencingSummary} -d fast5_pass ${fastqPass}
   """
}


process articMinIONNanopolish {
    tag { sampleName }

    label 'largecpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.*", mode: "copy"
    publishDir "${params.outdir}/climb_upload/${params.runPrefix}/${sampleName}", pattern: "${sampleName}.consensus.fasta", mode: 'copy'

    input:
    tuple file("nanopolish/*"), file(bcFastqPass), file("nanopolish/*"), file(schemeRepo), file(runDirectory)

    output:
    file("${sampleName}.*")
    tuple sampleName, file("${sampleName}.primertrimmed.sorted.bam"), emit: sorted_bam

    script:
    if ( bcFastqPass =~ /.*NB\d{2}.fastq$/ ) {
        sampleName = params.runPrefix + "-" + ( bcFastqPass =~ /.*(NB\d{2}).fastq$/ )[0][1]
    } else {
        sampleName = params.runPrefix
    }

    if ( params.normalise )
        if ( params.minimap )
            """
            ln -s ${runDirectory}/fast5_pass .                       

            artic minion --minimap \
            --normalise ${params.normalise} \
            --threads ${task.cpus} \
            --scheme-directory ${schemeRepo}/${params.schemeDir} \
            --read-file ${bcFastqPass} \
            --nanopolish-read-file nanopolish/${params.runPrefix}_fastq_pass.fastq \
            ${params.scheme}/${params.schemeVersion} \
            ${sampleName}
            """
        else 
            """
            ln -s ${runDirectory}/fast5_pass .            
 
            artic minion \
            --normalise ${params.normalise} \
            --threads ${task.cpus} \
            --scheme-directory ${schemeRepo}/${params.schemeDir} \
            --read-file ${bcFastqPass} \
            --nanopolish-read-file nanopolish/${params.runPrefix}_fastq_pass.fastq \
            ${params.scheme}/${params.schemeVersion} \
            ${sampleName}
            """
    else
        if ( params.minimap )
            """
            ln -s ${runDirectory}/fast5_pass .
 
            artic minion --minimap \
            --threads ${task.cpus} \
            --scheme-directory ${schemeRepo}/${params.schemeDir} \
            --read-file ${bcFastqPass} \
            --nanopolish-read-file nanopolish/${params.runPrefix}_fastq_pass.fastq \
            ${params.scheme}/${params.schemeVersion} \
            ${sampleName}
            """
        else
            """
            ln -s ${runDirectory}/fast5_pass .

            artic minion --threads ${task.cpus} \
            --scheme-directory ${schemeRepo}/${params.schemeDir} \
            --read-file ${bcFastqPass} \
            --nanopolish-read-file nanopolish/${params.runPrefix}_fastq_pass.fastq \
            ${params.scheme}/${params.schemeVersion} \
            ${sampleName}
            """
}

process articMinIONMedaka {
    tag { sampleName }

    label 'largecpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.*", mode: "copy"
    publishDir "${params.outdir}/climb_upload/${params.runPrefix}/${sampleName}", pattern: "${sampleName}.consensus.fasta", mode: 'copy'

    input:
    tuple file(bcFastqPass), file(schemeRepo), file(runDirectory)

    output:
    file("${sampleName}.*")
    tuple sampleName, file("${sampleName}.primertrimmed.sorted.bam"), emit: sorted_bam

    script:
    if ( bcFastqPass =~ /.*NB\d{2}.fastq$/ ) {
        sampleName = params.runPrefix + "-" + ( bcFastqPass =~ /.*(NB\d{2}).fastq$/ )[0][1]
    } else {
        sampleName = params.runPrefix
    }

    if ( params.normalise )
        if ( params.minimap )
            """
            artic minion --medaka \
            --minimap \
            --normalise ${params.normalise} \
            --threads ${task.cpus} \
            --scheme-directory ${schemeRepo}/${params.schemeDir} \
            --read-file ${bcFastqPass} \
            ${params.scheme}/${params.schemeVersion} \
            ${sampleName}
            """
        else
            """
            artic minion --medaka \
            --normalise ${params.normalise} \
            --threads ${task.cpus} \
            --scheme-directory ${schemeRepo}/${params.schemeDir} \
            --read-file ${bcFastqPass} \
            ${params.scheme}/${params.schemeVersion} \
            ${sampleName}
            """
    else
        if ( params.minimap )
            """
            artic minion --medaka \
            --minimap \
            --threads ${task.cpus} \
            --scheme-directory ${schemeRepo}/${params.schemeDir} \
            --read-file ${bcFastqPass} \
            ${params.scheme}/${params.schemeVersion} \
            ${sampleName}
            """
        else
            """
            artic minion --medaka \
            --threads ${task.cpus} \
            --scheme-directory ${schemeRepo}/${params.schemeDir} \
            --read-file ${bcFastqPass} \
            ${params.scheme}/${params.schemeVersion} \
            ${sampleName}
            """
}

process articRemoveUnmappedReads {
    tag { sampleName }

    publishDir "${params.outdir}/climb_upload/${params.runPrefix}/${sampleName}", pattern: "${sampleName}.mapped.primertrimmed.sorted.bam", mode: 'copy'

    cpus 1

    input:
    tuple(sampleName, path(bamfile))

    output:
    file("${sampleName}.mapped.primertrimmed.sorted.bam")

    script:
    """
    samtools view -F4 -o ${sampleName}.mapped.primertrimmed.sorted.bam ${bamfile} 
    """
}

