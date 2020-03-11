// ARTIC processes

process articDownloadScheme{
    tag params.schemeRepoURL

    label 'internet'

    publishDir "${params.outdir}/${task.process.split(":")[1]}", pattern: "scheme", mode: "copy"

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

    publishDir "${params.outdir}/${task.process.split(":")[1]}", pattern: "${params.runPrefix}_fastq_pass.fastq", mode: "copy"
    publishDir "${params.outdir}/${task.process.split(":")[1]}", pattern: "${params.runPrefix}_sequencing_summary.txt", mode: "copy"

    input:
    file(runDirectory)

    output:
    tuple file("${params.runPrefix}_fastq_pass.fastq"), file("${params.runPrefix}_sequencing_summary.txt"), emit: gathered
    path "${params.runPrefix}_fastq_pass.fastq", emit: fastq

    script:
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

    cpus 4

    publishDir "${params.outdir}/${task.process.split(":")[1]}", pattern: "${params.runPrefix}_pass_NB*.fastq", mode: "copy"

    input:
    tuple file(fastqPass), file(sequencingSummary)

    output:
    file "${params.runPrefix}_pass_NB*.fastq", emit: demultiplexed

    script:
    """
    artic demultiplex --threads ${task.cpus} ${fastqPass}
    """
}

process nanopolishIndex {
   tag params.runPrefix

   label 'largemem'

   cpus 1

   publishDir "${params.outdir}/${task.process.split(":")[1]}", pattern: "${fastqPass}.index*", mode: "copy"

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


process articMinION {
    tag { sampleName }

    cpus 10

    publishDir "${params.outdir}/${task.process.split(":")[1]}", pattern: "${sampleName}.*", mode: "copy"

    input:
    tuple file("nanopolish/*"), file(bcFastqPass), file("nanopolish/*"), file(schemeRepo), file(runDirectory)

    output:

    script:
    if ( bcFastqPass =~ /.*_NB\d{2}.fastq$/ ) {
        sampleName = ( bcFastqPass =~ /.*(NB\d{2}).fastq$/ )[0][0]
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
            --normalise ${params.normalise}
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
