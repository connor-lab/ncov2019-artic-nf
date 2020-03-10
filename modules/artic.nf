// ARTIC processes

process articDownloadScheme{
    tag params.schemeRepoURL

    publishDir "${params.outdir}/${task.process.split(":")[1]}", pattern: "scheme", mode: "copy"

    output:
    file("scheme")

    script:
    """
    git clone ${params.schemeRepoURL} scheme
    """
}

process articGather {
    tag params.runDirectory

    publishDir "${params.outdir}/${task.process.split(":")[1]}", pattern: "${params.runDirectory}_fastq_pass.fastq", mode: "copy"
    publishDir "${params.outdir}/${task.process.split(":")[1]}", pattern: "${params.runDirectory}_sequencing_summary.txt", mode: "copy"

    input:
    file(runDirectory)

    output:
    tuple file("${params.runDirectory}_fastq_pass.fastq"), file("${params.runDirectory}_sequencing_summary.txt"), emit: gathered
    path "${params.runDirectory}_fastq_pass.fastq", emit: fastq

    script:
    """
    artic gather \
    --min-length ${params.min_length} \
    --max-length ${params.max_length} \
    --prefix ${params.runDirectory} \
    --directory ${runDirectory} 
    """
}

process articDemultiplex {
    tag params.runDirectory

    cpus 4

    publishDir "${params.outdir}/${task.process.split(":")[1]}", pattern: "${params.runDirectory}_pass_NB*.fastq", mode: "copy"

    input:
    tuple file(fastqPass), file(sequencingSummary)

    output:
    file "${params.runDirectory}_pass_NB*.fastq", emit: demultiplexed

    script:
    """
    artic demultiplex --threads ${task.cpus} ${fastqPass}
    """
}

process nanopolishIndex {
   tag params.runDirectory

   cpus 1

   publishDir "${params.outdir}/${task.process.split(":")[1]}", pattern: "${fastqPass}.index*", mode: "copy"

   input:
   tuple file(fastqPass), file(sequencingSummary), file(runDirectory)

   output:
   file "${fastqPass}.index*"

   script:
   """
   nanopolish index -s ${sequencingSummary} -d ${runDirectory}/fast5_pass ${fastqPass}
   """
}


process articMinION {
    tag { sampleName }

    cpus 4

    publishDir "${params.outdir}/${task.process.split(":")[1]}", pattern: "${sampleName}*", mode: "copy"

    input:
    tuple file("nanopolish/*"), file(bcFastqPass), file("nanopolish/*"), file(schemeRepo), file(runDirectory)

    output:

    script:
    if ( bcFastqPass =~ /.*_NB\d{2}.fastq$/ ) {
        sampleName = ( bcFastqPass =~ /.*(NB\d{2}).fastq$/ )[0][0]
    } else {
        sampleName = params.runDirectory
    }

    if ( params.normalise )
        if ( params.minimap )
            """
            artic minion --minimap \
            --normalise ${params.normalise} \
            --threads ${task.cpus} \
            --scheme-directory ${schemeRepo}/${params.schemeDir} \
            --read-file ${bcFastqPass} \
            --nanopolish-read-file nanopolish/${params.runDirectory}_fastq_pass.fastq \
            ${params.scheme}/${params.schemeVersion} \
            ${sampleName}
            """
        else 
            """
            artic minion \
            --normalise ${params.normalise}
            --threads ${task.cpus} \
            --scheme-directory ${schemeRepo}/${params.schemeDir} \
            --read-file ${bcFastqPass} \
            --nanopolish-read-file nanopolish/${params.runDirectory}_fastq_pass.fastq \
            ${params.scheme}/${params.schemeVersion} \
            ${sampleName}
            """
    else
        if ( params.minimap )
            """
            artic minion --minimap \
            --threads ${task.cpus} \
            --scheme-directory ${schemeRepo}/${params.schemeDir} \
            --read-file ${bcFastqPass} \
            --nanopolish-read-file nanopolish/${params.runDirectory}_fastq_pass.fastq \
            ${params.scheme}/${params.schemeVersion} \
            ${sampleName}
            """
        else
            """
            artic minion --threads ${task.cpus} \
            --scheme-directory ${schemeRepo}/${params.schemeDir} \
            --read-file ${bcFastqPass} \
            --nanopolish-read-file nanopolish/${params.runDirectory}_fastq_pass.fastq \
            ${params.scheme}/${params.schemeVersion} \
            ${sampleName}
            """
}
