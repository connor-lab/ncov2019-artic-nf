// ARTIC processes

process articDownloadScheme{
    tag params.schemeRepoURL

    label 'internet'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "scheme", mode: "copy"

    output:
    path "scheme/**/${params.schemeVersion}/*.reference.fasta" , emit: reffasta
    path "scheme/**/${params.schemeVersion}/${params.scheme}.bed" , emit: bed
    path "scheme" , emit: scheme

    script:
    """
    git clone ${params.schemeRepoURL} scheme
    """
}

process articGuppyPlex {
    tag { samplePrefix + "-" + fastqDir }

    label 'largemem'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${params.prefix}*.fastq", mode: "copy"

    input:
    path(fastqDir)

    output:
    path "${params.prefix}*.fastq", emit: fastq

    script:
    """
    artic guppyplex \
    --min-length ${params.min_length} \
    --max-length ${params.max_length} \
    --prefix ${params.prefix} \
    --directory ${fastqDir}
    """
}

process articMinION {
    tag { sampleName }

    label 'largecpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.*", mode: "copy"

    input:
    tuple file(fastq), file(schemeRepo), file(fast5Pass), file(seqSummary)

    output:
    file("${sampleName}.*")
    
    tuple sampleName, file("${sampleName}.primertrimmed.rg.sorted.bam"), emit: ptrim
    tuple sampleName, file("${sampleName}.sorted.bam"), emit: mapped
    tuple sampleName, file("${sampleName}.consensus.fasta"), emit: consensus_fasta

    script:
    // Make an identifier from the fastq filename
    sampleName = fastq.getBaseName().replaceAll(~/\.fastq.*$/, '')

    // Configure artic minion pipeline
    minionRunConfigBuilder = []

    if (params.medaka) {
    minionRunConfigBuilder.add("--medaka")
    }

    if ( params.normalise ) {
    minionRunConfigBuilder.add("--normalise ${params.normalise}")
    }
    
    if ( params.bwa ) {
    minionRunConfigBuilder.add("--bwa")
    } else {
    minionRunConfigBuilder.add("--minimap2")
    }

    minionFinalConfig = minionRunConfigBuilder.join(" ")

    """
    artic minion ${minionFinalConfig} \
    --threads ${task.cpus} \
    --scheme-directory ${schemeRepo}/${params.schemeDir} \
    --read-file ${fastq} \
    --fast5-directory ${fast5Pass} \
    --sequencing-summary ${seqSummary} \
    ${params.scheme}/${params.schemeVersion} \
    ${sampleName}
    """
}

process articRemoveUnmappedReads {
    tag { sampleName }

    cpus 1

    input:
    tuple(sampleName, path(bamfile))

    output:
    tuple( sampleName, file("${sampleName}.mapped.sorted.bam"))

    script:
    """
    samtools view -F4 -o ${sampleName}.mapped.sorted.bam ${bamfile} 
    """
}

