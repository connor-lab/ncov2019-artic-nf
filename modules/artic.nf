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

process articGuppyPlex {
    tag { params.runPrefix + "_" + fastqDir }

    label 'largemem'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${params.runPrefix}*.fastq", mode: "copy"

    input:
    tuple( path(runDirectory), fastqDir )

    output:
    path "${params.runPrefix}*.fastq", emit: fastq

    script:
    
    if ( fastqDir =~ /barcode/ )   
    """
    ln -s ${runDirectory}/fast5_pass .

    artic guppyplex \
    --min-length ${params.min_length} \
    --max-length ${params.max_length} \
    --prefix ${params.runPrefix} \
    --directory ${runDirectory}/fastq_pass/${fastqDir}
    """
    
    else if ( fastqDir =~ /fastq_pass/ )
    """
    ln -s ${runDirectory}/fast5_pass .

    artic guppyplex \
    --min-length ${params.min_length} \
    --max-length ${params.max_length} \
    --prefix ${params.runPrefix} \
    --directory ${runDirectory}/fastq_pass
    """
 
}

process articMinION {
    tag { sampleName }

    label 'largecpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.*", mode: "copy"
    publishDir "${params.outdir}/climb_upload/${params.runPrefix}/${sampleName}", pattern: "${sampleName}.consensus.fasta", mode: 'copy'

    input:
    tuple file(fastq), file(schemeRepo), file(runDirectory)

    output:
    file("${sampleName}.*")
    tuple sampleName, file("${sampleName}.primertrimmed.rg.sorted.bam"), emit: sorted_bam
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
    --fast5-directory ${runDirectory}/fast5_pass \
    --sequencing-summary ${runDirectory}/*cing_summary*.txt \
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
    tuple( sampleName, file("${sampleName}.mapped.primertrimmed.sorted.bam"))

    script:
    """
    samtools view -F4 -o ${sampleName}.mapped.primertrimmed.sorted.bam ${bamfile} 
    """
}

