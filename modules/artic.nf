// ARTIC processes

process articDownloadScheme{
    tag 'reference'

    label 'internet'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "scheme", mode: "copy"

    output:
    path "${params.schemeDir}/${params.scheme}/${params.schemeVersion}/*.reference.fasta" , emit: reffasta
    path "${params.schemeDir}/${params.scheme}/${params.schemeVersion}/*.primer.bed" , emit: bed
    path "${params.schemeDir}" , emit: scheme

    script:
    """
    git clone ${params.schemeRepoURL} ${params.schemeDir}
    """
}

process articGuppyPlex {
    tag { prefix }

    label 'largemem'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${params.prefix}*.fastq", mode: "copy"

    input:
    tuple(prefix, path(fastq))

    output:
    tuple val(prefix), path("${prefix}_.fastq"), emit: fastq

    script:
    """
    artic guppyplex \
    --min-length ${params.min_length} \
    --max-length ${params.max_length} \
    --prefix ${prefix} \
    --directory ./
    """
}

process articMinIONMedaka {
    tag { sampleName }

    cpus 4

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}*", mode: "copy"

    input:
    tuple sampleName, file(fastq), file(schemeRepo)

    output:
    file("${sampleName}*")
    
    tuple sampleName, file("${sampleName}.primertrimmed.rg.sorted.bam"), emit: ptrim
    tuple sampleName, file("${sampleName}.sorted.bam"), emit: mapped
    tuple sampleName, file("${sampleName}.consensus.fasta"), emit: consensus_fasta
    tuple sampleName, file("${sampleName}.pass.vcf.gz"), emit: vcf

    script:
    // Make an identifier from the fastq filename
    //sampleName = fastq.getBaseName().replaceAll(~/\.fastq.*$/, '')

    // Configure artic minion pipeline
    minionRunConfigBuilder = []

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
    artic minion --medaka \
    ${minionFinalConfig} \
    --threads ${task.cpus} \
    --scheme-directory ${schemeRepo} \
    --read-file ${fastq} \
    ${params.scheme}/${params.schemeVersion} \
    ${sampleName}
    """
}

process splitSeqSum {
    tag 'splitSeqSum'

    cpus 4

    input:
    file(seqSummary)

    output:
    file("barcodes/*.txt")

    script:
    """
    split_summary_by_barcode.py ${seqSummary} barcodes 
    """

}

process articMinIONNanopolish {
    tag { sampleName }

    cpus 4
    memory '3 GB'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}*", mode: "copy"

    input:
    tuple barcode, file(fastq), file(seqSummary), file(schemeRepo), file(fast5Pass)

    output:
    file("${sampleName}*")
    
    tuple sampleName, file("${sampleName}.primertrimmed.rg.sorted.bam"), emit: ptrim
    tuple sampleName, file("${sampleName}.sorted.bam"), emit: mapped
    tuple sampleName, file("${sampleName}.consensus.fasta"), emit: consensus_fasta
    tuple sampleName, file("${sampleName}.pass.vcf.gz"), emit: vcf

    script:
    // Make an identifier from the fastq filename
    sampleName = fastq.getBaseName().replaceAll(~/\.fastq.*$/, '')

    // Configure artic minion pipeline
    minionRunConfigBuilder = []

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
    --scheme-directory ${schemeRepo} \
    --read-file ${fastq} \
    --fast5-directory ./ \
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

process getObjFilesONT {
    /**
    * fetches fastq files from object store using OCI bulk download (https://docs.oracle.com/en-us/iaas/tools/oci-cli/2.24.4/oci_cli_docs/cmdref/os/object/bulk-download.html)
    * @input
    * @output
    */

    tag { prefix }

    input:
        tuple bucket, prefix

    output:
        tuple prefix, path("${prefix}.fastq.gz")

    script:
	db=params.krakdb
        """
	oci os object bulk-download \
		-bn $bucket \
		--download-dir ./ \
		--overwrite \
		--auth instance_principal \
		--prefix $prefix

	kraken2 -db ${db} \
		--memory-mapping \
		--report ${prefix}_summary.txt \
		--output ${prefix}_read_classification \
        	*.fastq.gz 

        awk '\$3==\"9606\" { print \$2 }' ${prefix}_read_classification >> kraken2_human_read_list
        awk '\$3!=\"9606\" { print \$2 }' ${prefix}_read_classification >> kraken2_nonhuman_read_list

        seqs=*.fastq.gz
        for seq in \${seqs}
        do
	    seqtk subseq \${seq} kraken2_nonhuman_read_list | gzip >> "${prefix}.fastq.gz"
	done
	"""
}

process articMinIONViridian {
    /**
    * runs viridian workflow https://github.com/iqbal-lab-org/viridian_workflow
    * @input
    * @output
    */

    tag { prefix }

    publishDir "${params.outdir}/viridian"

    input:
        tuple prefix, path("${prefix}.fastq.gz"),path(schemeRepo)

    output:
        tuple prefix, path("${prefix}_outdir/viridian/consensus.final_assembly.fa")

    script:
        """
        wget https://raw.githubusercontent.com/iqbal-lab-org/viridian_workflow/master/data/nCoV-artic-v3.bed
        viridian_workflow run_one_sample_ont \
		${schemeRepo}/nCoV-2019/V3/nCoV-2019.reference.fasta \
		nCoV-artic-v3.bed \
		${prefix}.fastq.gz \
		${prefix}_outdir/
        """
}

