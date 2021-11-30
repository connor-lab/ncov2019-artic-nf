process viridianPrimers {
    /**
    * runs viridian workflow https://github.com/iqbal-lab-org/viridian_workflow
    * @input
    * @output
    */

    tag { prefix }

    publishDir "${params.outdir}/consensus_seqs/", mode: 'copy', pattern: "*.fasta"
    publishDir "${params.outdir}/VCF/", mode: 'copy', pattern: "*.vcf"
    publishDir "${params.outdir}/qc/", mode: 'copy', pattern: "*.json"
    if (params.TESToutputMODE){
        publishDir "${params.outdir}/bam/", mode: 'copy', pattern: "*.bam"
    }

    input:
        tuple prefix, path("${prefix}_1.fastq.gz"), path("${prefix}_2.fastq.gz"),path('primers'), path('ref.fa'),path("*")

    output:
        tuple prefix, path("${prefix}.fasta"), emit: consensus
        tuple prefix, path("${prefix}.viridian_log.json"), emit: coverage
        tuple prefix, path("${prefix}.vcf"), emit: vcfs
        tuple prefix, path("${prefix}.bam"), emit: bam

 
    script:
    """
    viridian_workflow run_one_sample \
            --tech illumina \
            --ref_fasta ref.fa \
            --amplicon_json primers \
            --reads1 ${prefix}_1.fastq.gz \
            --reads2 ${prefix}_2.fastq.gz \
            --outdir ${prefix}_outdir/ \
            --sample_name ${prefix} \
            --keep_bam
    cp ${prefix}_outdir/consensus.fa ${prefix}.fasta
    cp ${prefix}_outdir/log.json ${prefix}.viridian_log.json
    cp ${prefix}_outdir/variants.vcf ${prefix}.vcf
    cp ${prefix}_outdir/reference_mapped.bam ${prefix}.bam
    """ 
}


process viridianAuto {
    /** 
    * runs viridian workflow https://github.com/iqbal-lab-org/viridian_workflow
    * @input
    * @output
    */

    tag { prefix }

    publishDir "${params.outdir}/consensus_seqs/", mode: 'copy', pattern: "*.fasta"
    publishDir "${params.outdir}/VCF/", mode: 'copy', pattern: "*.vcf"
    publishDir "${params.outdir}/qc/", mode: 'copy', pattern: "*.json"
    if (params.TESToutputMODE){
        publishDir "${params.outdir}/bam/", mode: 'copy', pattern: "*.bam"
    }   

    input:
        tuple prefix, path("${prefix}_1.fastq.gz"), path("${prefix}_2.fastq.gz"), path('ref.fa'),path("*")

    output:
        tuple prefix, path("${prefix}.fasta"), emit: consensus
        tuple prefix, path("${prefix}.viridian_log.json"), emit: coverage
        tuple prefix, path("${prefix}.vcf"), emit: vcfs
        tuple prefix, path("${prefix}.bam"), emit: bam


    script:
    """
    viridian_workflow run_one_sample \
            --tech illumina \
            --ref_fasta ref.fa \
            --reads1 ${prefix}_1.fastq.gz \
            --reads2 ${prefix}_2.fastq.gz \
            --outdir ${prefix}_outdir/ \
            --sample_name ${prefix} \
            --keep_bam
    cp ${prefix}_outdir/consensus.fa ${prefix}.fasta
    cp ${prefix}_outdir/log.json ${prefix}.viridian_log.json
    cp ${prefix}_outdir/variants.vcf ${prefix}.vcf
    cp ${prefix}_outdir/reference_mapped.bam ${prefix}.bam
    """ 
}


process viridianONTPrimers {
    /**
    * runs viridian workflow https://github.com/iqbal-lab-org/viridian_workflow
    * @input
    * @output
    */

    tag { prefix }

    publishDir "${params.outdir}/consensus_seqs/", mode: 'copy', pattern: "*.fasta"
    publishDir "${params.outdir}/VCF/", mode: 'copy', pattern: "*.vcf"
    publishDir "${params.outdir}/qc/", mode: 'copy', pattern: "*.json"
    if (params.TESToutputMODE){
	publishDir "${params.outdir}/bam/", mode: 'copy', pattern: "*.bam"
    }

    input:
        tuple prefix, path("${prefix}.fastq.gz"),path(schemeRepo),path('primers')

    output:
        tuple prefix, path("${prefix}.fasta"), emit: consensus
        tuple prefix, path("${prefix}.viridian_log.json"), emit: coverage
        tuple prefix, path("${prefix}.vcf"), emit: vcfs
	tuple prefix, path("${prefix}.bam"), emit: bam

    script:
        """
        viridian_workflow run_one_sample \
		--tech ont \
		--ref_fasta ${schemeRepo}/nCoV-2019/V3/nCoV-2019.reference.fasta \
		--amplicon_json primers \
		--reads ${prefix}.fastq.gz \
		--outdir ${prefix}_outdir/ \
		--sample_name ${prefix} \
		--keep_bam
        cp ${prefix}_outdir/consensus.fa ${prefix}.fasta
        cp ${prefix}_outdir/log.json ${prefix}.viridian_log.json
        cp ${prefix}_outdir/variants.vcf ${prefix}.vcf
	cp ${prefix}_outdir/reference_mapped.bam ${prefix}.bam
        """
}

process viridianONTAuto {
    /**
    * runs viridian workflow https://github.com/iqbal-lab-org/viridian_workflow
    * @input
    * @output
    */

    tag { prefix }

    publishDir "${params.outdir}/consensus_seqs/", mode: 'copy', pattern: "*.fasta"
    publishDir "${params.outdir}/VCF/", mode: 'copy', pattern: "*.vcf"
    publishDir "${params.outdir}/qc/", mode: 'copy', pattern: "*.json"
    if (params.TESToutputMODE){
	publishDir "${params.outdir}/bam/", mode: 'copy', pattern: "*.bam"
    }

    input:
        tuple prefix, path("${prefix}.fastq.gz"),path(schemeRepo)

    output:
        tuple prefix, path("${prefix}.fasta"), emit: consensus
        tuple prefix, path("${prefix}.viridian_log.json"), emit: coverage
        tuple prefix, path("${prefix}.vcf"), emit: vcfs
	tuple prefix, path("${prefix}.bam"), emit: bam

    script:
        """
        viridian_workflow run_one_sample \
		--tech ont \
		--ref_fasta ${schemeRepo}/nCoV-2019/V3/nCoV-2019.reference.fasta \
		--reads ${prefix}.fastq.gz \
		--outdir ${prefix}_outdir/ \
		--sample_name ${prefix} \
		--keep_bam
        cp ${prefix}_outdir/consensus.fa ${prefix}.fasta
        cp ${prefix}_outdir/log.json ${prefix}.viridian_log.json
        cp ${prefix}_outdir/variants.vcf ${prefix}.vcf
	cp ${prefix}_outdir/reference_mapped.bam ${prefix}.bam
        """
}
