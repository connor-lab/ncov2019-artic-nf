process viridianPrimers {
    /**
    * runs viridian workflow https://github.com/iqbal-lab-org/viridian_workflow
    * @input
    * @output
    */

    tag { prefix }

    publishDir "${params.outdir}/consensus_seqs/", mode: 'copy', pattern: "*.fasta", saveAs: { filename -> "${prefix}.fasta"}
    publishDir "${params.outdir}/VCF/", mode: 'copy', pattern: "*.vcf", saveAs: { filename -> "${prefix}.vcf"}
    publishDir "${params.outdir}/qc/", mode: 'copy', pattern: "*.json", saveAs: { filename -> "${prefix}.viridian_log.json"}
    if (params.TESToutputMODE){
        publishDir "${params.outdir}/bam/", mode: 'copy', pattern: "*.bam", saveAs: { filename -> "${prefix}.bam"}
    }

    input:
        tuple prefix, path("${prefix}_1.fastq.gz"), path("${prefix}_2.fastq.gz"),path('primers'), path('ref.fa'),path("*")

    output:
        tuple prefix, path("${prefix}_outdir/consensus.fa"), emit: consensus
        tuple prefix, path("${prefix}_outdir/log.json"), emit: coverage
        tuple prefix, path("${prefix}_outdir/variants.vcf"), emit: vcfs
        tuple prefix, path("${prefix}_outdir/reference_mapped.bam"), emit: bam

 
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
    """ 
}


process viridianAuto {
    /** 
    * runs viridian workflow https://github.com/iqbal-lab-org/viridian_workflow
    * @input
    * @output
    */

    tag { prefix }

    publishDir "${params.outdir}/consensus_seqs/", mode: 'copy', pattern: "*.fasta", saveAs: { filename -> "${prefix}.fasta"}
    publishDir "${params.outdir}/VCF/", mode: 'copy', pattern: "*.vcf", saveAs: { filename -> "${prefix}.vcf"}
    publishDir "${params.outdir}/qc/", mode: 'copy', pattern: "*.json", saveAs: { filename -> "${prefix}.viridian_log.json"}
    if (params.TESToutputMODE){
        publishDir "${params.outdir}/bam/", mode: 'copy', pattern: "*.bam", saveAs: { filename -> "${prefix}.bam"}
    }

    input:
        tuple prefix, path("${prefix}_1.fastq.gz"), path("${prefix}_2.fastq.gz"), path('ref.fa'),path("*")


    output:
        tuple prefix, path("${prefix}_outdir/consensus.fa"), emit: consensus
        tuple prefix, path("${prefix}_outdir/log.json"), emit: coverage
        tuple prefix, path("${prefix}_outdir/variants.vcf"), emit: vcfs
        tuple prefix, path("${prefix}_outdir/reference_mapped.bam"), emit: bam

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
    """ 
}


process viridianONTPrimers {
    /**
    * runs viridian workflow https://github.com/iqbal-lab-org/viridian_workflow
    * @input
    * @output
    */

    tag { prefix }
    publishDir "${params.outdir}/consensus_seqs/", mode: 'copy', pattern: "*.fasta", saveAs: { filename -> "${prefix}.fasta"}
    publishDir "${params.outdir}/VCF/", mode: 'copy', pattern: "*.vcf", saveAs: { filename -> "${prefix}.vcf"}
    publishDir "${params.outdir}/qc/", mode: 'copy', pattern: "*.json", saveAs: { filename -> "${prefix}.viridian_log.json"}
    if (params.TESToutputMODE){
        publishDir "${params.outdir}/bam/", mode: 'copy', pattern: "*.bam", saveAs: { filename -> "${prefix}.bam"}
    }

    output:
        tuple prefix, path("${prefix}_outdir/consensus.fa"), emit: consensus
        tuple prefix, path("${prefix}_outdir/log.json"), emit: coverage
        tuple prefix, path("${prefix}_outdir/variants.vcf"), emit: vcfs
        tuple prefix, path("${prefix}_outdir/reference_mapped.bam"), emit: bam

    input:
        tuple prefix, path("${prefix}.fastq.gz"),path(schemeRepo),path('primers')


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
        """
}

process viridianONTAuto {
    /**
    * runs viridian workflow https://github.com/iqbal-lab-org/viridian_workflow
    * @input
    * @output
    */

    tag { prefix }
    publishDir "${params.outdir}/consensus_seqs/", mode: 'copy', pattern: "*.fasta", saveAs: { filename -> "${prefix}.fasta"}
    publishDir "${params.outdir}/VCF/", mode: 'copy', pattern: "*.vcf", saveAs: { filename -> "${prefix}.vcf"}
    publishDir "${params.outdir}/qc/", mode: 'copy', pattern: "*.json", saveAs: { filename -> "${prefix}.viridian_log.json"}
    if (params.TESToutputMODE){
        publishDir "${params.outdir}/bam/", mode: 'copy', pattern: "*.bam", saveAs: { filename -> "${prefix}.bam"}
    }


    input:
        tuple prefix, path("${prefix}.fastq.gz"),path(schemeRepo)

    output:
        tuple prefix, path("${prefix}_outdir/consensus.fa"), emit: consensus
        tuple prefix, path("${prefix}_outdir/log.json"), emit: coverage
        tuple prefix, path("${prefix}_outdir/variants.vcf"), emit: vcfs
        tuple prefix, path("${prefix}_outdir/reference_mapped.bam"), emit: bam

    script:
        """
        viridian_workflow run_one_sample \
		--tech ont \
		--ref_fasta ${schemeRepo}/nCoV-2019/V3/nCoV-2019.reference.fasta \
		--reads ${prefix}.fastq.gz \
		--outdir ${prefix}_outdir/ \
		--sample_name ${prefix} \
		--keep_bam
        """
}
