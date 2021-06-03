process readTrimming {
    /**
    * Trims paired fastq using trim_galore (https://github.com/FelixKrueger/TrimGalore)
    * @input tuple(sampleName, path(forward), path(reverse))
    * @output trimgalore_out tuple(sampleName, path("*_val_1.fq.gz"), path("*_val_2.fq.gz"))
    */

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: '*_val_{1,2}.fq.gz', mode: 'copy'

    cpus 2

    input:
    tuple(sampleName, path(forward), path(reverse))

    output:
    tuple(sampleName, path("*_val_1.fq.gz"), path("*_val_2.fq.gz")) optional true

    script:
    """
    if [[ \$(gunzip -c ${forward} | head -n4 | wc -l) -eq 0 ]]; then
      exit 0
    else
      trim_galore --paired $forward $reverse
    fi
    """
}

process indexReference {
    /**
    * Indexes reference fasta file in the scheme repo using bwa.
    */

    tag { ref }

    input:
        path(ref)

    output:
        tuple path('ref.fa'), path('ref.fa.*')

    script:
        """
        ln -s ${ref} ref.fa
        bwa index ref.fa
        """
}

process readMapping {
    /**
    * Maps trimmed paired fastq using BWA (http://bio-bwa.sourceforge.net/)
    * Uses samtools to convert to BAM, sort and index sorted BAM (http://www.htslib.org/doc/samtools.html)
    * @input 
    * @output 
    */

    tag { sampleName }

    label 'largecpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.sorted.bam", mode: 'copy'

    input:
        tuple sampleName, path(forward), path(reverse), path(ref), path("*")

    output:
        tuple(sampleName, path("${sampleName}.sorted.bam"))

    script:
      """
      bwa mem -t ${task.cpus} ${ref} ${forward} ${reverse} | \
      samtools sort -o ${sampleName}.sorted.bam
      """
}

process trimPrimerSequences {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.mapped.bam", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.mapped.primertrimmed.sorted.bam", mode: 'copy'

    input:
    tuple sampleName, path(bam), path(bedfile)

    output:
    tuple sampleName, path("${sampleName}.mapped.bam"), emit: mapped
    tuple sampleName, path("${sampleName}.mapped.primertrimmed.sorted.bam" ), emit: ptrim

    script:
    if (params.allowNoprimer){
        ivarCmd = "ivar trim -e"
    } else {
        ivarCmd = "ivar trim"
    }
   
    if ( params.cleanBamHeader )
        """
        samtools reheader --no-PG  -c 'sed "s/${sampleName}/sample/g"' ${bam} | \
        samtools view -F4 -o sample.mapped.bam

        mv sample.mapped.bam ${sampleName}.mapped.bam
        
        samtools index ${sampleName}.mapped.bam

        ${ivarCmd} -i ${sampleName}.mapped.bam -b ${bedfile} -m ${params.illuminaKeepLen} -q ${params.illuminaQualThreshold} -p ivar.out

        samtools reheader --no-PG  -c 'sed "s/${sampleName}/sample/g"' ivar.out.bam | \
        samtools sort -o sample.mapped.primertrimmed.sorted.bam

        mv sample.mapped.primertrimmed.sorted.bam ${sampleName}.mapped.primertrimmed.sorted.bam
        """

    else
        """
        samtools view -F4 -o ${sampleName}.mapped.bam ${bam}
        samtools index ${sampleName}.mapped.bam
        ${ivarCmd} -i ${sampleName}.mapped.bam -b ${bedfile} -m ${params.illuminaKeepLen} -q ${params.illuminaQualThreshold} -p ivar.out
        samtools sort -o ${sampleName}.mapped.primertrimmed.sorted.bam ivar.out.bam
        """
}

process callVariants {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.variants.tsv", mode: 'copy'

    input:
    tuple(sampleName, path(bam), path(ref))

    output:
    tuple sampleName, path("${sampleName}.variants.tsv"), emit: variants

    script:
        """
        samtools mpileup -A -d 0 --reference ${ref} -B -Q 0 ${bam} |\
        ivar variants -r ${ref} -m ${params.ivarMinDepth} -p ${sampleName}.variants -q ${params.ivarMinVariantQuality} -t ${params.ivarMinFreqThreshold}
        """
}

process makeConsensus {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.primertrimmed.consensus.fa", mode: 'copy'

    input:
        tuple(sampleName, path(bam))

    output:
        tuple(sampleName, path("${sampleName}.primertrimmed.consensus.fa"))

    script:
        """
        samtools mpileup -aa -A -B -d ${params.mpileupDepth} -Q0 ${bam} | \
        ivar consensus -t ${params.ivarFreqThreshold} -m ${params.ivarMinDepth} \
        -n N -p ${sampleName}.primertrimmed.consensus
        """
}

process cramToFastq {
    /**
    * Converts CRAM to fastq (http://bio-bwa.sourceforge.net/)
    * Uses samtools to convert to CRAM, to FastQ (http://www.htslib.org/doc/samtools.html)
    * @input
    * @output
    */

    input:
        tuple sampleName, file(cram)

    output:
        tuple sampleName, path("${sampleName}_1.fastq.gz"), path("${sampleName}_2.fastq.gz")

    script:
        """
        samtools collate -u ${cram} -o tmp.bam
        samtools fastq -1 ${sampleName}_1.fastq.gz -2 ${sampleName}_2.fastq.gz tmp.bam
        rm tmp.bam
        """
}

process getObjFiles {
    /**
    * fetches fastq files from object store using OCI bulk download (https://docs.oracle.com/en-us/iaas/tools/oci-cli/2.24.4/oci_cli_docs/cmdref/os/object/bulk-download.html)
    * @input
    * @output
    */

    input:
        tuple bucket, prefix

    output:
        tuple prefix, path("${prefix}_C1.fastq.gz"), path("${prefix}_C2.fastq.gz")

    script:
	db=params.krakdb
        """
	oci os object bulk-download \
		-bn $bucket \
		--download-dir ./ \
		--overwrite \
		--auth instance_principal \
		--prefix $prefix

	kraken2 --paired -db ${db} \
		--report ${prefix}_summary.txt \
		--output ${prefix}_read_classification \
        	${prefix}_1.fastq.gz ${prefix}_2.fastq.gz

        awk '\$3==\"9606\" { print \$2 }' ${prefix}_read_classification >> kraken2_human_read_list
        awk '\$3!=\"9606\" { print \$2 }' ${prefix}_read_classification >> kraken2_nonhuman_read_list

	seqtk subseq ${prefix}_1.fastq.gz kraken2_nonhuman_read_list | gzip > "${prefix}_C1.fastq.gz"
	seqtk subseq ${prefix}_2.fastq.gz kraken2_nonhuman_read_list | gzip > "${prefix}_C2.fastq.gz"
	"""
}

process viridian {
    /**
    * runs viridian workflow https://github.com/iqbal-lab-org/viridian_workflow
    * @input
    * @output
    */

    publishDir "${params.outdir}/viridian"

    input:
        tuple prefix, path("${prefix}_1.fastq.gz"), path("${prefix}_2.fastq.gz"),path(bedfile), path('ref.fa'),path("*")

    output:
        tuple prefix, path("${prefix}_outdir/viridian/consensus.final_assembly.fa")

    script:
        """
        wget https://raw.githubusercontent.com/iqbal-lab-org/viridian_workflow/master/data/nCoV-artic-v3.bed
        viridian_workflow run_one_sample \
		ref.fa \
		nCoV-artic-v3.bed \
		${prefix}_1.fastq.gz \
		${prefix}_2.fastq.gz \
		${prefix}_outdir/
        """
}
