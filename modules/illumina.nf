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
    tuple(val(sampleName), path(forward), path(reverse))

    output:
    tuple(val(sampleName), path("*_val_1.fq.gz"), path("*_val_2.fq.gz")) optional true

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
        tuple(path('ref.fa'), path('ref.fa.*'))

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
        tuple(val(sampleName), path(forward), path(reverse), path(ref), path("*"))

    output:
        tuple(val(sampleName), path("${sampleName}.sorted.bam"))

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
    tuple(val(sampleName), path(bam), path(bedfile))

    output:
    tuple(val(sampleName), path("${sampleName}.mapped.bam"), emit: mapped)
    tuple(val(sampleName), path("${sampleName}.mapped.primertrimmed.sorted.bam" ), emit: ptrim)

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

process getDepths {
    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.depths.tsv", mode: 'copy'

    input:
    tuple(val(sampleName), path(bam), path(ref))

    output:
    tuple(val(sampleName), path("${sampleName}.depths.tsv"), emit: depths)

    script:
        """
        samtools mpileup -aa -A -d 0 -Q 0 -q ${params.ivarMinVariantQuality} -B -f ${ref} ${bam} | cut -f1-4 > "${sampleName}.depths.tsv"
        """
}
process callVariants {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.variants.tsv", mode: 'copy'

    errorStrategy { sleep(Math.pow(2, task.attempt) * 1 as long); return 'retry' }
    maxRetries 9999999999

    input:
    tuple(val(sampleName), path(bam), path(ref))

    output:
    tuple(val(sampleName), path("${sampleName}.variants.tsv"), emit: variants)

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
        tuple(val(sampleName), path(bam))

    output:
        tuple(val(sampleName), path("${sampleName}.primertrimmed.consensus.fa"), emit: consensus)

    script:
        """
        samtools mpileup -aa -A -B -d ${params.mpileupDepth} -Q0 ${bam} | \
        ivar consensus -t ${params.ivarFreqThreshold} -m ${params.ivarMinDepth} \
        -n N -p ${sampleName}.primertrimmed.consensus
        """
}

process callLineage {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.pangolin.csv", mode: 'copy'

    errorStrategy { sleep(Math.pow(2, task.attempt) * 1 as long); return 'retry' }
    maxRetries -1

    input:
        tuple(val(sampleName), path(consensus))

    output:
        tuple(val(sampleName), path("${sampleName}.pangolin.csv"))

    script:
        """
        pangolin ${consensus} -t 1 --outfile ${sampleName}.pangolin.csv --verbose
        """
}

process freyjaDemix {
    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.freyja.demix.tsv", mode: 'copy'

    errorStrategy 'ignore'

    input:
        tuple(val(sampleName), path(variants), path(depths))

    output:
        tuple(val(sampleName), path("${sampleName}.freyja.demix.tsv"))

    script:
        """
        freyja demix ${variants} ${depths} --output ${sampleName}.freyja.demix.tsv
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
        tuple(val(sampleName), file(cram))

    output:
        tuple(val(sampleName), path("${sampleName}_1.fastq.gz"), path("${sampleName}_2.fastq.gz"))

    script:
        """
        samtools collate -u ${cram} -o tmp.bam
        samtools fastq -1 ${sampleName}_1.fastq.gz -2 ${sampleName}_2.fastq.gz tmp.bam
        rm tmp.bam
        """
}

