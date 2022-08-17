process readTrimming {
    /**
    * Trims paired fastq using trim_galore (https://github.com/FelixKrueger/TrimGalore)
    * @input tuple(sampleName, path(forward), path(reverse))
    * @output trimgalore_out tuple(sampleName, path("*_val_1.fq.gz"), path("*_val_2.fq.gz"))
    */

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: '*_val_{1,2}.fq.gz', mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: '*_trimming_report.txt', mode: 'copy'

    cpus 2

    input:
    tuple(sampleName, path(forward), path(reverse))

    output:
    path "*_trimming_report.txt", optional: true, emit: log
    tuple sampleName, path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), optional: true, emit: trim

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
        tuple sampleName, path("${sampleName}.sorted.bam")

    script:
      """
        bwa mem -t ${task.cpus} ${ref} ${forward} ${reverse} | \
        samtools sort -o ${sampleName}.sorted.bam
      """
}

process flagStat {
    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.flagstat", mode: 'copy'

    input:
    //tuple sampleName, path(forward), path(reverse), path(ref), path("*")
    tuple sampleName, path(bam)

    output:
    tuple sampleName, path("${sampleName}.flagstat")

    script:
    """
    samtools flagstat ${sampleName}.sorted.bam  >${sampleName}.flagstat
    """
}

process trimPrimerSequences {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.mapped.bam", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.mapped.primertrimmed.sorted.bam", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.mapped.primertrimmed.sorted.bam.bai", mode: 'copy'

    input:
    tuple sampleName, path(bam), path(bedfile)


    output:
    tuple sampleName, path("${sampleName}.mapped.bam"), emit: mapped
    tuple sampleName, path("${sampleName}.mapped.primertrimmed.sorted.bam" ), emit: ptrim
    tuple sampleName, path("${sampleName}.mapped.primertrimmed.sorted.bam.bai"), emit: bai

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
        samtools index sample.mapped.primertrimmed.sorted.bam

        mv sample.mapped.primertrimmed.sorted.bam ${sampleName}.mapped.primertrimmed.sorted.bam
        mv sample.mapped.primertrimmed.sorted.bam.bai ${sampleName}.mapped.primertrimmed.sorted.bam.bai
        """

    else
        """
        samtools view -F4 -o ${sampleName}.mapped.bam ${bam}
        samtools index ${sampleName}.mapped.bam
        ${ivarCmd} -i ${sampleName}.mapped.bam -b ${bedfile} -m ${params.illuminaKeepLen} -q ${params.illuminaQualThreshold} -p ivar.out
        samtools sort -o ${sampleName}.mapped.primertrimmed.sorted.bam ivar.out.bam
        samtools index ${sampleName}.mapped.primertrimmed.sorted.bam
        """
}

process depth {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.depth", mode: 'copy'

    input:
    tuple(sampleName, path(bam), path(bai))

    output:
    tuple sampleName, path("${sampleName}.depth"), emit: depth

    script:
        """
        samtools index ${sampleName}.mapped.primertrimmed.sorted.bam
        sambamba depth base -c0 ${sampleName}.mapped.primertrimmed.sorted.bam -o ${sampleName}.depth
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
        samtools faidx ${ref}
        samtools mpileup -A -d 0 --reference ${ref} -B -Q 0 ${bam} |\
        ivar variants -r ${ref} -m ${params.varMinDepth} -p ${sampleName}.variants -q ${params.varMinVariantQuality} -t ${params.varMinFreqThreshold}
        """
}

process makeConsensus {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.primertrimmed.consensus.fa", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.primertrimmed.consensus.qual.txt", mode: 'copy'

    input:
        tuple(sampleName, path(bam))

    output:
        tuple sampleName, path("${sampleName}.primertrimmed.consensus.qual.txt")
        tuple sampleName, path("${sampleName}.primertrimmed.consensus.fa"), emit: consensus_fasta

    script:
        """
        samtools mpileup -aa -A -B -d ${params.mpileupDepth} -Q0 ${bam} | \
        ivar consensus -t ${params.varFreqThreshold} -m ${params.varMinDepth} \
        -n N -p ${sampleName}.primertrimmed.consensus
        """
}

process freebayes {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.variants.norm.vcf", mode: 'copy'

    input:
    tuple(sampleName, path(bam), path(ref))

    output:
    tuple sampleName, path("${sampleName}.variants.norm.vcf"), emit:vcf

    script:
        """
        freebayes -p 1 \
                  --min-coverage ${params.freeMinCov} \
                  --min-base-quality ${params.freeMinBaseQual} \
                  -f ${ref} \
                  -F ${params.freeMinAltFrac} \
                  -m ${params.freeMinMapQual} ${bam} > ${sampleName}.freebayes.raw.vcf
        bcftools norm ${sampleName}.freebayes.raw.vcf \
                  -f ${ref} \
                  -o ${sampleName}.variants.norm.vcf
        """
}

process annotationVEP {
    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.freebayes.vep.vcf", mode: 'copy'

    input:
        tuple sampleName, path(vcf),path(ref)

    output:
        tuple sampleName, path("${sampleName}.freebayes.vep.vcf")

    script:
        """
        bgzip -i -f -c ${baseDir}/typing/MN908947.3.gff >MN908947.3.gff.gz
        tabix -f MN908947.3.gff.gz

        bgzip -i -f -c ${sampleName}.variants.norm.vcf > ${sampleName}.variants.norm.vcf.gz

        vep -i ${sampleName}.variants.norm.vcf.gz --format vcf --gff MN908947.3.gff.gz --fasta ${ref} -o ${sampleName}.freebayes.vep.vcf --vcf --force_overwrite --no_stats --distance 10 --hgvs
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
