process illuminaDownloadScheme {
    tag params.scriptsRepoURL

    label 'internet'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "scripts", mode: "copy"

    output:
        path "scripts/${params.scriptsDir}/make_depth_mask.py" , emit: depthmask
        path "scripts/${params.scriptsDir}/vcftagprimersites.py" , emit: vcftagprimersites
        path "scripts/${params.scriptsDir}/*"

    script:
        """
        git clone ${params.scriptsRepoURL} scripts
        """
}

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
        bwa mem -t ${task.cpus} ${ref} ${forward} ${reverse} | samtools view -bS | \
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
        """
        samtools view -F4 -o ${sampleName}.mapped.bam ${bam}
        samtools index ${sampleName}.mapped.bam
        ${ivarCmd} -i ${sampleName}.mapped.bam -b ${bedfile} -m ${params.illuminaKeepLen} -q ${params.illuminaQualThreshold} -p ivar.out
        samtools sort -o ${sampleName}.mapped.primertrimmed.sorted.bam ivar.out.bam
        """
}

process callVariantsLofreq {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.lofreq.vcf", mode: 'copy'

    input:
        tuple(sampleName, path(bam), path(ref))

    output:
        tuple sampleName, path("${sampleName}.lofreq.vcf"), emit: variants

    script:
        """
        # Index bam file
        samtools index ${bam}

        lofreq indelqual --dindel ${bam} -f ${ref} |\
        lofreq call --call-indels --min-bq ${params.lofreqMinBaseQuality} --min-alt-bq ${params.lofreqMinBaseQuality} \
        --min-mq ${params.lofreqMinMapQuality} --no-default-filter --use-orphan --max-depth 1000000 \
        --min-cov ${params.lofreqMinCov} -f ${ref} -o ${sampleName}.lofreq.vcf -
        """
}

process findLowCoverageRegions {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.lowcoverage.txt", mode: 'copy'

    input:
        tuple(sampleName, path(bam), path(ref), path(depthmask), path(vcftagprimersites))

    output:
        tuple sampleName, path("${sampleName}.lowcoverage.txt"), emit: lowCoverageRegions

    script:
        """
        # was broken due to relative import- sed to patch this
        sed 's/from .vcftagprimersites import read_bed_file/from vcftagprimersites import read_bed_file/g' ${depthmask} > edited_depth_mask.py
        
        # Create bam index (required for pysam)
        samtools index ${bam}

        python edited_depth_mask.py --depth ${params.minDepthThreshold} ${ref} \
        ${bam} ${sampleName}.lowcoverage.txt
        """
}

process filterLowAlleleFrequencyVariants {
    
    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}*.vcf", mode: 'copy'

    input:
        tuple(sampleName, path(vcf))

    output:
        tuple sampleName, path("${sampleName}_filterlowaf"), emit: afFilteredVcf

    script:
        """
        # Switch off default filters
        lofreq filter --no-defaults --cov-min 20 --af-min 0.25 \
        --in ${vcf} \
        --out ${sampleName}_filterlowaf.vcf --print-all
        """
}

process splitPrimerSiteVariants {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}*.vcf", mode: 'copy'

    input:
        tuple(sampleName, path(vcf), path(schemeRepo))

    output:
        //tuple sampleName, path("${sampleName}.primersite.vcf"), emit: primerVariants
        //tuple sampleName, path("${sampleName}.notprimersite.vcf"), emit: nonPrimerVariants
        tuple sampleName, path("${sampleName}.primersite.vcf"), path("${sampleName}.notprimersite.vcf")

    script:
        """
        bedtools sort -i ${schemeRepo}/${params.schemeDir}/${params.scheme}/${params.schemeVersion}/nCoV-2019.scheme.bed |\
        bedtools merge > primers.merged.sorted.bed 
        bedtools intersect -header -v -a ${vcf} -b primers.merged.sorted.bed >${sampleName}.notprimersite.vcf
        bedtools intersect -header -a ${vcf} -b primers.merged.sorted.bed >${sampleName}.primersite.vcf
        """
}

process lofreqVariantFilters {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}*.vcf", mode: 'copy'

    input:
        tuple(sampleName, path(primerVcf), path(outsidePrimerVcf))

    output:
        tuple sampleName, path("${sampleName}.lofreqfiltered.vcf"), emit: lofreqFilteredVcf

    script:
        """
        # Default fdr strand bias and quality filter applied
        lofreq filter --sb-incl-indels --snvqual-mtc bonf --indelqual-mtc bonf \
        --in ${outsidePrimerVcf} \
        --out ${sampleName}.strandAndQualfiltered.vcf --print-all
        # Switch off default filters
        lofreq filter --no-defaults --snvqual-mtc bonf --indelqual-mtc bonf \
        --in ${primerVcf} \
        --out ${sampleName}.qualfiltered.vcf --print-all
        
        # Merge variants post filtering
        bgzip ${sampleName}.strandAndQualFiltered.vcf 
        bcftools index ${sampleName}.strandAndQualFiltered.vcf.gz
        bgzip ${sampleName}.qualFiltered.vcf 
        bcftools index ${sampleName}.qualFiltered.vcf.gz
        bcftools concat ${sampleName}.strandAndQualFiltered.vcf.gz ${sampleName}.qualFiltered.vcf.gz |\
        bcftools sort - > ${sampleName}.lofreqfiltered.vcf
        """
}

process customVariantFilters {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}*.vcf", mode: 'copy'

    input:
        tuple(sampleName, path(vcf))

    output:
        //tuple sampleName, path("${sampleName}.custompass.vcf"), emit: customPassVcf
        //tuple sampleName, path("${sampleName}.customfail.vcf"), emit: customFailVcf
        tuple sampleName, path("${sampleName}.customfiltered.vcf"), emit: customFilteredVcf

    script:
        """
        vcf_filter.py --illumina ${vcf} ${sampleName}.custompass.vcf ${dataset_id}.customfail.vcf

        # Merge variants post filtering
        bgzip ${sampleName}.custompass.vcf 
        bcftools index ${sampleName}.custompass.vcf.gz
        bgzip ${dataset_id}.customfail.vcf 
        bcftools index ${dataset_id}.customfail.vcf.gz
        bcftools concat ${sampleName}.custompass.vcf.gz ${dataset_id}.customfail.vcf.gz |\
        bcftools sort - > ${sampleName}.customfiltered.vcf
        """
}

process mergeCustomFilteredVcfs {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}*.vcf", mode: 'copy'

    input:
        tuple(sampleName, path(vcf))

    output:
        tuple sampleName, path("${sampleName}.lofreqfiltered.combined.vcf"), emit: mergedFilteredVcf

    script:
        """
        bgzip ${first_vcf}
        bgzip ${second_vcf}
        bcftools index ${first_vcf}.gz
        bcftools index ${second_vcf}.gz
        bcftools concat ${first_vcf}.gz ${second_vcf}.gz |\
        bcftools sort - > ${sampleName}.lofreqfiltered.combined.vcf
        """
}

process applyIupac {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}*.vcf", mode: 'copy'

    input:
        tuple(sampleName, path(vcf))

    output:
        tuple dataset_id, path("${sampleName}.edited.vcf"), emit: iupacVcf

    script:
        """
        vcf_edit.py -v ${vcf} -c lofreq -o ${sampleName}.edited.vcf
        """
}

process removeFilteredVariants {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}*.vcf", mode: 'copy'

    input:
        tuple(sampleName, path(vcf)) 

    output:
        tuple dataset_id, path("${sampleName}.pass.vcf"), emit: passVcf
        tuple dataset_id, path("${sampleName}.fail.vcf"), emit: failVcf


    script:
        """
        bgzip ${vcf}
        bcftools index ${vcf}.gz
        bcftools view -f .,PASS ${vcf}.gz -O v -o ${sampleName}.pass.vcf
        bcftools view -f min_af_0.250000,min_dp_20,sb_fdr,snvqual_bonf,indelqual_bonf,lowaf_indel ${vcf}.gz -O v -o ${sampleName}.fail.vcf
        """
}

process createConsensus {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}*.fasta", mode: 'copy'

    input:
        tuple(sampleName, path(passVcf), path(coverageMask), path(ref))

    output:
        tuple sampleName, path("${sampleName}.interim.consensus.fasta"), emit: interimConsensus

    script:
        """
        bgzip ${passVcf}
        bcftools index ${passVcf}.gz
        bcftools consensus -f ${ref} ${passVcf}.gz -m ${coverageMask} -o ${sampleName}.interim.consensus.fasta
        """
}

process addSampleNameToConsensusHeader {
    
    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}*.fasta", mode: 'copy'

    input:
        tuple(sampleName, path(fasta))

    output:
        tuple sampleName, path("${sampleName}.consensus.fasta"), emit: consensus

    script:
        """
        original_header=\$(grep ^">" ${fasta})
        sed "s/\$original_header/\$original_header|${sampleName}/g" ${fasta} > ${sampleName}.consensus.fasta
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

