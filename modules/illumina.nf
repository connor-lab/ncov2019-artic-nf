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
    if [[ \$(zcat ${forward} | head -n4 | wc -l) -eq 0 ]]; then
      exit 0
    else
      trim_galore --paired $forward $reverse
    fi
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
        tuple(path(schemeRepo), sampleName, path(forward), path(reverse))

    output:
        tuple(sampleName, path("ref.fa"), path("${sampleName}.sorted.bam"))

    script:
        """
        ln -s ${schemeRepo}/${params.schemeDir}/${params.scheme}/${params.schemeVersion}/*.reference.fasta ref.fa
        bwa index ref.fa
        bwa mem -t ${task.cpus} ref.fa ${forward} ${reverse} | samtools view -bS | \
        samtools sort -o ${sampleName}.sorted.bam
        """
}

process makeIvarBedfile {

    tag { schemeRepo }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "ivar.bed", mode: 'copy'

    input:
    file(schemeRepo)

    output:
    file("ivar.bed")

    script:
    """
    #!/usr/bin/env python3
  
    import csv

    bedrows = []
    with open("${schemeRepo}/${params.schemeDir}/${params.scheme}/${params.schemeVersion}/nCoV-2019.scheme.bed", newline='') as bedfile:
        reader = csv.reader(bedfile, delimiter='\t')
        for row in reader:
            row[4] = '60'
            if row[3].endswith('LEFT'):
                 row.append('+')
            else: 
                row.append('-')
            bedrows.append(row)

    with open('ivar.bed', 'w', newline='') as bedfile:
        writer = csv.writer(bedfile, delimiter='\t')
        for row in bedrows:
            writer.writerow(row)
    """
}

process trimPrimerSequences {

    tag { sampleName }

    publishDir "${params.outdir}/climb_upload/${params.runPrefix}/${sampleName}", pattern: "${sampleName}.mapped.primertrimmed.sorted.bam", mode: 'copy'

    input:
        tuple(path(bedfile), sampleName, path(ref), path(bam))

    output:
        tuple(sampleName, path("${sampleName}.mapped.primertrimmed.sorted.bam"))

    script:
    if (params.allowNoprimer){
        ivarCmd = "ivar trim -e"
    } else {
        ivarCmd = "ivar trim"
    }
        """
        samtools view -F4 -o ivar.bam ${bam}
        samtools index ivar.bam
        ${ivarCmd} -i ivar.bam -b ${bedfile} ${params.illuminaKeepLen} -q ${params.illuminaQualThreshold} -p ivar.out
        samtools sort -o ${sampleName}.mapped.primertrimmed.sorted.bam ivar.out.bam
        """
}

process makeConsensus {

    tag { sampleName }

    publishDir "${params.outdir}/climb_upload/${params.runPrefix}/${sampleName}", pattern: "${sampleName}.primertrimmed.consensus.fa", mode: 'copy'

    input:
        tuple(sampleName, path(bam))

    output:
        tuple(sampleName, path("${sampleName}.primertrimmed.consensus.fa"))

    script:
        """
        samtools mpileup -A -d ${params.mpileupDepth} -Q0 ${bam} | \
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

    //publishDir "${params.outdir}/${task.process.replaceAll(":","_")}"

    input:
        file(cram)

    output:
        tuple(val("${cram.toString().replaceFirst(/\.cram/, "")}"), path("*_1.fastq.gz"),path("*_2.fastq.gz"))

    script:
        """
        samtools collate -u ${cram} -o tmp.bam
        samtools fastq -1 ${cram.toString().replaceFirst(/\.cram/, "_1.fastq.gz")} -2 ${cram.toString().replaceFirst(/\.cram/, "_2.fastq.gz")} tmp.bam
        rm tmp.bam
        """
}

