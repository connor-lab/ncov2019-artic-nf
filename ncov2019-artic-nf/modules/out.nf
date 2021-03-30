process bamToCram {
    /**
    * Convert bam to cram using samtools with embeded reference, indexes the CRAM (http://www.htslib.org/doc/samtools.html)
    * @input 
    * @output 
    */

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${bam.baseName}.cram", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${bam.baseName}.cram.crai", mode: 'copy'

    input:
        tuple(sampleName, path(bam), path(ref))

    output:
        tuple sampleName, path("${bam.baseName}.cram"), emit: cramed
        tuple sampleName, path("${bam.baseName}.cram.crai"), emit: cramedidx

    script:
        """
         samtools view --write-index -C -O cram,embed_ref -T ${ref} -o ${bam.baseName}.cram ${bam}
        """
}
