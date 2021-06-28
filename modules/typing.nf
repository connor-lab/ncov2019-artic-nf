
process typeVariants {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/variants", pattern: "${sampleName}.variants.csv", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/vcf", pattern: "${sampleName}.csq.vcf", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/typing", pattern: "${sampleName}.typing.csv", mode: 'copy'

    input:
    tuple val(sampleName), path(variants), path(gff), path(ref), path(yaml)

    output:
    path "${sampleName}.variants.csv", optional: true, emit: variants_csv
    path "${sampleName}.typing.csv", optional: true, emit: typing_csv
    path "${sampleName}.csq.vcf", emit: csq_vcf

    script:
    if( params.illumina )
        """
        type_vcf.py \
        -i ${sampleName} \
        -y ${yaml} \
        -ov ${sampleName}.csq.vcf \
        -ot ${sampleName}.typing.csv \
        -os ${sampleName}.variants.csv \
        -dp ${params.csqDpThreshold} \
        -af ${params.csqAfThreshold} \
        -t ${variants} \
        ${gff} ${ref}
        """
    else
        """
        type_vcf.py \
        -i ${sampleName} \
        -y ${yaml} \
        -ov ${sampleName}.csq.vcf \
        -ot ${sampleName}.typing.csv \
        -os ${sampleName}.variants.csv \
        -dp ${params.csqDpThreshold} \
        -af ${params.csqAfThreshold} \
        -v ${variants} \
        ${gff} ${ref}
        """
}

process mergeTypingCSVs {

    tag { params.prefix }

    publishDir "${params.outdir}", pattern: "${params.prefix}.typing_summary.csv", mode: 'copy'
    publishDir "${params.outdir}", pattern: "${params.prefix}.variant_summary.csv", mode: 'copy'

    input:
    tuple path('typing/*'), path('variant/*')

    output:
    path "${params.prefix}.typing_summary.csv", emit: typing_summary_csv
    path "${params.prefix}.variant_summary.csv", emit: variant_summary_csv

    script:
    """
    #!/usr/bin/env python3
    import glob
    import csv

    dirs = ['typing', 'variant']

    for dir in dirs:
        globstring = dir + '/*.csv'
        files = glob.glob(globstring)

        header_written = False
        out_fn = "${params.prefix}." +dir+ '_summary.csv'
        with open(out_fn, 'w') as outfile:
            for fl in files:
                with open(fl, 'r' ) as csvfile:
                    csvreader = csv.DictReader(csvfile)
                    for row in csvreader:
                        if not header_written:
                            writer = csv.DictWriter(outfile, fieldnames=list(row.keys()))
                            writer.writeheader()
                            header_written = True

                        writer.writerow(row)

    """
}

