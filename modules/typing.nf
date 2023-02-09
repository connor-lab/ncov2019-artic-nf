
process typeVariants {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/variants", pattern: "${sampleName}.variants.csv", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/vcf", pattern: "${sampleName}.csq.vcf", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/typing", pattern: "${sampleName}.typing.csv", mode: 'copy'

    input:
    tuple sampleName, path(variants), path(gff), path(ref), path(yaml)

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

process pangolinTyping {
    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: 'copy', pattern: "${sampleName}.pangolin.csv"


    input:
    tuple val(sampleName), path(consensus_fasta)

    output:
    tuple val(sampleName), path("${sampleName}.pangolin.csv"), emit: pangolin

    script:
    """
    pangolin ${consensus_fasta} --max-ambig 0.2 --outfile ${sampleName}.pangolin.csv
    """
}

process nextclade {
    tag { sampleName }

    label 'nextclade'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: 'copy', pattern: "${sampleName}.tsv"

    input:
    tuple val(sampleName), path(consensus_fasta)

    output:
    tuple(sampleName, path("${sampleName}_tree.json"),
    path("${sampleName}.tsv"),path("${sampleName}.json"))

    script:
    """
    echo \$(nextclade --version 2>&1) > nextclade_version.txt
    nextclade dataset get --name 'sars-cov-2' --output-dir 'data/sars-cov-2'
    nextclade run \
        --input-dataset data/sars-cov-2 \
        --output-tree ${sampleName}_tree.json \
        --output-tsv ${sampleName}.tsv \
        --output-json ${sampleName}.json \
        ${consensus_fasta}

    """

}
