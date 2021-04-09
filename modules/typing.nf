process alignSeqs {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/msa", mode: 'copy', overwrite: false, pattern: "${sampleName}.muscle.aln"
 
    tag { sampleName }

    input:
      tuple sampleName, path(sample), path(reference)

    output:
      tuple sampleName, path("${sampleName}.muscle.aln")

    script:
      """
      sed "s/>.*/>${sampleName}/g" $sample > ${sampleName}.clean.fa
      cat $reference ${sampleName}.clean.fa > pre.aln
      muscle -in pre.aln -out ${sampleName}.muscle.aln
      """
}

process typeVariants {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/typing_json", mode: 'copy', overwrite: true, pattern: "${sampleName}.json.gz"
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}/variant_csv", mode: 'copy', overwrite: true, pattern: "${sampleName}.csv"
 
    tag { sampleName }

    input:
      tuple sampleName, refName, path(msa), path(yaml_dir), path(gb)

    output:
      path("aln2type.${sampleName}.csv"), emit: typing_csv optional true
      path("${sampleName}.csv"), emit: variants_csv optional true
      path("${sampleName}.json.gz") optional true

    script:
      if ( gb.getBaseName() != 'dummyfile' ){
        """
        aln2type --gb ${gb} --output_unclassified . . aln2type.${sampleName}.csv ${refName} ${msa} ${yaml_dir}/*.yml
        """
      } else {
        """
        aln2type --output_unclassified . . aln2type.${sampleName}.csv ${refName} ${msa} ${yaml_dir}/*.yml
        """
      }
      
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

