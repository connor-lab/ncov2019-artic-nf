process pango {
    tag { sampleName }

    publishDir "${params.outdir}/analysis/pango/${params.prefix}", mode: 'copy'

    input:
    tuple(sampleName,  path(fasta))

    output:
    tuple(sampleName, path("${sampleName}_lineage_report.csv"))

    script:
    """
    pangolin ${fasta}
    mv lineage_report.csv ${sampleName}_lineage_report.csv
    """
}


process nextclade {
    tag { sampleName }

    publishDir "${params.outdir}/analysis/nextclade/${params.prefix}", mode: 'copy'

    input:
    tuple(sampleName,  path(fasta))

    output:
    tuple(sampleName, path("${sampleName}.tsv"))

    script:
    """
    nextclade --input-fasta ${fasta} \
        --output-tree ${sampleName}_tree.json \
        --output-tsv ${sampleName}.tsv \
        --output-json ${sampleName}.json
    """
}

process getVariantDefinitions {
    output:
    path('variant_definitions') 

    script:
    """
    git clone https://github.com/phe-genomics/variant_definitions
    """
}



process aln2type {
    tag { sampleName }

    publishDir "${params.outdir}/analysis/aln2type/${params.prefix}", mode: 'copy'

    input:
    tuple(sampleName,  path(fasta),path(variant_definitions), path(reffasta), path("*"))

    output:
    tuple(sampleName, path("${sampleName}.csv")) optional true

    script:
    """
    cat $reffasta  ${fasta} > unaligned.fasta
    mafft --auto unaligned.fasta > aln.fasta
    aln2type sample_json_out \
	sample_csv_out \
	--output_unclassified \
	${sampleName}.csv \
	MN908947.3 \
	aln.fasta \
	variant_definitions/variant_yaml/*.yml

    """
}


process makeReport {
    tag { sampleName }

    publishDir "${params.outdir}/analysis/report/${params.prefix}", mode: 'copy'

    input:
    tuple(sampleName, path('pango.csv'), path('aln2type.csv'), path('nextclade.tsv'))

    output:
    path("${sampleName}_report.tsv"), emit: tsv
    path("${sampleName}_report.json"), emit: json

    script:
    """
    makeReport.py ${sampleName}
    """
}

process FN4_upload {
    tag { sampleName }

    input:
    tuple(sampleName,  path(fasta),path(variant_definitions), path(reffasta), path("*"))

    script:
    """
    mafft --auto \
        --thread 1 \
        --addfull ${fasta} \
        --keeplength $reffasta \
        > ${sampleName}_wuhan.fa

    FN4ormater.py -i ${sampleName}_wuhan.fa -r MN908947.3 -s ${sampleName} -o ${sampleName}.fasta

    oci os object put \
	-bn FN4-queue \
	--force \
        --auth instance_principal \
	--file ${sampleName}.fasta \
	--metadata "{\\"sampleID\\":\\"$sampleName\\"}"

    """
}
