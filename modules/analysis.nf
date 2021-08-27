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


process download_nextclade_files {
    publishDir "${params.outdir}/analysis/nextclade/${params.prefix}", mode: 'copy'
    output:
    file('nextclade_files')

    script:
    """
    nextclade_ver=`(nextclade -v)`
    wget -P nextclade_files https://raw.githubusercontent.com/nextstrain/nextclade/\${nextclade_ver}/data/sars-cov-2/tree.json
    wget -P nextclade_files https://raw.githubusercontent.com/nextstrain/nextclade/\${nextclade_ver}/data/sars-cov-2/genemap.gff
    wget -P nextclade_files https://raw.githubusercontent.com/nextstrain/nextclade/\${nextclade_ver}/data/sars-cov-2/qc.json 
    echo \$nextclade_ver > nextclade_files/version.txt
    """

}


process nextclade {
    tag { sampleName }

    publishDir "${params.outdir}/analysis/nextclade/${params.prefix}", mode: 'copy'

    input:
    tuple(sampleName,  path(fasta), path(reffasta), path(nextclade_files)) 

    output:
    tuple sampleName, path("${sampleName}.tsv"), emit: tsv
    tuple sampleName, path("${sampleName}.json"), emit: json

    script:
    """
    nextclade --input-fasta ${fasta} \
	--input-root-seq ${reffasta} \
	--input-gene-map=nextclade_files/genemap.gff \
	--input-tree=nextclade_files/tree.json \
	--input-qc-config=nextclade_files/qc.json \
	--output-json=${sampleName}.json \
	--output-tsv=${sampleName}.tsv \
	--output-basename=${sampleName}
    nextclade_ver=`(nextclade -v)`
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
    tuple(sampleName,  file('consensus.fasta'), path(variant_definitions), path(reffasta), path("*"))

    script:
    """
    mafft --auto \
        --thread 1 \
        --addfull 'consensus.fasta' \
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
