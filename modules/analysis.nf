process getVariantDefinitions {
    
    output:
    path('variant_definitions') 

    script:
    """
    git clone https://github.com/phe-genomics/variant_definitions
    """
}

process makeReport {
    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: 'copy', pattern: "${sampleName}_report.tsv"

    input:
    tuple(sampleName, path('pangolinTyping.csv'), path('nextclade_tree.json'), path('nextclade.tsv'),
		path('nextclade.json'))

    output:
    path "${sampleName}_report.tsv", emit: tsv

    script:
    """
    makeReport.py ${sampleName}
    """
}
<<<<<<< Updated upstream
=======

process versions {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "*", mode: 'copy'

    output:
    file "*versions.csv"

    script:
    if ( params.illumina )
        """
        bcftools -v > version_bcftools.txt
        bwa > version_bwa.txt 2>&1 || true
        fastqc -v > version_fastqc.txt
        gofasta -v > version_gofasta.txt
        freebayes --version > version_freebayes.txt
        ivar version > version_ivar.txt
        multiqc --version > version_multiqc.txt
        nextclade --version > version_nextclade.txt
        picard CollectWgsMetrics -version > version_picard.txt 2>&1 || true
        python --version > version_python.txt
        samtools --version > version_samtools.txt
        snakemake -v > version_snakemake.txt
        trim_galore -v > version_trim_galore.txt
        get_versions.py ${params.prefix}_versions.csv
        """
    else if ( params.medaka || params.nanopolish )
        """
        artic --version > version_artic.txt
        fastqc -v > version_fastqc.txt
        multiqc --version > version_multiqc.txt
        bwa > version_bwa.txt 2>&1 || true
        mafft --version > version_mafft.txt 2>&1
        medaka --version > version_medaka.txt
        minimap2 --version > version_minimap2.txt
        muscle -version > version_muscle.txt
        nanopolish --version > version_nanopolish.txt
        nextclade --version > version_nextclade.txt
        porechop --version > version_porechop.txt
        python --version > version_python.txt
        samtools --version > version_samtools.txt
<<<<<<< Updated upstream
        snakemake -v > version_snakemake.txt
        get_versions.py ${params.prefix}_versions.csv
        """
}


process pangoversions {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "*", mode: 'copy'

    output:
    file "*versions.csv"

    script:
    if ( params.illumina )
        """
        scorpio -cv > version_constellations.txt
        gofasta -v > version_gofasta.txt
        pangolin -pv > version_pangolin-data.txt
        pangolin -v > version_pangolin.txt
        python --version > version_python.txt
        scorpio -v > version_scorpio.txt
        snakemake -v > version_snakemake.txt
        usher --version > version_usher.txt
        get_versions.py ${params.prefix}_versions.csv
        """
    else if ( params.medaka || params.nanopolish )
        """
        scorpio -cv > version_constellations.txt
        gofasta -v > version_gofasta.txt
        pangolin -pv > version_pangolin-data.txt
        pangolin -v > version_pangolin.txt
        python --version > version_python.txt
        scorpio -v > version_scorpio.txt
=======
>>>>>>> Stashed changes
        snakemake -v > version_snakemake.txt
        get_versions.py ${params.prefix}_versions.csv
        """
}
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======


process pangoversions {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "*", mode: 'copy'

    output:
    file "*versions.csv"

    script:
    if ( params.illumina )
        """
        scorpio -cv > version_constellations.txt
        gofasta -v > version_gofasta.txt
        pangolin -pv > version_pangolin-data.txt
        pangolin -v > version_pangolin.txt
        python --version > version_python.txt
        scorpio -v > version_scorpio.txt
        snakemake -v > version_snakemake.txt
        usher --version > version_usher.txt
        get_versions.py ${params.prefix}_versions.csv
        """
    else if ( params.medaka || params.nanopolish )
        """
        scorpio -cv > version_constellations.txt
        gofasta -v > version_gofasta.txt
        pangolin -pv > version_pangolin-data.txt
        pangolin -v > version_pangolin.txt
        python --version > version_python.txt
        scorpio -v > version_scorpio.txt
        snakemake -v > version_snakemake.txt
        get_versions.py ${params.prefix}_versions.csv
        """
}
>>>>>>> Stashed changes
