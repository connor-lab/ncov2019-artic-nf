process readsTrimming {
    /**
    * Trims paired fastq using trim_galore (https://github.com/FelixKrueger/TrimGalore)
    * @input tuple(dataset_id, path(forward), path(reverse))
    * @output trimgalore_out tuple(dataset_id, path("*_val_1.fq.gz"), path("*_val_2.fq.gz"))
    */

    tag { dataset_id }

    publishDir "${params.output}/${task.process}", pattern: '*_val_{1,2}.fq.gz', mode: 'copy'

    cpus 2

    input:
    tuple(dataset_id, path(forward), path(reverse))

    output:
    tuple(dataset_id, path("*_val_1.fq.gz"), path("*_val_2.fq.gz")) optional true

    script:
    """
    if [[ \$(zcat ${forward} | head -n4 | wc -l) -eq 0 ]]; then
      exit 0
    else
      trim_galore --paired $forward $reverse
    fi
    """
}


