process getObjCsv {
    /**
    * fetches CSV file (sp3data.csv) from object store using OCI bulk download (https://docs.oracle.com/en-us/iaas/tools/oci-cli/2.24.4/oci_cli_docs/cmdref/os/object/bulk-download.html)
    * @input
    * @output
    */

    //tag { prefix }

    input:
        tuple bucket, path

    output:
        path("*.csv")

    script:
    """
    oci os object get \
        -bn $bucket \
        --auth instance_principal \
        --file sp3data.csv \
        --name $path
    """
}