import java.nio.file.Paths

def makeFastqSearchPath ( illuminaPrefixes, illuminaSuffixes, fastq_exts ) {

    def fastqSearchPath = []

    for (suff in illuminaSuffixes){
        for(ext in fastq_exts){
          if ( illuminaPrefixes ) {
            for (prefix in illuminaPrefixes) {
              
              // Make a glob to recurse directories
              dirNameGlob = params.directory.replaceAll(/\/+$/, "") + '**'
 
              // Make a filename glob
              fileNameGlob = prefix + suff + ext

              // Build a path
              searchPath = Paths.get(dirNameGlob, fileNameGlob )

              fastqSearchPath.add(searchPath.toString())
              }
          } else {

              // Make a glob to recurse directories
              dirNameGlob = params.directory.replaceAll(/\/+$/, "") + '**'

              // Make a glob for filenames
              fileNameGlob = suff + ext

              // Build a path
              searchPath = Paths.get(dirNameGlob, fileNameGlob)

              fastqSearchPath.add(searchPath.toString())
          }
        }
    }

    return fastqSearchPath
}

def getCsvBucketPath ( filePath  ) {
    objPath = java.nio.file.Paths.get("/data/inputs/s3")

    try {
        csvPath = objPath.relativize(java.nio.file.Paths.get(filePath))
        bucket = csvPath.getName(0)
        path = bucket.relativize(csvPath)
        return [ bucket, path ]

    } catch (Exception e) {
        println("Path to CSV with --catsup does not exists on this system and couldn't find a potential bucket path. Quitting.")
        println(e)
        System.exit(1)
    }
}
