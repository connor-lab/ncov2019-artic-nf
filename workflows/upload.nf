
include {prepareUploadDirectory} from '../modules/upload.nf'
include {uploadToCLIMB} from '../modules/upload.nf'

workflow CLIMBrsync {
    take:
      ch_uploadDirectories
      ch_CLIMBkey

    main:
      prepareUploadDirectory(ch_uploadDirectories.collect())
      uploadToCLIMB(ch_CLIMBkey.combine(prepareUploadDirectory.out))
}

