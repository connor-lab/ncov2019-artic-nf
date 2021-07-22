# NOTES
## INSTRUCTIONS TO DATA UPLOADER
1. Upload SAMPLE_FILES to the SAMPLE_DIRECTORY/
2. Requires - \
    i. For Illumina: illumina/<SAMPLE_NAME> \
    ii. For Nanopore \
        a. artic: artic/<SAMPLE_NAME> \
        b. midnight: midnight/<SAMPLE_NAME> 
3. The result directory will contains $SAMPLE_NAME_Results_$date \
e.g., SAMPLE_NAME = SRR001 \
The result directory name: SRR001_Results_2021-07-21
4. The SAMPLE files will be copied to another DIRECTORY 
5. The SAMPLE_FILES will be removed from the SAMPLE_DIRECTORY/


## Update Singularity containers
1. bi-weekly build new containers
2. check for updates.

## important info
1. to run pipeline after downloading the latest from git, need to change permission for files on bin/ folder.
