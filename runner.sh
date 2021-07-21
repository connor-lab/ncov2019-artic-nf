#!/bin/bash
module load nextflow-20.10.0
module load singularity-3.7.1

INPUT_DIRS=data/
OUT_DIRS=$(cd $INPUT_DIRS; find * -type d -prune -exec ls -d {} \; |head -1)
echo $OUT_DIRS

ILLUM_DIRS=data/illumina/
ILLUM_RES=$(cd $ILLUM_DIRS; ls -td * |head -1)

ARTIC_DIRS=data/artic/
ARTIC_RES=$(cd $ARTIC_DIRS; find * -type d -prune -exec ls -d {} \; |head -1)

MIDNIGHT_DIRS=data/midnight/
MIDNIGHT_RES=$(cd $MIDNIGHT_DIRS; find * -type d -prune -exec ls -d {} \; |head -1)

if [ $OUT_DIRS == "illumina" ]; then
echo "
     ========================================== 
     FOUND ILLUMINA SAMPLE DIRECTORY: $ILLUM_RES  
     ==========================================="
echo "
      ====================================== 
      START ANALYSIS WITH ILLUMINA PIPELINE 
      ======================================"
#nextflow run main.nf -profile singularity \
#    --illumina --prefix "${ILLUM_RES}_illumina"     \
#    --directory data/   \
#    --outdir ${ILLUM_RES}_Results_`date +%Y-%m-%d`
wait
echo "
     ====================================
       ILLUMINA ANALYSIS COMPLETED!!!
     ===================================="
echo "
      =======================================
      DELETING SAMPLE DIRECTORY $ILLUM_RES 
      ======================================="

rm -rf $PWD/data/illumina/$ILLUM_RES

elif [ $OUT_DIRS == "artic" ]; then
echo "
     =================================================
     FOUND NANOPORE ARTIC SAMPLE DIRECTORY: $ARTIC_RES 
     ================================================="
echo "
      ===========================================
      START ANALYSIS WITH NANOPORE ARTIC PIPELINE 
      ==========================================="
#nextflow run main.nf -profile singularity \
#    --nanopolish --prefix "${ARTIC_RES}_nanopore" \
#    --basecalled_fastq data/artic/${ARTIC_RES}/fastq_pass/ \
#    --fast5_pass data/artic/${ARTIC_RES}/fast5_pass/ \
#    --sequencing_summary data/artic/${ARTIC_RES}/*.txt \
#    --outdir ${ARTIC_RES}_Results_`date +%Y-%m-%d`
wait
echo "
     ======================================
       NANOPORE ARTIC ANALYSIS COMPLETED!!!  
     ======================================"
echo "
      ====================================
      DELETING SAMPLE DIRECTORY $ARTIC_RES 
      ===================================="

rm -rf $PWD/data/artic/$ARTIC_RES

else [ $OUT_DIRS == "midnight" ]
echo "
     ==============================================
     FOUND NANOPORE SAMPLE DIRECTORY: $MIDNIGHT_RES 
     =============================================="
echo "
      ==============================================
      START ANALYSIS WITH NANOPORE MIDNIGHT PIPELINE 
      =============================================="
#nextflow run main.nf -profile singularity \
#    --nanopolish --prefix "${MIDNIGHT_RES}_nanopore" \
#    --basecalled_fastq data/midnight/${MIDNIGHT_RES}/fastq_pass/ \
#    --fast5_pass data/midnight/${MIDNIGHT_RES}/fast5_pass/ \
#    --sequencing_summary data/midnight/${MIDNIGHT_RES}/*.txt \
#    --scheme-directory primer_schemes/midnight/nCoV-2019/V1/ \
#    --outdir ${MIDNIGHT_RES}_Results_`date +%Y-%m-%d`
wait
echo "
     =========================================
       NANOPORE MIDNIGHT ANALYSIS COMPLETED!!!  
     ========================================="
echo "
      =======================================
      DELETING SAMPLE DIRECTORY $MIDNIGHT_RES 
      ======================================="

rm -rf $PWD/data/midnight/$MIDNIGHT_RES
fi