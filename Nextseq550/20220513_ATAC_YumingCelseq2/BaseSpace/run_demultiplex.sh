#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=120000]
#BSUB -q long
#BSUB -W 21:00
#BSUB -o "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20220513_ATAC_YumingCelseq2/BaseSpace/%J%I.demultiplex.out"
#BSUB -e "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20220513_ATAC_YumingCelseq2/BaseSpace/%J%I.demultiplex.err"
#########################################

## split library into samples in Celseq2 by CS2 barcode
module load java/1.8.0_171
module load condas/2018-05-11
source activate fastQTL
Dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20220513_ATAC_YumingCelseq2/BaseSpace

for s in Celseq2_pool1 Celseq2_pool2 Celseq2_pool3 Celseq2_pool4 Celseq2_pool5 Celseq2_pool6 Celseq2_pool7; do
  echo "working on" ${s}; date

  cd $Dir/fastqs_Celseq2/skin_eQTL_CelSeq2
  barcodeMap=$Dir/data_info_Celseq2/barcodeMap/${s}
  gunzip ${s}_R1.fastq.gz && gunzip ${s}_R2.fastq.gz
  java -Xmx100g -jar /home/ed70w/bin/splitter_09.21.15_13.57.jar F1=${s}_R1.fastq F2=${s}_R2.fastq B=NNNNNNSSSSSS M=$barcodeMap HD=1 O=${s}
  gzip ${s}_R1.fastq && gzip ${s}_R2.fastq

  cd $Dir/fastqs_Celseq2/skin_eQTL_CelSeq2/${s}
  for i in $(ls *fq);do
    fname=$(echo ${i} | sed 's/.fq//g')
    cat ${i} | awk '{if (NR%4==1) gsub("_",":",$1); print}' | gzip > $Dir/fastqs_Celseq2/skin_eQTL_CelSeq2/${fname}.fq.gz
    rm ${i}
    echo ${i}; date
  done
  echo ${s}; date
done

