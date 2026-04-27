#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=120000]
#BSUB -q long
#BSUB -W 21:00
#BSUB -o "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210812_ChIPseq_Celseq2/%J%I.demultiplex.out"
#BSUB -e "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210812_ChIPseq_Celseq2/%J%I.demultiplex.err"
#########################################

## split library into samples in Celseq2 by CS2 barcode
module load java/1.8.0_171
module load condas/2018-05-11
source activate sshan_isoform
Dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210812_ChIPseq_Celseq2
for s in Celseq2_1_S1 Celseq2_2_S2 Celseq2_3_S3 Celseq2_4_S4 Celseq2_5_S5 Celseq2_6_S6;do
  echo "working on" ${s}
  date
  cd $Dir/basespace/fastqs/laneMerged
  barcodeMap=$Dir/data_info_Celseq2/barcodeMap/${s}
  gunzip ${s}_R1.fastq.gz && gunzip ${s}_R2.fastq.gz
  java -Xmx100g -jar /home/ed70w/bin/splitter_09.21.15_13.57.jar F1=${s}_R1.fastq F2=${s}_R2.fastq B=NNNNNNSSSSSS M=$barcodeMap HD=1 O=${s}
  gzip ${s}_R1.fastq && gzip ${s}_R2.fastq
  cd $Dir/basespace/fastqs/laneMerged/${s}
  for i in $(ls *fq);do
    fname=$(echo ${i} | sed 's/.fq//g')
    cat ${i} | awk '{if (NR%4==1) gsub("_",":",$1); print}' | gzip > $Dir/basespace/fastqs/laneMerged/${fname}.fq.gz
    rm ${i}
    echo ${i}
    date
  done
  echo ${s}
  date
done

