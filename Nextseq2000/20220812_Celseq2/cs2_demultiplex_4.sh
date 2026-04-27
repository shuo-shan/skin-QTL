#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=120000]
#BSUB -q long
#BSUB -W 21:00
#BSUB -o "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20220812_Celseq2/log/%J%I.demultiplex.out"
#BSUB -e "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20220812_Celseq2/log/%J%I.demultiplex.err"

########################################
# split library into samples in Celseq2 by CS2 barcode
module load java/1.8.0_171
Dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20220812_Celseq2
for s in CS_RNA_4_S4;do
#for s in CS_RNA_1_S1  CS_RNA_2_S2  CS_RNA_3_S3  CS_RNA_4_S4  CS_RNA_5_S5  CS_RNA_6_S6  CS_RNA_7_S7  CS_RNA_8_S8;do
  echo "working on" ${s}; date
  cd /nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20220812_Celseq2/fastq_RNAseq/2022_08_09_skin-eQTL_RNA
  barcodeMap=$Dir/data_info_Celseq2/barcodeMap/${s}
  #zcat ${s}_L001_R1_001.fastq.gz ${s}_L002_R1_001.fastq.gz > ${s}_R1.fastq;
  #zcat ${s}_L001_R2_001.fastq.gz ${s}_L002_R2_001.fastq.gz > ${s}_R2.fastq 
  gunzip ${s}_R1.fastq.gz && gunzip ${s}_R2.fastq.gz
  java -Xmx100g -jar /home/ed70w/bin/splitter_09.21.15_13.57.jar F1=${s}_R1.fastq F2=${s}_R2.fastq B=NNNNNNSSSSSS M=$barcodeMap HD=1 O=${s}
  echo "done with splitter program";date;
  gzip ${s}_R1.fastq && gzip ${s}_R2.fastq
  #rm ${s}_R1.fastq ${s}_R2.fastq
  echo "done with gzip";date;
  cd $Dir/fastq_RNAseq/2022_08_09_skin-eQTL_RNA/${s}
  for i in $(ls *fq);do
    fname=$(echo ${i} | sed 's/.fq//g')
    cat ${i} | awk '{if (NR%4==1) gsub("_",":",$1); print}' | gzip > $Dir/fastq_barcodeModified/${fname}.fq.gz
    rm ${i}
    echo ${i}
    date
  done
  echo ${s}
  date
done
