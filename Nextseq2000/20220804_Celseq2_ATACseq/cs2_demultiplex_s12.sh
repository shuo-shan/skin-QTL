#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=120000]
#BSUB -q long
#BSUB -W 21:00
#BSUB -o "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20220804_Celseq2_ATACseq/log/%J%I.demultiplex.out"
#BSUB -e "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20220804_Celseq2_ATACseq/log/%J%I.demultiplex.err"

########################################
# split library into samples in Celseq2 by CS2 barcode
module load java/1.8.0_171
Dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20220804_Celseq2_ATACseq
for s in RNA_12_S12;do
  echo "working on" ${s}; date
  cd /nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20220804_Celseq2_ATACseq/fastqs_RNAseq/20220804_Celseq2_ATACseq/
  barcodeMap=$Dir/data_info_Celseq2/barcodeMap/${s}
  zcat ${s}_L001_R1_001.fastq.gz ${s}_L002_R1_001.fastq.gz > ${s}_R1.fastq;
  zcat ${s}_L001_R2_001.fastq.gz ${s}_L002_R2_001.fastq.gz > ${s}_R2.fastq 
  #gunzip ${s}_R1.fastq.gz && gunzip ${s}_R2.fastq.gz
  java -Xmx100g -jar /home/ed70w/bin/splitter_09.21.15_13.57.jar F1=${s}_R1.fastq F2=${s}_R2.fastq B=NNNNNNSSSSSS M=$barcodeMap HD=1 O=${s}
  echo "done with splitter program";date;
  gzip ${s}_R1.fastq && gzip ${s}_R2.fastq
  #rm ${s}_R1.fastq ${s}_R2.fastq
  echo "done with gzip";date;
  cd $Dir/fastqs_RNAseq/20220804_Celseq2_ATACseq/${s}
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
