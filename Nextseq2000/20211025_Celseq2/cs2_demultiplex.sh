#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=120000]
#BSUB -q long
#BSUB -W 21:00
#BSUB -o "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20211025_Celseq2/%J%I.splitfq_by_cs2.out"
#BSUB -e "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20211025_Celseq2/%J%I.splitfq_by_cs2.err"

########################################
# split library into samples in Celseq2 by CS2 barcode
module load java/1.8.0_171
module load condas/2018-05-11
source activate sshan_isoform
Dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20211025_Celseq2
for s in F36F_F37F_F38F_F39F_S8 F40F_F41F_F42F_F44F_S9 F45F_F46F_F47F_F48F_S10 F46K_F47K_S2 F48K_F49K_F23F_F24F_S3 F49F_F50F_F51F_F52F_S11 F50K_F51K_F27F_F28F_S4 F52K_F53K_F30F_F31F_S5 F55K_F56K_F32F_F33F_S6 F57K_F58K_F34F_F35F_S7 F59K_F60K_S1;do
  echo "working on" ${s}
  date
  cd /nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20211025_Celseq2/fastq
  barcodeMap=$Dir/data_info_Celseq2/barcodeMap/${s}
  gunzip ${s}_L001_R1_001.fastq.gz && gunzip ${s}_L001_R2_001.fastq.gz
  java -Xmx100g -jar /home/ed70w/bin/splitter_09.21.15_13.57.jar F1=${s}_L001_R1_001.fastq F2=${s}_L001_R2_001.fastq B=NNNNNNSSSSSS M=$barcodeMap HD=1 O=${s}
  gzip ${s}_L001_R1_001.fastq && gzip ${s}_L001_R2_001.fastq
  rm ${s}_L001_R1_001.fastq ${s}_L001_R2_001.fastq
  cd $Dir/fastq/${s}
  for i in $(ls *fq);do
    gzip ${i}
    echo ${i}
    date
  done
  echo ${s}
  date
done
