#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=120000]
#BSUB -q long
#BSUB -W 21:00
#BSUB -o "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20211208_Celseq2_ATACseq/log/%J%I.demultiplex.out"
#BSUB -e "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20211208_Celseq2_ATACseq/log/%J%I.demultiplex.err"

########################################
# split library into samples in Celseq2 by CS2 barcode
module load java/1.8.0_171
module load condas/2018-05-11
source activate sshan_isoform
Dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20211208_Celseq2_ATACseq
for s in F28K_F30K F30MP_F50MP_F62MP_F62MI F33K_F34K F38MP_F38MI_F39MP_F42MP F53F_F55F_F56F_F57F F58F_F59F_F60F_F61F F62F_F63F;do
  echo "working on" ${s}
  date
  cd /nl/umw_manuel_garber/human/skin/eQTLs/Nextseq2000/20211208_Celseq2_ATACseq/fastq
  barcodeMap=$Dir/data_info_Celseq2/barcodeMap/${s}
  gunzip ${s}_R1.fastq.gz && gunzip ${s}_R2.fastq.gz
  java -Xmx100g -jar /home/ed70w/bin/splitter_09.21.15_13.57.jar F1=${s}_R1.fastq F2=${s}_R2.fastq B=NNNNNNSSSSSS M=$barcodeMap HD=1 O=${s}
  gzip ${s}_R1.fastq && gzip ${s}_R2.fastq
  rm ${s}_R1.fastq ${s}_R2.fastq
  cd $Dir/fastq/${s}
  for i in $(ls *fq);do
    fname=$(echo ${i} | sed 's/.fq//g')
    cat ${i} | awk '{if (NR%4==1) gsub("_",":",$1); print}' | gzip > $Dir/fastq/Celseq2/${fname}.fq.gz
    rm ${i}
    echo ${i}
    date
  done
  echo ${s}
  date
done
