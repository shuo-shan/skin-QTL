#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=204000]
#BSUB -q long
#BSUB -W 121:00
#BSUB -e "./%J%I.err"
#BSUB -o "./%J%I.out"

for file in Celseq2_pool1 Celseq2_pool2 Celseq2_pool3 Celseq2_pool4 Celseq2_pool5 Celseq2_pool6 Celseq2_pool7;do
	zcat ${file}*R1_001.fastq.gz > ${file}_R1.fastq; gzip ${file}_R1.fastq; zcat ${file}*R2_001.fastq.gz > ${file}_R2.fastq; gzip ${file}_R2.fastq;
	echo done with ${file}
done
