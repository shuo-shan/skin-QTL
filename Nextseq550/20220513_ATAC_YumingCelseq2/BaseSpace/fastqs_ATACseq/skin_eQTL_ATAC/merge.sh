#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=204000]
#BSUB -q long
#BSUB -W 121:00
#BSUB -e "./%J%I.err"
#BSUB -o "./%J%I.out"

for file in ATAC_F25_KRT_IFN_24h_100ng ATAC_F49_KRT_IFN_24h_100ng ATAC_F49_KRT_PBS_24h_0ng ATAC_F55_KRT_PBS_24h_0ng;do
	zcat ${file}*R1_001.fastq.gz > ${file}_R1.fastq; gzip ${file}_R1.fastq; zcat ${file}*R2_001.fastq.gz > ${file}_R2.fastq; gzip ${file}_R2.fastq;
	echo done with ${file}
done
