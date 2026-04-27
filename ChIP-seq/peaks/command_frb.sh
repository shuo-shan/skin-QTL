#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=40800]
#BSUB -q long
#BSUB -W 121:00
#BSUB -e "./%J%I.err"
#BSUB -o "./%J%I.out"

module load samtools/1.16.1
module load macs2/2.2.7.1

cd /pi/manuel.garber-umw/human/skin/eQTLs/ChIP-seq/peaks
samtools merge mergedBamFile.bam *.bam
samtools merge mergedBamFile_FRB.bam ChIP_F25_FRB_IFN.bam ChIP_F25_FRB_PBS.bam ChIP_F49_FRB_IFN.bam ChIP_F49_FRB_PBS.bam ChIP_F55_FRB_IFN.bam ChIP_F55_FRB_PBS.bam
echo "done merging all bam files"; date
macs2 callpeak --bw 300 -t mergedBamFile.bam -n merged_all_FRB_skin-eQTL_H3K27ac_ChIPseq_files -g hs
echo "done calling summit and narrowpeaks for all merged bam files";date



