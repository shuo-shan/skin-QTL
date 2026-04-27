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

cd /pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/peaks
samtools merge mergedBamFile.bam *.bam
echo "done merging all bam files"; date
macs2 callpeak --bw 300 -t mergedBamFile.bam -n merged_all_skin-eQTL_ATACseq_files -g hs
echo "done calling summit and narrowpeaks for all merged bam files";date

