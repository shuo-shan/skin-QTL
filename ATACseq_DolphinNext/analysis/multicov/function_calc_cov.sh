#!/bin/bash
# coverage of each sample at merged MACS2 called peaks
# shuo.shan@umassmed.edu
# Nov2021

#### user input
sample=$1

#### packages & env
module load bedtools/2.29.2
module load samtools/1.9

wd=/nl/umw_manuel_garber/human/skin/eQTLs/ATACseq_DolphinNext/analysis/multicov
bamDir=/nl/umw_manuel_garber/human/skin/eQTLs/ATACseq_DolphinNext/softlinks/bam
sampleF=${bamDir}/${sample}.bam
bedF=/nl/umw_manuel_garber/human/skin/eQTLs/ATACseq_DolphinNext/analysis/merged.bed
#### only run if bam hasn't been indexed 
#cd ${bamDir}
#samtools index ${sampleF}

#### calculate coverage
cd ${wd}
bedtools multicov -bams ${sampleF} -bed ${bedF} > ${sample}.merged.bg
date

