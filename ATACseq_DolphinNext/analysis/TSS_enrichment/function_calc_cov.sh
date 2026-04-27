#!/bin/bash
# TSS enrichment score
# shuo.shan@umassmed.edu
# Nov2021

#### user input
sample=$1
#### packages & env
module load bedtools/2.29.2
module load samtools/1.9

wd=/nl/umw_manuel_garber/human/skin/eQTLs/ATACseq_DolphinNext/analysis/TSS_enrichment
bamDir=/nl/umw_manuel_garber/human/skin/eQTLs/ATACseq_DolphinNext/softlinks/bam
sampleF=${bamDir}/${sample}.bam

#### 
cd ${bamDir}
samtools index ${sampleF}
cd ${wd}
bedtools multicov -bams ${sampleF} -bed TSS.binned.bed > TSS.binned.cov.${sample}.bg
date

