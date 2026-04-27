#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=10000]
#BSUB -q long
#BSUB -W 8:00
#BSUB -e "./step3.%J%I.err"
#BSUB -o "./step3.%J%I.out"
# script from https://raw.githubusercontent.com/wiki/GuttmanLab/sprite-pipeline/3.-Filtering.md
# remove reads that fall in the repeatmasker regions and hg38 blacklist regions

module load bedtools/2.29.2
module load samtools/1.9

prefix=SPRITE-F37-KRT-PBS
nope_regions=/nl/umw_manuel_garber/human/skin/eQTLs/SPRITE/pipeline/sprite-pipeline-master/hg38_blacklist_rmsk.milliDivLessThan140.bed
dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq500/20220506_SPRITE_SMARTseq
cd ${dir}

bedtools intersect -v -a ${dir}/${prefix}.DNA.bowtie2.mapq20.bam -b ${nope_regions} > ${dir}/${prefix}.DNA.bowtie2.mapq20.masked.bam
samtools sort ${dir}/${prefix}.DNA.bowtie2.mapq20.masked.bam -o ${dir}/${prefix}.DNA.bowtie2.mapq20.masked.sorted.bam
samtools index ${dir}/${prefix}.DNA.bowtie2.mapq20.masked.sorted.bam
