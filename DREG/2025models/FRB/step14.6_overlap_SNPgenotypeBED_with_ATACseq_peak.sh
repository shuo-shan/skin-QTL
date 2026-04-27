#!/bin/bash

this_bed=$1
dir=$(dirname ${this_bed})

module load bedtools
head -1 /pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks/ATACseq_peak_annotation.tsv > ${dir}/temp_output
awk 'NR>1' /pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks/ATACseq_peak_annotation.tsv > /pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks/ATACseq_peak_annotation.bed

bedtools intersect \
	-a ${this_bed} \
	-b /pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks/ATACseq_peak_annotation.bed \
	-wa -wb \
	>> ${dir}/temp_output
mv ${dir}/temp_output ${this_bed}.new.bed
