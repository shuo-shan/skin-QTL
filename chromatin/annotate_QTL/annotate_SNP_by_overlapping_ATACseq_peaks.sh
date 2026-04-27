#!/bin/bash
module load bedtools

workingdir=$1
SNPbed=$2
prefix=$3

########## fetch SNPs that overlap ATACseq peaks
genome=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
atacPeaks=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks/ATACseq_peak_annotation.tsv
cd ${workingdir}
awk 'NR>1' ${atacPeaks} | bedtools sort -g ${genome} > ${workingdir}/ATACseq_peak_annotation.bed

awk 'NR==1' ${atacPeaks} | awk '{OFS=FS="\t"}{print "snp_chr","snp_start","snp_end","snp","REF","ALT",$0}' > ${workingdir}/${prefix}_SNP_ATAC_overlap.txt
bedtools intersect -a ${SNPbed} -b ${workingdir}/ATACseq_peak_annotation.bed -wa -wb >> ${workingdir}/${prefix}_SNP_ATAC_overlap.txt


echo "Wrote result to ${workingdir}/${prefix}_SNP_ATAC_overlap.txt"; date
