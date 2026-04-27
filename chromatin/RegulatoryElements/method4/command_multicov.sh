#!/bin/bash
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=2040]
#BSUB -q long
#BSUB -W 121:00
#BSUB -e "./%J%I.err"
#BSUB -o "./%J%I.out"

# script overview: define promoters and enhancers based on H3K27ac peaks
# method: ELisz 2018 paper. For details see Crystal Dropbox log on promoter_enhancer_annotation.docx

module load bedtools/2.30.0
dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3
cd ${dir}

genome_sizes=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
atac_window=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/peaks/merged_all_skin-eQTL_ATACseq_files_summits_300bp_flanking_window.bed
dict=${dir}/dictionary_atac_peaks_annotated_H3K27ac_activity.txt

# 1. create a bed file of filtered ATACseq peaks for bedtools
cat ${dict} | cut -f1 > temp.peaks 
bash /pi/manuel.garber-umw/sshan/scripts/function_rapid_fgrep.sh ${dir}/temp.peaks ${atac_window} this.filtered.peaks hehe 5000
cat this.filtered.peaks | cut -f1-4 > this.filtered.peaks.bed

##### calculate coverage in each 200bp region around summit
while read c; do
	echo ${c} | bsub -W 8:00 -n 1 -R "span[hosts=1]" -R rusage[mem=2040] -q long -o “./log/multicov_%J%I.out” -e “./log/multicov_%J%I.err”
done < commands_multicov_ATAC.txt

while read c; do
	echo ${c} | bsub -W 8:00 -n 1 -R "span[hosts=1]" -R rusage[mem=1020] -q long -o “./log/multicov_%J%I.out” -e “./log/multicov_%J%I.err”
done < commands_multicov_ChIP.txt

##### merge all coverage tables together
mkdir coverage
mv coverage*.bed coverage
cd coverage
# first, ATACseq coverage table
rm header.txt
echo -e "chr\nstart\nend\nname" > header.txt
for f in coverage_ATAC*.bed; do
	sample=$(echo ${f} | cut -d'_' -f3,4,5)
	echo ${sample} >> header.txt
done 
cat coverage_ATAC_F25_FRB_IFN_merged_allcts_regulatory_region.bed | cut -f1,2,3,4 > rows.txt
for f in coverage_ATAC*.bed; do
	cat ${f} | cut -f5 > tempcov_${f}
	echo "done with "${f}
done
cat header.txt | tr '\n' '\t' | sed 's/$/\n/g' > merged_coverage_ATACseq_allcts_regulatory_region.bed
paste rows.txt tempcov_*.bed >> merged_coverage_ATACseq_allcts_regulatory_region.bed
rm header.txt rows.txt tempcov_*

# next, ChIPseq coverage table
rm header.txt
echo -e "chr\nstart\nend\nname" > header.txt
for f in coverage_ChIP*.bed; do
	sample=$(echo ${f} | cut -d'_' -f3,4,5)
	echo ${sample} >> header.txt
done
cat coverage_ChIP_F25_FRB_IFN_merged_allcts_regulatory_region.bed | cut -f1,2,3,4 > rows.txt
for f in coverage_ChIP*.bed; do
	cat ${f} | cut -f5 > tempcov_${f}
	echo "done with "${f}
done
cat header.txt | tr '\n' '\t' | sed 's/$/\n/g' > merged_coverage_ChIPseq_allcts_regulatory_region.bed
paste rows.txt tempcov_*.bed >> merged_coverage_ChIPseq_allcts_regulatory_region.bed
rm header.txt rows.txt tempcov_* 

##### move the merged tables to main directory
mv merged_coverage_ATACseq_allcts_regulatory_region.bed ../
mv merged_coverage_ChIPseq_allcts_regulatory_region.bed ../



