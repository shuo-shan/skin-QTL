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
dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements
cd ${dir}

genome_sizes=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes

##### find midpoint for each region and expand by 100 on both sides
cat allcts_promoters.bed | awk '{OFS="\t"}{if (($3-$2)%2==1) print $1,($3-$2+1)/2,($3-$2+1)/2+1,$4; else print $1,($3-$2)/2,($3-$2)/2+1,$4}' > allcts_promoters_summit.bed
cat allcts_enhancers.bed | awk '{OFS="\t"}{if (($3-$2)%2==1) print $1,($3-$2+1)/2,($3-$2+1)/2+1,$4; else print $1,($3-$2)/2,($3-$2)/2+1,$4}' > allcts_enhancers_summit.bed

cat allcts_promoters_summit.bed | bedtools sort -i stdin | bedtools slop -i stdin -g ${genome_sizes} -b 100 | bedtools sort -i stdin > allcts_promoters_200bp.bed
cat allcts_enhancers_summit.bed | bedtools sort -i stdin | bedtools slop -i stdin -g ${genome_sizes} -b 100 | bedtools sort -i stdin > allcts_enhancers_200bp.bed

##### calculate coverage in each 200bp region around summit
while read c; do
	echo ${c} | bsub -W 8:00 -n 1 -R "span[hosts=1]" -R rusage[mem=1020] -q long -o “./%J%I.out” -e “./%J%I.err”
done < commands_multicov.txt

##### merge all coverage tables together
rm header.txt
echo -e "chr\nstart\nend\nname" > header.txt
for f in coverage_ATAC_*_enhancers_200bp.bed;do	sample=$(echo ${f} | cut -d'_' -f3,4,5); echo ${sample} >> header.txt; done

cat coverage_ATAC_F25_FRB_IFN_merged_allcts_enhancers_200bp.bed | cut -f1,2,3,4 > rows.txt
for f in coverage_ATAC_*_enhancers_200bp.bed;do cat ${f} | cut -f5 > tempcov_${f}; echo "done with "${f}; done
paste rows.txt tempcov_*.bed > merged_coverage_ATAC_enhancers_200bp.bed
rm tempcov_* rows.txt

cat coverage_ATAC_F25_FRB_IFN_merged_allcts_promoters_200bp.bed | cut -f1,2,3,4 > rows.txt
for f in coverage_ATAC_*_promoters_200bp.bed;do cat ${f} | cut -f5 > tempcov_${f}; echo "done with "${f}; done
paste rows.txt tempcov_*.bed > merged_coverage_ATAC_promoters_200bp.bed
rm tempcov_* rows.txt

cat coverage_ATAC_F25_FRB_IFN_merged_allcts_enhancers.bed | cut -f1,2,3,4 > rows.txt
for f in coverage_ATAC_*_enhancers.bed;do cat ${f} | cut -f5 > tempcov_${f}; echo "done with "${f}; done
paste rows.txt tempcov_*.bed > merged_coverage_ATAC_enhancers.bed
rm tempcov_* rows.txt

cat coverage_ATAC_F25_FRB_IFN_merged_allcts_promoters.bed | cut -f1,2,3,4 > rows.txt
for f in coverage_ATAC_*_promoters.bed;do cat ${f} | cut -f5 > tempcov_${f}; echo "done with "${f}; done
paste rows.txt tempcov_*.bed > merged_coverage_ATAC_promoters.bed
rm tempcov_* rows.txt



