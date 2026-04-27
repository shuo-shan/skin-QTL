#!/bin/bash
# written by shuo.shan@umassmed.edu, 03/2024

conda activate fastQTL
module load bedtools
dir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/QTL_region_TF_enrichment
cd ${dir}

# compile region of interest
region_type=promoter
region=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method4/promoters_and_enhancers_surrounding_genes_all3cts.bed
cat ${region} | grep promoter | grep -w FRB > region.bed
# resources for celltype
#f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all_with_colnames.txt
#bed=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/snps_near_expressed_genes.bed
#f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/KRT_minimal/masteroutput_all_with_colnames.txt
#bed=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/KRT_minimal/snps_near_expressed_genes.bed
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/FRB_new/masteroutput_all_with_colnames.txt
bed=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/FRB_new/snps_near_expressed_genes.bed

# reQTL
QTL_type=FRB_reQTL
cat ${f} | awk '{OFS=FS="\t"}NR>1{if ($7 != "." && $7 < 0.000001) print $1}' | sort | uniq > QTL.txt # reQTL
grep -w -f QTL.txt ${bed} > QTL.bed
# list of non-QTLs
cat ${f} | grep -v -w -f QTL.txt | cut -f1 | sort | uniq > nonQTL.txt
grep -w -f nonQTL.txt ${bed} > nonQTL.bed
# overlap QTLs with regions
bedtools intersect -a <(bedtools sort -i QTL.bed) -b <(bedtools sort -i region.bed) | cut -f4 | sort | uniq > QTL_in_region.txt
# QTLs NOT in regions
bedtools intersect -v -a <(bedtools sort -i QTL.bed) -b <(bedtools sort -i region.bed) | cut -f4 | sort | uniq > QTL_not_in_region.txt
# overlap nonQTLs with regions
bedtools intersect -a <(bedtools sort -i nonQTL.bed) -b <(bedtools sort -i region.bed) | cut -f4 | sort | uniq > nonQTL_in_region.txt
# nonQTLs NOT in regions
bedtools intersect -v -a <(bedtools sort -i nonQTL.bed) -b <(bedtools sort -i region.bed) | cut -f4 | sort | uniq > nonQTL_not_in_region.txt
# calculate OR and 95% CI
QTL_in_region=$(wc -l QTL_in_region.txt | cut -d' ' -f1)
QTL_not_in_region=$(wc -l QTL_not_in_region.txt | cut -d' ' -f1)
nonQTL_in_region=$(wc -l nonQTL_in_region.txt | cut -d' ' -f1)
nonQTL_not_in_region=$(wc -l nonQTL_not_in_region.txt | cut -d' ' -f1)
Rscript calculate_OddsRatio_and_95CI.R ${QTL_in_region} ${QTL_not_in_region} ${nonQTL_in_region} ${nonQTL_not_in_region} ${QTL_type} ${region_type} ${dir}/ODDS_RATIO_95CI.txt

# PBSeQTL
QTL_type=FRB_PBSeQTL
cat ${f} | awk '{OFS=FS="\t"}NR>1{if ($15 != "." && $15 < 0.000001) print $1}' | sort | uniq > QTL.txt # PBSeQTL
grep -w -f QTL.txt ${bed} > QTL.bed
cat ${f} | grep -v -w -f QTL.txt | cut -f1 | sort | uniq > nonQTL.txt
grep -w -f nonQTL.txt ${bed} > nonQTL.bed
bedtools intersect -a <(bedtools sort -i QTL.bed) -b <(bedtools sort -i region.bed) | cut -f4 | sort | uniq > QTL_in_region.txt
bedtools intersect -v -a <(bedtools sort -i QTL.bed) -b <(bedtools sort -i region.bed) | cut -f4 | sort | uniq > QTL_not_in_region.txt
bedtools intersect -a <(bedtools sort -i nonQTL.bed) -b <(bedtools sort -i region.bed) | cut -f4 | sort | uniq > nonQTL_in_region.txt
bedtools intersect -v -a <(bedtools sort -i nonQTL.bed) -b <(bedtools sort -i region.bed) | cut -f4 | sort | uniq > nonQTL_not_in_region.txt
# calculate OR and 95% CI
QTL_in_region=$(wc -l QTL_in_region.txt | cut -d' ' -f1)
QTL_not_in_region=$(wc -l QTL_not_in_region.txt | cut -d' ' -f1)
nonQTL_in_region=$(wc -l nonQTL_in_region.txt | cut -d' ' -f1)
nonQTL_not_in_region=$(wc -l nonQTL_not_in_region.txt | cut -d' ' -f1)
Rscript calculate_OddsRatio_and_95CI.R ${QTL_in_region} ${QTL_not_in_region} ${nonQTL_in_region} ${nonQTL_not_in_region} ${QTL_type} ${region_type} ${dir}/ODDS_RATIO_95CI.txt

# IFNeQTL
QTL_type=FRB_IFNeQTL
cat ${f} | awk '{OFS=FS="\t"}NR>1{if ($23 != "." && $23 < 0.000001) print $1}' | sort | uniq > QTL.txt # IFNeQTL
grep -w -f QTL.txt ${bed} > QTL.bed
cat ${f} | grep -v -w -f QTL.txt | cut -f1 | sort | uniq > nonQTL.txt
grep -w -f nonQTL.txt ${bed} > nonQTL.bed
bedtools intersect -a <(bedtools sort -i QTL.bed) -b <(bedtools sort -i region.bed) | cut -f4 | sort | uniq > QTL_in_region.txt
bedtools intersect -v -a <(bedtools sort -i QTL.bed) -b <(bedtools sort -i region.bed) | cut -f4 | sort | uniq > QTL_not_in_region.txt
bedtools intersect -a <(bedtools sort -i nonQTL.bed) -b <(bedtools sort -i region.bed) | cut -f4 | sort | uniq > nonQTL_in_region.txt
bedtools intersect -v -a <(bedtools sort -i nonQTL.bed) -b <(bedtools sort -i region.bed) | cut -f4 | sort | uniq > nonQTL_not_in_region.txt
# calculate OR and 95% CI
QTL_in_region=$(wc -l QTL_in_region.txt | cut -d' ' -f1)
QTL_not_in_region=$(wc -l QTL_not_in_region.txt | cut -d' ' -f1)
nonQTL_in_region=$(wc -l nonQTL_in_region.txt | cut -d' ' -f1)
nonQTL_not_in_region=$(wc -l nonQTL_not_in_region.txt | cut -d' ' -f1)
Rscript calculate_OddsRatio_and_95CI.R ${QTL_in_region} ${QTL_not_in_region} ${nonQTL_in_region} ${nonQTL_not_in_region} ${QTL_type} ${region_type} ${dir}/ODDS_RATIO_95CI.txt

# clean-up
rm QTL.txt QTL.bed nonQTL.txt nonQTL.bed *QTL*region.txt
rm region.bed


