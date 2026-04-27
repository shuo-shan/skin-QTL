#!/bin/bash
# written by shuo.shan@umassmed.edu, 03/2024
source activate fastQTL
module load bedtools
dir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/QTL_region_TF_enrichment

# region of interest
TF=$1 # from /pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/QTL_region_TF_enrichment/TF_list.txt
region_type=${TF}_peak
mkdir -p ${dir}/${TF}

cd /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENCODE/humanTFChIP/
less filtered_metadata.txt | grep ${TF} | cut -d';' -f3 > temp_id_${TF}.txt
ls | grep -f temp_id_${TF}.txt  > temp_files_${TF}.txt
while read f;do
	zcat ${f} >> all_files_${TF}.bed
done < temp_files_${TF}.txt
bedtools sort -i all_files_${TF}.bed > all_files_sorted_${TF}.bed
bedtools merge -i all_files_sorted_${TF}.bed > ${dir}/${TF}/region.bed
rm temp_id_${TF}.txt temp_files_${TF}.txt all_files_sorted_${TF}.bed
cd ${dir}/${TF}

######################### FRB
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
Rscript ${dir}/calculate_OddsRatio_and_95CI.R ${QTL_in_region} ${QTL_not_in_region} ${nonQTL_in_region} ${nonQTL_not_in_region} ${QTL_type} ${region_type} ${dir}/ODDS_RATIO_95CI.txt

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
Rscript ${dir}/calculate_OddsRatio_and_95CI.R ${QTL_in_region} ${QTL_not_in_region} ${nonQTL_in_region} ${nonQTL_not_in_region} ${QTL_type} ${region_type} ${dir}/ODDS_RATIO_95CI.txt

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
Rscript ${dir}/calculate_OddsRatio_and_95CI.R ${QTL_in_region} ${QTL_not_in_region} ${nonQTL_in_region} ${nonQTL_not_in_region} ${QTL_type} ${region_type} ${dir}/ODDS_RATIO_95CI.txt

# clean-up
rm QTL.txt QTL.bed nonQTL.txt nonQTL.bed *QTL*region.txt

######################### KRT
# resources for celltype
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/KRT_minimal/masteroutput_all_with_colnames.txt
bed=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/KRT_minimal/snps_near_expressed_genes.bed

# reQTL
QTL_type=KRT_reQTL
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
Rscript ${dir}/calculate_OddsRatio_and_95CI.R ${QTL_in_region} ${QTL_not_in_region} ${nonQTL_in_region} ${nonQTL_not_in_region} ${QTL_type} ${region_type} ${dir}/ODDS_RATIO_95CI.txt

# PBSeQTL
QTL_type=KRT_PBSeQTL
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
Rscript ${dir}/calculate_OddsRatio_and_95CI.R ${QTL_in_region} ${QTL_not_in_region} ${nonQTL_in_region} ${nonQTL_not_in_region} ${QTL_type} ${region_type} ${dir}/ODDS_RATIO_95CI.txt

# IFNeQTL
QTL_type=KRT_IFNeQTL
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
Rscript ${dir}/calculate_OddsRatio_and_95CI.R ${QTL_in_region} ${QTL_not_in_region} ${nonQTL_in_region} ${nonQTL_not_in_region} ${QTL_type} ${region_type} ${dir}/ODDS_RATIO_95CI.txt

# clean-up
rm QTL.txt QTL.bed nonQTL.txt nonQTL.bed *QTL*region.txt

######################### MEL
# resources for celltype
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all_with_colnames.txt
bed=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/snps_near_expressed_genes.bed

# reQTL
QTL_type=MEL_reQTL
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
Rscript ${dir}/calculate_OddsRatio_and_95CI.R ${QTL_in_region} ${QTL_not_in_region} ${nonQTL_in_region} ${nonQTL_not_in_region} ${QTL_type} ${region_type} ${dir}/ODDS_RATIO_95CI.txt

# PBSeQTL
QTL_type=MEL_PBSeQTL
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
Rscript ${dir}/calculate_OddsRatio_and_95CI.R ${QTL_in_region} ${QTL_not_in_region} ${nonQTL_in_region} ${nonQTL_not_in_region} ${QTL_type} ${region_type} ${dir}/ODDS_RATIO_95CI.txt

# IFNeQTL
QTL_type=MEL_IFNeQTL
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
Rscript ${dir}/calculate_OddsRatio_and_95CI.R ${QTL_in_region} ${QTL_not_in_region} ${nonQTL_in_region} ${nonQTL_not_in_region} ${QTL_type} ${region_type} ${dir}/ODDS_RATIO_95CI.txt

# clean-up
rm QTL.txt QTL.bed nonQTL.txt nonQTL.bed *QTL*region.txt
rm region.bed

