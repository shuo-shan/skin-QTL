#!/bin/bash

reQTL=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/QTL_region_TF_enrichment/reQTL.bed
outDir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/QTL_region_TF_enrichment/

##########################################################################################
# CHUNK 1
# overlap TF CHIPseq peak with 1E-5 thresholded reQTLs
cd /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENCODE/humanTFChIP
module load bedtools
for f in *.bed.gz;do TF=$(echo ${f} | cut -d'_' -f2 | sed 's/-human//g'); ID=$(echo ${f} | cut -d'_' -f1); bedtools intersect -a ${reQTL} -b ${f} | cut -f1-6 | awk -v TF=${TF} -v ID=${ID} '{OFS=FS="\t"}{print $0,TF,ID}' >> MELreQTL1E-06_overlapping_ChIPpeaks.txt; echo ${f}; done
mv MELreQTL1E-06_overlapping_ChIPpeaks.txt ${outDir}
# ran intersect_reQTL_with_TFChIPseqPeaks.R to calculate whether any pair of TF always show-up together. found some interesting results. see word doc "QTL-region-TF-enrichment". 

##########################################################################################
# CHUNK 2
# expand each CHIPseq peak file by 500 up and down stream, overlap with 1E-5 thresholded reQTLs
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all_with_colnames.txt
cat ${f} | awk 'NR>1{if ($7 != "." && $7 < 0.00001) print $1}' | sort | uniq > reQTL.txt
genome=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
cd /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENCODE/humanTFChIP
module load bedtools
for f in *.bed.gz;do
	TF=$(echo ${f} | cut -d'_' -f2 | sed 's/-human//g')
	ID=$(echo ${f} | cut -d'_' -f1)
	bedtools slop -b 500 -i ${f} -g ${genome} > temp_expanded_${f}.bed
	bedtools intersect -a ${reQTL} -b temp_expanded_${f}.bed | cut -f1-6 | awk -v TF=${TF} -v ID=${ID} '{OFS=FS="\t"}{print $0,TF,ID}' >> MELreQTL1E-05_overlapping_expandedChIPpeaks.txt
	echo ${f}
	rm temp_expanded_${f}.bed
done

mv MELreQTL1E-05_overlapping_expandedChIPpeaks.txt ${outDir}
cd ${outDir}

##########################################################################################
# CHUNK 3
# overlap with 1E-6 thresholded reQTLs and prune the overlapped ones
outDir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/QTL_region_TF_enrichment/
cd ${outDir}
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all_with_colnames.txt
cat ${f} | awk 'NR>1{if ($7 != "." && $7 < 0.000001) print $1}' | sort | uniq > reQTL.txt
bed=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/master_filtered_genotype.bed
grep -w -f reQTL.txt ${bed} > reQTL.bed
reQTL=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/QTL_region_TF_enrichment/reQTL.bed
genome=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
module load bedtools

cd /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENCODE/humanTFChIP
rm MELreQTL1E-6_overlapping_ChIPpeaks.txt
for f in *.bed.gz;do
	TF=$(echo ${f} | cut -d'_' -f2 | sed 's/-human//g')
	ID=$(echo ${f} | cut -d'_' -f1)
	bedtools intersect -a ${reQTL} -b ${f} | cut -f1-6 | awk -v TF=${TF} -v ID=${ID} '{OFS=FS="\t"}{print $0,TF,ID}' >> MELreQTL1E-6_overlapping_ChIPpeaks.txt
	echo ${f}
done
mv MELreQTL1E-6_overlapping_ChIPpeaks.txt ${outDir}
cd ${outDir}

genotype_table=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
module load bcftools
cat MELreQTL1E-6_overlapping_ChIPpeaks.txt | cut -f4 | sort | uniq > MELreQTL1E-6_overlapping_ChIPpeaks_SNPID.txt
bcftools view --include ID==@MELreQTL1E-6_overlapping_ChIPpeaks_SNPID.txt ${genotype_table} -Oz -o MELreQTL1E-6_overlapping_ChIPpeaks_SNPID.vcf.gz
bcftools +prune -m 0.2 -w 1000  MELreQTL1E-6_overlapping_ChIPpeaks_SNPID.vcf.gz -Oz -o pruned_MELreQTL1E-6_overlapping_ChIPpeaks_SNPID.vcf.gz
bcftools index pruned_MELreQTL1E-6_overlapping_ChIPpeaks_SNPID.vcf.gz
bcftools view -H pruned_MELreQTL1E-6_overlapping_ChIPpeaks_SNPID.vcf.gz | cut -f3 > pruned_MELreQTL1E-6_overlapping_ChIPpeaks_SNPID.txt

grep -w -f pruned_MELreQTL1E-6_overlapping_ChIPpeaks_SNPID.txt MELreQTL1E-6_overlapping_ChIPpeaks.txt > pruned_MELreQTL1E-6_overlapping_ChIPpeaks.txt
rm pruned_MELreQTL1E-6_overlapping_ChIPpeaks_SNPID.txt pruned_MELreQTL1E-6_overlapping_ChIPpeaks_SNPID.vcf.gz

# run intersect_reQTL_with_TFChIPseqPeaks.R to calculate whether any pair of TF always show-up together. found some interesting results. see word doc "QTL-region-TF-enrichment".



##########################################################################################
# CHUNK 4
# check which reQTLs overlap the TF of interest
TF=CTCF
cd ${outDir}
cat pruned_MELreQTL1E-6_overlapping_ChIPpeaks.txt | cut -f4 | sort | uniq > temp_reQTL.txt
grep -w -f temp_reQTL.txt reQTL.bed > temp_reQTL.bed
cd /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENCODE/humanTFChIP/
less filtered_metadata.txt | grep ${TF} | cut -d';' -f3 > temp_id.txt
ls | grep -f temp_id.txt > temp_files.txt
while read f;do
	fname=$(echo ${f} | cut -d'_' -f1)
	bedtools intersect -a ${outDir}/temp_reQTL.bed -b ${f} | awk -v ID=${fname} '{OFS=FS="\t"}{print $0,ID}' >> temp_reQTL_TF-of-interest_overlap.bed
	echo ${f}
done < temp_files.txt
cat temp_reQTL_TF-of-interest_overlap.bed | cut -f4,44  | sort -k1 # this prints the SNPs that overlap with which file.




#################
# CHUNK 5
# overlap with 1E-5 thresholded reQTLs and prune the overlapped ones
# 1. fetch 1E-5 thresholded MEL reQTLs:
outDir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/QTL_region_TF_enrichment/
cd ${outDir}
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all_with_colnames.txt
cat ${f} | awk 'NR>1{if ($7 != "." && $7 < 0.00001) print $1}' | sort | uniq > reQTL.txt
bed=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/master_filtered_genotype.bed
grep -w -f reQTL.txt ${bed} > reQTL.bed
# 2. fetch reQTLs that overlap ChIPseq peaks
reQTL=${outDir}/reQTL.bed
genome=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
module load bedtools
cd /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENCODE/humanTFChIP
rm MELreQTL1E-05_overlapping_ChIPpeaks.bed
for f in *.bed.gz;do
        TF=$(echo ${f} | cut -d'_' -f2 | sed 's/-human//g')
        ID=$(echo ${f} | cut -d'_' -f1)
        bedtools intersect -a ${reQTL} -b ${f} | cut -f1-6 | awk -v TF=${TF} -v ID=${ID} '{OFS=FS="\t"}{print $0,TF,ID}' >> MELreQTL1E-05_overlapping_ChIPpeaks.bed
        echo ${f}
done
mv MELreQTL1E-05_overlapping_ChIPpeaks.bed ${outDir}
cd ${outDir}
# 3. check which reQTLs overlap the TF of interest
TF=MAX
cd ${outDir}
cat MELreQTL1E-05_overlapping_ChIPpeaks.bed | cut -f4 | sort | uniq > temp_reQTL.txt
grep -w -f temp_reQTL.txt reQTL.bed > temp_reQTL.bed
cd /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENCODE/humanTFChIP/
less filtered_metadata.txt | grep ${TF} | cut -d';' -f3 > temp_id.txt
ls | grep -f temp_id.txt > temp_files.txt
rm temp_reQTL_${TF}_overlap.bed
while read f;do
        fname=$(echo ${f} | cut -d'_' -f1)
        bedtools intersect -a ${outDir}/temp_reQTL.bed -b ${f} | awk -v ID=${fname} '{OFS=FS="\t"}{print $0,ID}' >> temp_reQTL_${TF}_overlap.bed
        echo ${f}
done < temp_files.txt
cat temp_reQTL_${TF}_overlap.bed | cut -f4 | sort | uniq > ${outDir}/reQTL_${TF}_overlap.txt
rm temp_id.txt temp_files.txt temp_reQTL_${TF}_overlap.bed
# do this for MAX and MYC. then MAX:MYC
cd ${outDir}
comm -12 reQTL_MAX_overlap.txt reQTL_MYC_overlap.txt > reQTL_MAX-MYC_overlap.txt
# 4. get the linked gene for each reQTL
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all_with_colnames.txt
TF=MYC
query=reQTL_${TF}_overlap.txt
cat ${f} | awk 'NR>1{OFS="\t"}{if ($7 != "." && $7 < 0.00001) print $1,$2,$7,$9}' | grep -w -f ${query} | cut -f2 | sort | uniq > reQTL_genes_${TF}_overlap.txt
# check the best PBS eQTL
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all_with_colnames.txt
gene=HSPA12A
cat ${f} | awk 'NR>1{OFS="\t"}{if ($15 != "." && $15 < 0.00001) print $1,$2,$15,$17}' | grep -w ${gene} | sort -gk3
# check the best reQTL
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all_with_colnames.txt
gene=SMIM3
cat ${f} | awk 'NR>1{OFS="\t"}{if ($7 != "." && $7 < 0.00001) print $1,$2,$7,$9}' | grep -w ${gene} | sort -gk3
# check the best reQTL that falls in the TF peak
cd ${outDir}
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all_with_colnames.txt
TF=MAX
gene=SMIM3
query=reQTL_${TF}_overlap.txt
cat ${f} | awk 'NR>1{OFS="\t"}{if ($7 != "." && $7 < 0.00001) print $1,$2,$7,$9}' | grep -w -f ${query} | grep -w ${gene}





#################
##### CHUNK 6
### check which reQTLs overlap the TF of interest
TF=CTCF
cd ${outDir}
cat MELreQTL1E-05_overlapping_ChIPpeaks.bed | cut -f4 | sort | uniq > temp_reQTL.txt
grep -w -f temp_reQTL.txt reQTL.bed > temp_reQTL.bed
cd /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENCODE/humanTFChIP/
less filtered_metadata.txt | grep ${TF} | cut -d';' -f3 > temp_id.txt
ls | grep -f temp_id.txt > temp_files.txt
while read f;do cp ${f} ../IGV; echo ${f}; done < temp_files.txt
rm temp_reQTL_${TF}_overlap.bed
while read f;do
        fname=$(echo ${f} | cut -d'_' -f1)
        bedtools intersect -a ${outDir}/temp_reQTL.bed -b ${f} | awk -v ID=${fname} '{OFS=FS="\t"}{print $0,ID}' >> temp_reQTL_${TF}_overlap.bed
        echo ${f}
done < temp_files.txt
cat temp_reQTL_${TF}_overlap.bed | cut -f4 | sort | uniq > ${outDir}/reQTL_${TF}_overlap.txt
rm temp_id.txt temp_files.txt temp_reQTL_${TF}_overlap.bed
cd ${outDir}
### get the linked gene for each reQTL
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all_with_colnames.txt
query=reQTL_${TF}_overlap.txt
cat ${f} | awk 'NR>1{OFS="\t"}{if ($7 != "." && $7 < 0.00001) print $1,$2,$7,$9}' | grep -w -f ${query} | sort -gk3 | cut -f2 | awk '!seen[$1]++' > reQTL_genes_${TF}_overlap.txt
# check the best PBS eQTL
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all_with_colnames.txt
gene=HSPA12A
cat ${f} | awk 'NR>1{OFS="\t"}{if ($15 != "." && $15 < 0.00001) print $1,$2,$15,$17}' | grep -w ${gene} | sort -gk3
# check the best reQTL
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all_with_colnames.txt
gene=LINC00920
cat ${f} | awk 'NR>1{OFS="\t"}{if ($7 != "." && $7 < 0.00001) print $1,$2,$7,$9}' | grep -w ${gene} | sort -gk3
# check the best reQTL that falls in the TF peak
cd ${outDir}
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all_with_colnames.txt
TF=CTCF
gene=AMN1
query=reQTL_${TF}_overlap.txt
cat ${f} | awk 'NR>1{OFS="\t"}{if ($7 != "." && $7 < 0.00001) print $1,$2,$7,$9}' | grep -w -f ${query} | grep -w ${gene}


#################
##### CHUNK 7
### overlap ALL TFs
cd ${outDir}
cat MELreQTL1E-05_overlapping_ChIPpeaks.bed | cut -f4,7 | sort -k1 | awk '!seen[$1"_"$2]++'> temp1.txt
awk '{if ($1 in snp_to_tf) snp_to_tf[$1] = snp_to_tf[$1] "," $2; else snp_to_tf[$1] = $2;} END {for (snp in snp_to_tf) print snp, snp_to_tf[snp];}' temp1.txt | tr ' ' '\t' > temp2.txt
cat temp2.txt | awk '{FS=OFS="\t"}{n = split($2, a, ","); print $0, n}' | sort -r -nk3 > MELreQTL1E-05_overlapping_TFpeaks.txt
rm temp1.txt temp2.txt
# more TF overlapping = most interesting SNPs? let's look at what genes they map to and modeling stats
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all_with_colnames.txt
cat MELreQTL1E-05_overlapping_TFpeaks.txt | cut -f1 | sort | uniq > temp1.txt
cat ${f} | awk 'NR>1{OFS="\t"}{if ($7 != "." && $7 < 0.00001) print $1,$2,$7,$9}' | grep -w -f temp1.txt > temp2.txt
echo -e "SNP\tGENE\tpval\tbeta\tTF\tNO_TF" > MELreQTL1E-05_overlapping_TFpeaks_with_modelingresults.txt 
join -1 1 -2 1 <(sort -k1 temp2.txt) <(sort -k1 MELreQTL1E-05_overlapping_TFpeaks.txt) | sort -r -nk6 >> MELreQTL1E-05_overlapping_TFpeaks_with_modelingresults.txt 
rm temp1.txt temp2.txt
