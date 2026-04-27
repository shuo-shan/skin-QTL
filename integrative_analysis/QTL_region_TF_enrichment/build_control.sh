#!/bin/bash
# written by shuo.shan@umassmed.edu, 03/2024

snp=$1
seed=$2

snp=rs55770741
seed=1
dir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/QTL_region_TF_enrichment
cd ${dir}

# for the given QTL, fetch the linked gene
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all_with_colnames.txt
grep -w ${snp} ${f} | awk '{OFS=FS="\t"}NR>1{if ($7 != "." && $7 < 0.000001) print $1,$2,$7}' > pairs.txt
cat pairs.txt | cut -f1 | sort | uniq > temp_snps.txt

# for the given QTL, fetch its MAF
module load bcftools
bcftools query -f'%ID AF=%AF\n' pruned_reQTL.vcf.gz | sed 's/AF=//g' | awk '{OFS="\t"}{print $1,1-$2}' > pruned_reQTL_MAF.txt
bcftools query 

# get the list of SNPs that fall within 1Mb window of the QTL-gene
genome=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
zcat /pi/manuel.garber-umw/human/skin/eQTLs/literature/UCSC_tracks/Homo_sapiens.GRCh38.105.gtf.gz | awk '{OFS=FS="\t"}{if ($3=="gene") print "chr"$0}' > temp_gene_model.txt
cat temp_gene_model.txt | grep gene_name | awk '{OFS=FS="\t"}{print $1,$4,$5,$7}' > temp_ver1_left.txt
cat temp_gene_model.txt | grep gene_name | sed 's/.*gene_name "//g' | sed 's/";.*//g' > temp_ver1_right.txt
paste temp_ver1_left.txt temp_ver1_right.txt > temp_ver1.txt
cat temp_ver1.txt | awk '{OFS=FS="\t"}{if ($4=="+") print $1,$2,$2+1,$4,$5}' > temp_ver2_plus.txt
cat temp_ver1.txt | awk '{OFS=FS="\t"}{if ($4=="-") print $1,$3-1,$3,$4,$5}' > temp_ver2_minus.txt
cat temp_ver2_plus.txt temp_ver2_minus.txt | bedtools sort -i stdin > tss.bed
cat ${genome} | cut -f1 > chr.txt
grep -w -f chr.txt tss.bed | bedtools slop -b 500000 -i stdin -g ${genome} > tss_window.bed
rm temp_gene_model.txt temp_ver1_left.txt temp_ver1_right.txt temp_ver1.txt temp_ver2_plus.txt temp_ver2_minus.txt tss.bed chr.txt


# prune genotype table to only include LD < 0.2
genotype_table=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
module load bcftools
bcftools +prune -m 0.2 -w 1000 ${genotype_table} -Oz -o pruned_genotype_mastertable.vcf.gz
bcftools index pruned_genotype_mastertable.vcf.gz


# for the given QTL, fetch the linked gene
gene=$(grep -w ${snp} pairs.txt | cut -f2)
## for the given QTL, fetch its MAF
MAF=$(grep -w ${snp} pruned_reQTL_MAF.txt | cut -f2) 
# get the list of SNPs that fall within 1Mb window of the QTL-gene
grep -w ${gene} tss_window.bed > region.bed
module load bcftools
vcf=pruned_genotype_mastertable.vcf.gz
bcftools view -R region.bed ${vcf} -Oz -o region.vcf.gz 
# get the list of SNPs that fall within 10% of the MAF of the QTL
bcftools +fill-tags region.vcf.gz -Oz -o annotated_region.vcf.gz -- -t AF
conda activate fastQTL
Rscript get_MAF_range ${MAF}
