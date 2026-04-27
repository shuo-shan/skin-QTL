# script written by Crystal SHan 04/2022
# goal: filter donor vcf file to contain only SNPs that fall in ATAC-seq peaks that are within 500kbp of any gene
# generated file: 37donors_gene_proximal_open_region_SNPs.vcf.gz
# generated file: snps_in_gene_proximal_ATAC_merged_peaks.vcf.gz

# In this folder:
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/GWAS_SNPs/open_region_SNPs
module load condas/2018-05-11
source activate fastQTL
module load bedtools/2.29.2
module load bcftools/1.9  
genotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
Ensembl_gene_anno=/nl/umw_manuel_garber/human/skin/eQTLs/literature/UCSC_tracks/Ensembl_GRCh38.105_biotype_annotation.txt

##### KRT
celltype=KRT
openRegions=/nl/umw_manuel_garber/human/skin/eQTLs/ATACseq_DolphinNext/softlinks/narrowPeak/ATAC_${celltype}_merged_peaks.bed
# 1. find all snps that fall in ATAC-seq peaks --> 227,212 snps
bcftools view --regions-file ${openRegions} ${genotypeF} -Oz -o snps_in_ATAC_merged_peaks.vcf.gz
bcftools index snps_in_ATAC_merged_peaks.vcf.gz
# 2. find all open region snps that are proximal to any gene --> 226,806 snps 
cat ${Ensembl_gene_anno} | awk '{OFS="\t"}{print $1,$2-500000,$2+500000,$5}' | awk '{OFS="\t"}{if ($2~/-/) print $1,0,$3,$4; else print $0}' > temp_expanded_gene.bed
bedtools sort -i temp_expanded_gene.bed > temp_sorted_expanded_gene.bed
bedtools merge -i temp_sorted_expanded_gene.bed > temp_merged_expanded_gene.bed
bcftools view --regions-file temp_merged_expanded_gene.bed snps_in_ATAC_merged_peaks.vcf.gz -Oz -o snps_in_gene_proximal_${celltype}_ATAC_merged_peaks.vcf.gz
bcftools index snps_in_gene_proximal_${celltype}_ATAC_merged_peaks.vcf.gz
rm temp_expanded_gene.bed temp_sorted_expanded_gene.bed temp_merged_expanded_gene.bed
rm snps_in_ATAC_merged_peaks.vcf.gz snps_in_ATAC_merged_peaks.vcf.gz.csi

##### FRB
celltype=FRB
openRegions=/nl/umw_manuel_garber/human/skin/eQTLs/ATACseq_DolphinNext/softlinks/narrowPeak/ATAC_${celltype}_merged_peaks.bed
# 1. find all snps that fall in ATAC-seq peaks --> 227,212 snps
bcftools view --regions-file ${openRegions} ${genotypeF} -Oz -o snps_in_ATAC_merged_peaks.vcf.gz
wait
bcftools index snps_in_ATAC_merged_peaks.vcf.gz
# 2. find all open region snps that are proximal to any gene --> 226,806 snps 
cat ${Ensembl_gene_anno} | awk '{OFS="\t"}{print $1,$2-500000,$2+500000,$5}' | awk '{OFS="\t"}{if ($2~/-/) print $1,0,$3,$4; else print $0}' > temp_expanded_gene.bed
bedtools sort -i temp_expanded_gene.bed > temp_sorted_expanded_gene.bed
bedtools merge -i temp_sorted_expanded_gene.bed > temp_merged_expanded_gene.bed
bcftools view --regions-file temp_merged_expanded_gene.bed snps_in_ATAC_merged_peaks.vcf.gz -Oz -o snps_in_gene_proximal_${celltype}_ATAC_merged_peaks.vcf.gz
wait
bcftools index snps_in_gene_proximal_${celltype}_ATAC_merged_peaks.vcf.gz
rm temp_expanded_gene.bed temp_sorted_expanded_gene.bed temp_merged_expanded_gene.bed
rm snps_in_ATAC_merged_peaks.vcf.gz snps_in_ATAC_merged_peaks.vcf.gz.csi

###
celltype=MEL
genotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
# 1. find all snps that fall in ATAC-seq peaks --> 227,212 snps
bcftools view --regions-file ${openRegions} ${genotypeF} -Oz -o snps_in_ATAC_merged_peaks.vcf.gz
wait
bcftools index snps_in_ATAC_merged_peaks.vcf.gz
# 2. find all open region snps that are proximal to any gene --> 226,806 snps 
cat ${Ensembl_gene_anno} | awk '{OFS="\t"}{print $1,$2-500000,$2+500000,$5}' | awk '{OFS="\t"}{if ($2~/-/) print $1,0,$3,$4; else print $0}' > temp_expanded_gene.bed
bedtools sort -i temp_expanded_gene.bed > temp_sorted_expanded_gene.bed
bedtools merge -i temp_sorted_expanded_gene.bed > temp_merged_expanded_gene.bed
bcftools view --regions-file temp_merged_expanded_gene.bed snps_in_ATAC_merged_peaks.vcf.gz -Oz -o snps_in_gene_proximal_${celltype}_ATAC_merged_peaks.vcf.gz
bcftools index snps_in_gene_proximal_${celltype}_ATAC_merged_peaks.vcf.gz
rm temp_expanded_gene.bed temp_sorted_expanded_gene.bed temp_merged_expanded_gene.bed
rm snps_in_ATAC_merged_peaks.vcf.gz snps_in_ATAC_merged_peaks.vcf.gz.csi

