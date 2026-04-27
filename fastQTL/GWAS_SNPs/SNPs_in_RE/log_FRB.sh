#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=100000]
#BSUB -q long
#BSUB -W 08:00
#BSUB -e "./%J%I.err"
#BSUB -o "./%J%I.out"

# script written by Crystal SHan 06/2022
# goal: filter donor vcf file to contain only SNPs that fall in ATAC-seq peaks that are within 500kbp of any gene
# generated file: 37donors_gene_proximal_regulatory_region_SNPs.vcf.gz


# In this folder:
dir=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/GWAS_SNPs/SNPs_in_RE
cd $dir
module load condas/2018-05-11
source activate fastQTL
module load bedtools/2.29.2
module load bcftools/1.9  
genotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
Ensembl_gene_anno=/nl/umw_manuel_garber/human/skin/eQTLs/literature/UCSC_tracks/Ensembl_GRCh38.105_biotype_annotation.txt
geneBody=/nl/umw_manuel_garber/human/skin/eQTLs/literature/UCSC_tracks/Ensembl_GRCh38.105_biotype_annotation.txt
geneExprs_krt=/nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq/keratinocytes/pipeline_01142022/analysis/highly_expressed_genes.txt
geneExprs_frb=/nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq/fibroblasts/pipeline_01152022/analysis/highly_expressed_genes.txt
geneExprs_mel=/nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq/melanocytes/pipeline_01122022/analysis/highly_expressed_genes.txt
promoters_frb=/nl/umw_manuel_garber/human/skin/eQTLs/chromatin/RegulatoryElements/FRB_promoters_annotated.txt
enhancers_frb=/nl/umw_manuel_garber/human/skin/eQTLs/chromatin/RegulatoryElements/FRB_enhancers.txt
promoters_krt=/nl/umw_manuel_garber/human/skin/eQTLs/chromatin/RegulatoryElements/KRT_promoters_annotated.txt
enhancers_krt=/nl/umw_manuel_garber/human/skin/eQTLs/chromatin/RegulatoryElements/KRT_enhancers.txt
promoters_mel=/nl/umw_manuel_garber/human/skin/eQTLs/chromatin/RegulatoryElements/MEL_promoters_annotated.txt
enhancers_mel=/nl/umw_manuel_garber/human/skin/eQTLs/chromatin/RegulatoryElements/MEL_enhancers.txt
expressedGenes=${geneExprs_frb}
promoters=${promoters_frb} # promoter name is column#4
enhancers=${enhancers_frb} # enhancer name is column#21
celltype=FRB # KRT, FRB, MEL

# generate BED file for gene body, promoter, and enhancers. Sort and merge.
cat ${geneBody} | awk '{OFS="\t"}{print $1,$2,$3,"gene_"$5"_"$4"_"$1"_"$2"_"$3}' > ${celltype}_gene_body_expressed_genes.bed #append gene name and gene Ensembl ID by "_"
cat ${promoters} | awk '{OFS="\t"}{print $1,$2,$3,$4}' > ${celltype}_promoters.bed
cat ${enhancers} | awk '{OFS="\t"}{print $1,$2,$3,$21}' > ${celltype}_enhancers.bed
module load bedtools/2.29.2
cat ${celltype}_gene_body_expressed_genes.bed ${celltype}_promoters.bed ${celltype}_enhancers.bed | bedtools sort -i stdin | bedtools merge -i stdin -c 4 -o collapse > ${celltype}_regulatory_regions.bed

regRegionsBED=${dir}/${celltype}_regulatory_regions.bed
cat ${regRegionsBED} | cut -f1,2,3,4 > ${celltype}.regions.bed
regRegions=${dir}/${celltype}.regions.bed
# 1. find all snps that fall in gene body or promoter or enhancers --> __ snps
bcftools view --regions-file ${regRegions} ${genotypeF} -Oz -o ${celltype}_snps_in_regions.vcf.gz
bcftools index ${celltype}_snps_in_regions.vcf.gz
# ^ DON"T DELETE. Time consuming step!

# 2. find all regulatory region snps that are proximal to any gene --> __ snps 
cat ${Ensembl_gene_anno} | awk '{OFS="\t"}{print $1,$2-500000,$2+500000,$5}' | awk '{OFS="\t"}{if ($2~/-/) print $1,0,$3,$4; else print $0}' > ${celltype}_temp_expanded_gene.bed
bedtools sort -i ${celltype}_temp_expanded_gene.bed > ${celltype}_temp_sorted_expanded_gene.bed
bedtools merge -i ${celltype}_temp_sorted_expanded_gene.bed > ${celltype}_temp_merged_expanded_gene.bed
bcftools view --regions-file ${celltype}_temp_merged_expanded_gene.bed ${celltype}_snps_in_regions.vcf.gz -Oz -o snps_in_gene_proximal_${celltype}_in_regulatory_region.vcf.gz
bcftools index snps_in_gene_proximal_${celltype}_in_regulatory_region.vcf.gz
#rm ${celltype}_gene_body_expressed_genes.txt ${celltype}_gene_body_expressed_genes.bed ${celltype}_promoters.bed ${celltype}_enhancers.bed
#rm temp_expanded_gene.bed temp_sorted_expanded_gene.bed temp_merged_expanded_gene.bed

