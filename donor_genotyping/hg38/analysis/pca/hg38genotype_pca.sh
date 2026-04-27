#!/bin/bash
# goal: perform PCA on list of SNPs in our gencove hg38 genotyping dta
# written by Crystal Shan: 10/2022

### 0. set-up working space
bsub -Is -q interactive -W 8:00 -n1 -R rusage[mem=450000] -R "span[hosts=1]" /bin/bash
module load condas/2018-05-11
source activate sshan_isoform
module load bcftools/1.9

### 1. load datasets
wd=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/analysis/pca
vcfF=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
ancestry=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/ancestry.txt
cd ${wd}

### 2. filter out correlated SNPs by LD pruning (exclude r2 >= 0.02)
### ended up with 129,082 SNPs
bcftools +prune /nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz --max-LD 0.2 --window 2000kb -Ov -o /nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/analysis/pca/pruned_vcf.vcf.gz

### 3. filter out SNPs with HWE<1e-30 (recommended by PLINK, to detect any serious genotyping errors)
# https://www.cog-genomics.org/plink/1.9/filter
# fetch northern european donors
cat ${ancestry} | grep Northern_and_Central_European | cut -f1 > temp.txt
paste temp.txt temp.txt > european_donors.txt; rm temp.txt
# fetch list of SNPs that pass hwe threshold of 1e-30 # 7,812,613
plink --vcf ${wd}/pruned_vcf.vcf.gz --snps-only --keep ${wd}/european_donors.txt --hwe 1e-50 midp --recode --out plink_hwe_filtered
cat plink_hwe_filtered.map | cut -f2 > snps_pass_hwe.txt # all 129,082 SNPs pass hwe filter

### 4. filter out SNPs if: [1] not autosomal; [2] have missing genotype calls in more than 10% of the samples; [3] MAF<0.05;
# log file is in: /nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/analysis/pca/plink_for_pca.log
# 19591 variants removed due to missing genotype data
# 21909 variants removed due to minor allele threshold
# 84596 variants and 37 people pass filters and QC.
plink --vcf ${wd}/pruned_vcf.vcf.gz \
      --autosome \
      --biallelic-only \
      --extract ${wd}/snps_pass_hwe.txt \
      --geno 0.1 \
      --maf 0.05 \
      --recode vcf \
      --out plink_for_pca
bcftools view ${wd}/plink_for_pca.vcf -Oz -o ${wd}/plink_for_pca.vcf.gz

### 5. PCA usig the 84,956 variants
# compute genetic distances, contains a square, symmetric matrix of the IBS distances for all pairs of individuals. These values range, in theory, from 0 to 1. In practice, one would never expect to observe values near 0 -- even completely unrelated individuals would be expected to share a very large proportion of the genome identical by state by chance alone (i.e. as opposed to identity by descent). A value of 1 would indicate a MZ twin pair, or a sample duplication
plink --vcf ${wd}/plink_for_pca.vcf \
      --distance-matrix \
      --out dataForPCA
# PCA and PCA scatter plot is in ${wd}/hg38genotype_pca.R
