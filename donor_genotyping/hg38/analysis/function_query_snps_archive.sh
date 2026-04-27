#!/bin/bash
module load bcftools/1.9
dir=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38
donor=$1

cd ${dir}/analysis
f=${dir}/vcf/organized/${donor}.vcf.gz
#snp=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/analysis/vitiligo-common-snps.txt
snp=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/analysis/skin_disease_GWAS_SNPs.txt
echo "step 1: extract all snps of interest" ; date;
#bcftools view --include ID==@vitiligo-common-snps.txt ${f} -Oz -o ${dir}/analysis/${donor}.snp.vcf.gz
bcftools view --include ID==@skin_disease_GWAS_SNPs.txt ${f} -Oz -o ${dir}/analysis/${donor}.snp.vcf.gz

echo "step 2: annotate snp genotype calling based on observed/imputed. if imputed, annotated passed or lowconf."; date;
bcftools query -f "%ID\t%REF\t%ALT\t%FILTER[\t%GP\t%GT\t%DS]\n" ${dir}/analysis/${donor}.snp.vcf.gz > ${dir}/analysis/${donor}.snp.temp.txt
cat ${dir}/analysis/${donor}.snp.temp.txt | awk '{OFS="\t"}{if ( $5~/1,1e-10,1e-10/ || $5~/1e-10,1,1e-10/ || $5~/1e-10,1e-10,1/ ) print $0,"observed"; else if ($4=="PASS") print $0,"imputed_passed"; else print $0,"imputed_failed" }' > ${dir}/analysis/${donor}.snp.txt

rm ${dir}/analysis/${donor}.snp.vcf.gz ${dir}/analysis/${donor}.snp.temp.txt
