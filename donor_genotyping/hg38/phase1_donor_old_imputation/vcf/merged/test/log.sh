#!/bin/bash
#BSUB -n 6
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=30000]
#BSUB -q long
#BSUB -W 121:00
### script to organize vcf data quality metrics

### crystal shan 07/2021
### bsub -Is -q interactive -W 8:00 -n1 -R rusage[mem=450000] -R "span[hosts=1]" /bin/bash

############################
#### download gencove data in hg38
module load condas/2018-05-11
source activate sshan_isoform
module load bcftools/1.9
module load vcftools/0.1.16
#dir=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged
#f=${dir}/bd841628-fcc2-487a-8460-f5428237f0c9.merged.filtered.2.AFtagged.vcf.gz
#cd ${dir}
############################
## 1. separate vcf into snps and indels
#vcftools --gzvcf ${f} --keep-only-indels --recode --recode-INFO-all --out all.indels
#bcftools view all.indels.recode.vcf -Oz -o all.indels.vcf.gz
#vcftools --gzvcf ${f} --remove-indels --recode --recode-INFO-all --out all.snps
#bcftools view all.snps.recode.vcf -Oz -o all.snps.vcf.gz
## 2. allele frequency histogram
#bcftools query -f '%AF\n' all.snps.vcf.gz > all.snps.AF.txt
#bcftools query -f '%AF\n' all.indels.vcf.gz > all.indels.AF.txt 
## 3. indel histogram
## 4. vcf count overview
## 5. coverage depth
## plotted locally, check skin_eQTL dropbox
## 6. substitution types in percentage
#bcftools query -f '%REF\t%ALT\t[\t%GT]]\n' all.snps.vcf.gz > all.snps.genotypes.txt
#bcftools view -h all.snps.vcf.gz | tail -1 | cut -f 10- | tr '\t' '\n' > all.snps.vcf.colnames
## analyze in R
# 7. Het/HomoAlt ratio
dir=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/organized/
cd ${dir}
#F31 F41 F22 F23 F24 F25 F27 F28 F30 
#F32 F33 F34 F35 F36 F37 F38 F39 F40 
#F42 F44 F45 F46 F48 F49 F50 F51 F52 
#F53 F55 F56 F57 F58 F59 F60 F61 F62 F63 
#VB126 VB150 VB172 VB173 CB032 CB043 CB045
for donor in VB126 VB150 VB172 VB173 CB032 CB043 CB045;do
  f=${donor}.vcf.gz
  echo "working on "${f}
  date +"[%d-%m-%y] %T"
  bcftools view --threads 5 ${f} -f .,PASS -Oz -o this.${donor}.vcf.gz 
  bcftools +counts this.${donor}.vcf.gz --threads 4 > counts.${donor}.txt
  vcftools --gzvcf this.${donor}.vcf.gz --remove-indels --recode --recode-INFO-all --out this.${donor}.snps
  rm this.${donor}.vcf.gz 
  bcftools query -f '%REF%ALT\t[\t%GT]]\n' this.${donor}.snps.recode.vcf > genotype.${donor}.txt
  rm this.${donor}.snps*
  Ts=$(cat genotype.${donor}.txt | grep -w 'AG\|GA\|CT\|TC' | wc -l)
  Tv=$(cat genotype.${donor}.txt | grep -w 'AC\|CA\|GT\|TG' | wc -l)
  HomoRef=$(cat genotype.${donor}.txt | grep -w '0/0' | wc -l)
  Het=$(cat genotype.${donor}.txt | grep -w '0/1' | wc -l)
  HomoAlt=$(cat genotype.${donor}.txt | grep -w '1/1' | wc -l)
  echo ${donor} ${Ts} ${Tv} ${HomoRef} ${Het} ${HomoAlt} > summary.${donor}.txt
  rm genotypes.${donor}.txt
  echo "done with "${f}
  date +"[%d-%m-%y] %T"
done




